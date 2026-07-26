# BWA-Capture: Targeted Alignment with Capture Sub-Index

BWA-Capture extends [BWA](README_bwa.md) with a capture sub-index that accelerates alignment for targeted sequencing (e.g., exome, panel) by restricting the search space to regions of interest.

## Overview

Three modules work together:

```
Read
 │
 ├─ 1. Dense SA for target regions ──→ sub-index SMEM search
 │      (capt_smem1, capt_extend)         │
 │                                         ├─ hits target → fast path
 │                                         └─ no hit ──→ fallback to genome BWT
 │
 ├─ 2. Cached genome SA lookup ──→ position table
 │      (capt_lookup_pos)              │
 │       converts sub-SA rank → genome SA rank in O(1)
 │
 └─ 3. Perfect-match shortcut ──→ precomputed alignment table
        (capt_perfect_get)            │
        read_len ∈ [rlen_lo, rlen_hi]  │
        + SMEM covers [0, len)         │
        → return cached mem_alnreg_t directly
        → skip chaining, SW, dedup entirely
```

## Module 1: Dense SA for Target Regions

Builds a standalone suffix array over BED target regions (plus padding), with dense BWT, OCC, and LCP. All suffixes are sorted by lexicographic order of the concatenated region sequences (forward + reverse complement).

- **SMEM search**: `capt_smem1` runs entirely on the sub-index for reads that hit the target
- **Extension**: `capt_extend` / `capt_extend_bail` for backward/forward extension
- **Fallback**: reads that escape the target region fall back to the genome BWT seamlessly via interval translation

## Module 2: Cached Genome SA Lookup

The position table maps sub-SA ranks to genome PAC coordinates in O(1). Built via multi-string BWT merge:

- Sub-SA suffixes and genome BWT suffixes are merged in a single pass
- Each sub-SA entry records its corresponding genome SA rank
- At runtime: `capt_lookup_pos(capt, rank)` returns the genome position directly, no binary search

## Module 3: Perfect-Match Shortcut

For reads in a precomputed length range that exactly match a suffix in the sub-SA, the full alignment result (`mem_alnreg_t`) is cached in the index.

**Index time** (Step 9): For each SA suffix with ≥ `rlen_lo` bases, for each read length in `[rlen_lo, rlen_hi]`, run `capt_mem_align1_core` and store all resulting `mem_alnreg_t` entries.

**Runtime**: In `mem_collect_intv`, the first `capt_smem1` call checks if the SMEM covers `[0, len)`. If so and `len ∈ [rlen_lo, rlen_hi]`, return the SA index immediately. `capt_mem_align1_core` then looks up `capt_perfect_get(capt, sa_idx, len, &count)` and returns the cached results — skipping `mem_chain_flt`, `mem_flt_chained_seeds`, `mem_chain2aln`, and `mem_sort_dedup_patch` entirely.

## Building the Capture Index

```bash
# Build genome BWT (standard BWA)
bwa index ref.fa

# Build capture sub-index
bwa index-capture ref.fa targets.bed -p 200 -r 100,150 -o prefix

# Options:
#   -p INT     padding around each region [0]
#   -r LO,HI   read length range for perfect-match precomputation [100,150]
#   -k INT     k-mer size, > max read length [200]
#   -o STR     output prefix (writes prefix.capt)
#   -d STR     dump merged regions to directory
```

This produces `prefix.capt` containing:
- Dense SA + BWT + OCC + LCP
- RLE-encoded per-suffix metadata
- Position table (sub-SA rank → genome PAC)
- Perfect-match alignment table (indexed by SA_idx × read_len)

## Alignment

```bash
bwa mem -C prefix ref.fa reads.fq > aln.sam
```

The `-C` flag loads the capture index and activates all three modules.

## Files

| File | Role |
|------|------|
| `capt.h` | Data structures, constants, lookup API |
| `capt.c` | Runtime: load, query, verify, extend |
| `bwa_index_capt.c` | Index construction (Steps 1-9) |
| `bwamem_capt.c` | Alignment core with perfect-match shortcut |
| `bwamem.c` | Modified `mem_collect_intv` for shortcut detection |

## License

Same as BWA (GPLv3). See [README_bwa.md](README_bwa.md) for the original BWA documentation.
