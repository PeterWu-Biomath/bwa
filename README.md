# BWA-Capture: Targeted Alignment with Capture Sub-Index

BWA-Capture extends [BWA](README_bwa.md) with a capture sub-index that accelerates alignment for targeted sequencing (e.g., exome, panel) by starting search from regions of interest.

## Overview

In panel sequencing, only a small fraction of the genome is targeted for enrichment and sequencing. A typical WES panel has ~50M target length and diagnostic panel would span at 100K~2M. BWA-Capture exploits this sparsity: reads that originate from target regions are enriched through a dense sub-index built exclusively over those regions, while off-target reads fall back to the full genome BWT with minimal overhead.

BWA-Capture is designed to produce results identical to standard BWA, so it can replace BWA in existing pipelines without downstream validation.

To achieve this, BWA-Capture introduces three modules. The first is a dense suffix array over the target regions with an uncompressed OCC table — trading index size for fast rank queries during SMEM search. When SMEM extension reaches the boundary of a target region or fails to extend, the sub-SA interval is translated to a genome SA interval via the position table, and alignment continues seamlessly on the full genome BWT. A lookup table maps (sub-SA rank, query length) to genome SA rank, compressed via run-length encoding (RLE).

The second module is a cached genome SA lookup — a table of (SA rank, genome position) entries covering only positions that fall within regions of interest. On a cache hit, the position is returned directly; on a miss, the query falls back to the standard SA lookup.

The third module introduces a shortcut in `capt_mem_align1_core` for reads that align perfectly within a target region. During index construction, the full `mem_alnreg_t` result is precomputed for every suffix in the sub-SA across a configurable read-length range `[rlen_lo, rlen_hi]`. At runtime, when a read's SMEM covers its full length and its length falls within that range, the cached result is returned via `capt_perfect_get` — skipping chaining, Smith-Waterman, and deduplication entirely.

Three modules work together:

```
Read
 │
 ├─ 1. SA for target regions ──→ sub-index SMEM search
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
bwa mem_capture <capt_prefix> [mem options] <ref.fa> <in.fq>
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
