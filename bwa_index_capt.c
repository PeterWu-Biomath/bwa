/*
 * bwa_index_capt.c — Build capture sub-index for accelerated bwa-mem.
 *
 * Usage: bwa index-capture <ref.fa> <target.bed> [-p padding] [-k 200] -o <capt_name>
 *
 * Construction pipeline:
 *   1. Parse BED → map chr names to PAC offsets, merge overlaps
 *   2. Extract padded region sequences from genome PAC (no concat)
 *   3. Build dense SA by direct k-mer extraction + integer sort
 *   4. Build dense BWT + OCC in one pass over SA
 *   5. Build LCP (Kasai)
 *   6. Load genome BWT → compute tail-RLE per suffix
 *   7. Compute position table via genome bwt_sa()
 *   8. Serialize .capt file
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>
#include "bwa.h"
#include "bwt.h"
#include "bntseq.h"
#include "utils.h"
#include "kvec.h"
#include "kseq.h"
#include "bwamem.h"
#include "capt.h"

extern int bwa_verbose;

/* ── Step 1: BED parsing ───────────────────────────────────────────── */

/*
 * Find chromosome index in bns by name.
 * Returns -1 if not found.
 */
static int
capt_bns_name2id(const bntseq_t *bns, const char *name)
{
    for (int i = 0; i < bns->n_seqs; i++) {
        if (strcmp(bns->anns[i].name, name) == 0) return i;
    }
    return -1;
}

/*
 * Parse a BED file into a list of regions.
 *
 * For each BED line, maps chromosome name to its PAC offset + length
 * via bns.  Stores chrom_offset, chrom_len, start, end.
 *
 * Returns a capt_regions_v, caller frees.
 */
static capt_regions_v *
capt_bed_read(const char *bed_path, const bntseq_t *bns)
{
    capt_regions_v *regs = calloc(1, sizeof(*regs));
    kv_init(*regs);

    FILE *fp = fopen(bed_path, "r");
    if (!fp) {
        fprintf(stderr, "[index-capture] ERROR: cannot open BED '%s'\n", bed_path);
        return regs;
    }

    char line[4096];
    while (fgets(line, sizeof(line), fp)) {
        /* skip empty, comments, track lines */
        if (line[0] == '\n' || line[0] == '#' || strncmp(line, "track", 5) == 0)
            continue;

        char chrom[256];
        int64_t start, end;
        if (sscanf(line, "%255s %ld %ld", chrom, &start, &end) < 3)
            continue;

        int chr_id = capt_bns_name2id(bns, chrom);
        if (chr_id < 0) {
            fprintf(stderr, "[index-capture] WARNING: chromosome '%s' not in reference, skipping\n", chrom);
            continue;
        }

        capt_region_t r;
        r.chrom_offset = bns->anns[chr_id].offset;
        r.chrom_len    = bns->anns[chr_id].len;
        r.start = start;
        r.end   = end;
        if (r.start > r.end) { int64_t tmp = r.start; r.start = r.end; r.end = tmp; }
        if (r.start < 0) r.start = 0;
        if (r.end > r.chrom_len) r.end = r.chrom_len;
        kv_push(capt_region_t, *regs, r);
    }
    fclose(fp);

    fprintf(stderr, "[index-capture]   parsed %zu regions from %s\n", regs->n, bed_path);
    return regs;
}


/* ── Step 2: Pad regions ───────────────────────────────────────────── */

static void
capt_bed_pad(capt_regions_v *regs, int padding)
{
    for (size_t i = 0; i < regs->n; i++) {
        capt_region_t *r = &regs->a[i];
        r->start -= padding;
        r->end   += padding;
        if (r->start < 0) r->start = 0;
        if (r->end > r->chrom_len) r->end = r->chrom_len;
    }
}


/* ── Step 3: Region merge ──────────────────────────────────────────── */

static int
capt_region_cmp(const void *a, const void *b)
{
    const capt_region_t *ra = a, *rb = b;
    if (ra->chrom_offset != rb->chrom_offset)
        return (ra->chrom_offset > rb->chrom_offset) ? 1 : -1;
    if (ra->start != rb->start)
        return (ra->start > rb->start) ? 1 : -1;
    return 0;
}

/*
 * Sort regions by (chrom_offset, start) and merge overlapping intervals
 * on the same chromosome.
 */
static void
capt_bed_merge(capt_regions_v *regs)
{
    if (regs->n == 0) return;

    /* Sort */
    qsort(regs->a, regs->n, sizeof(capt_region_t), capt_region_cmp);

    /* Single-pass merge */
    size_t j = 0;
    for (size_t i = 1; i < regs->n; i++) {
        capt_region_t *cur = &regs->a[j];
        capt_region_t *nxt = &regs->a[i];

        if (nxt->chrom_offset == cur->chrom_offset &&
            nxt->start <= cur->end) {
            /* Overlap: merge */
            if (nxt->end > cur->end) cur->end = nxt->end;
        } else {
            /* No overlap: advance */
            j++;
            if (i != j) regs->a[j] = *nxt;
        }
    }
    regs->n = j + 1;

    fprintf(stderr, "[index-capture]   merged to %zu regions\n", regs->n);
}


/* ── Region printing (verbose) ─────────────────────────────────────── */

/*
 * Print merged regions as PAC-space intervals.
 * PAC pos = chrom_offset + start (forward strand, 0..l_pac-1).
 */
static void
capt_regions_print(const capt_regions_v *regs)
{
    for (size_t i = 0; i < regs->n; i++) {
        const capt_region_t *r = &regs->a[i];
        fprintf(stderr, "  [%zu] chrom_off=%lld start=%lld end=%lld  "
                "→ pac=[%lld, %lld)  len=%lld\n",
                i,
                (long long)r->chrom_offset,
                (long long)r->start,
                (long long)r->end,
                (long long)(r->chrom_offset + r->start),
                (long long)(r->chrom_offset + r->end),
                (long long)(r->end - r->start));
    }
}


/* ── Dump merged regions as BED ────────────────────────────────────── */

static void capt_dump_regions(const capt_regions_v *regs, const bntseq_t *bns,
                               const char *dir)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/region.txt", dir);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "[index-capture] ERROR: cannot write '%s'\n", path);
        return;
    }
    for (size_t i = 0; i < regs->n; i++) {
        const capt_region_t *r = &regs->a[i];
        /* find chr name from offset */
        const char *chr = "unknown";
        for (int j = 0; j < bns->n_seqs; j++) {
            if (bns->anns[j].offset == r->chrom_offset) {
                chr = bns->anns[j].name;
                break;
            }
        }
        fprintf(fp, "%s\t%ld\t%ld\tregion_%zu\n",
                chr, (long)r->start, (long)r->end, i);
    }
    fclose(fp);
    fprintf(stderr, "[index-capture]   dumped %zu merged regions to %s\n",
            regs->n, path);
}

/* ── Step 4: Multi-string SA ───────────────────────────────────────── */

/*
 * Gather region sequences (fwd+RC) as uint8_t**, each terminated by sentinel=4.
 * SA entries are packed: (seq_idx << 32) | offset.
 * Returns packed SA (sorted), n_sa, and end_cnt[4].
 */

static uint8_t **g_ms_seqs;       /* region sequences for comparator */
static int64_t  *g_ms_seq_lens;   /* length of each sequence (excluding $) */
static int       g_ms_k;          /* max comparison depth */
static int capt_ms_cmp(const void *a, const void *b)
{
    uint64_t pa = *(const uint64_t *)a;
    uint64_t pb = *(const uint64_t *)b;
    int si_a = (int)(pa >> 32), off_a = (int)(pa & 0xffffffff);
    int si_b = (int)(pb >> 32), off_b = (int)(pb & 0xffffffff);
    const uint8_t *A = g_ms_seqs[si_a] + off_a;
    const uint8_t *B = g_ms_seqs[si_b] + off_b;
    int rem_a = (int)g_ms_seq_lens[si_a] - off_a;
    int rem_b = (int)g_ms_seq_lens[si_b] - off_b;
    int maxd = rem_a < rem_b ? rem_a : rem_b;
    if (maxd > g_ms_k) maxd = g_ms_k;

    for (int d = 0; d < maxd; d++) {
        int va = (int)A[d], vb = (int)B[d];  /* $=4 > T, sorts last */
        if (va != vb) return (va < vb) ? -1 : 1;
        if (va == 4) break;  /* both hit sentinel */
    }
    if (rem_a != rem_b) return (rem_a < rem_b) ? -1 : 1;
    return 0;
}

static uint64_t *
capt_build_sa_multistr(const uint8_t *pac, int64_t l_pac,
                        const capt_regions_v *regs, int k,
                        uint64_t *n_sa, uint64_t end_cnt[4],
                        int64_t **t2p_out)
{
    int n_regions = (int)regs->n;
    int n_seqs = n_regions * 2;  /* fwd + RC */

    /* gather sequences and count suffixes */
    g_ms_seqs     = calloc(n_seqs, sizeof(uint8_t *));
    g_ms_seq_lens = calloc(n_seqs, sizeof(int64_t));
    int64_t *t2p  = NULL;  /* total text2pac (used later for pos table) */
    uint64_t total_suf = 0;

    for (int i = 0; i < n_regions; i++) {
        const capt_region_t *r = &regs->a[i];
        int64_t ps = r->chrom_offset + r->start;
        int64_t pe = r->chrom_offset + r->end;
        int64_t len = pe - ps;
        if (len <= 0) continue;

        /* forward */
        uint8_t *fwd = malloc((size_t)len + 1);
        for (int64_t p = 0; p < len; p++) {
            int64_t pp = ps + p;
            fwd[p] = (uint8_t)((pac[pp >> 2] >> ((~pp & 3) << 1)) & 3);
        }
        fwd[len] = 4;
        g_ms_seqs[i * 2]     = fwd;
        g_ms_seq_lens[i * 2]  = len + 1;   /* +1 for sentinel */
        total_suf += (uint64_t)(len + 1);
        end_cnt[fwd[len - 1]]++;

        /* RC */
        uint8_t *rc = malloc((size_t)len + 1);
        for (int64_t p = 0; p < len; p++) {
            int64_t pp = pe - 1 - p;
            uint8_t b = (uint8_t)((pac[pp >> 2] >> ((~pp & 3) << 1)) & 3);
            rc[p] = 3 - b;
        }
        rc[len] = 4;
        g_ms_seqs[i * 2 + 1]     = rc;
        g_ms_seq_lens[i * 2 + 1]  = len + 1;   /* +1 for sentinel */
        total_suf += (uint64_t)(len + 1);
        end_cnt[rc[len - 1]]++;
    }
    *n_sa = total_suf;

    if (bwa_verbose >= 3) {
        fprintf(stderr, "[index-capture]   %d seqs, %llu suffixes\n",
                n_seqs, (unsigned long long)total_suf);
        fprintf(stderr, "[index-capture]   end_cnt: A=%llu C=%llu G=%llu T=%llu\n",
                (unsigned long long)end_cnt[0], (unsigned long long)end_cnt[1],
                (unsigned long long)end_cnt[2], (unsigned long long)end_cnt[3]);
    }

    /* build suffix positions (packed) */
    uint64_t *sa_pos = malloc((size_t)total_suf * sizeof(uint64_t));
    uint64_t idx = 0;
    for (int si = 0; si < n_seqs; si++) {
        if (!g_ms_seqs[si]) continue;
        int64_t slen = g_ms_seq_lens[si];
        for (int p = 0; p < slen; p++) {
            sa_pos[idx++] = ((uint64_t)(uint32_t)si << 32) | (uint64_t)(uint32_t)p;
        }
    }

    /* sort */
    g_ms_k = k;
    if (bwa_verbose >= 3)
        fprintf(stderr, "[index-capture]   sorting...\n");
    qsort(sa_pos, (size_t)total_suf, sizeof(uint64_t), capt_ms_cmp);
    g_ms_k = 0;

    /* build t2p array for position table later: for each SA entry, map to PAC position */
    t2p = malloc((size_t)total_suf * sizeof(int64_t));
    for (uint64_t i = 0; i < total_suf; i++) {
        int si = (int)(sa_pos[i] >> 32);
        int off = (int)(sa_pos[i] & 0xffffffff);
        /* reconstruct PAC position from sequence index and offset */
        int ri = si / 2;
        int is_rc = si & 1;
        const capt_region_t *r = &regs->a[ri];
        int64_t ps = r->chrom_offset + r->start;
        if (is_rc) {
            int64_t pe = r->chrom_offset + r->end;
            t2p[i] = l_pac + (pe - 1 - off);  /* RC PAC position */
        } else {
            t2p[i] = ps + off;  /* forward PAC position */
        }
    }

    *t2p_out = t2p;
    g_ms_seq_lens = NULL;  /* keep g_ms_seqs alive for BWT build */
    return sa_pos;  /* sorted, packed SA */
}

/* ── Dump SA to file (150 chars) ───────────────────────────────────── */

static void capt_dump_sa_ms(const uint64_t *sa, uint64_t n_sa,
                             const char *dir)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/sa_dump.txt", dir);
    FILE *fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "[index-capture] ERROR: cannot write '%s'\n", path); return; }
    for (uint64_t i = 0; i < n_sa; i++) {
        int si = (int)(sa[i] >> 32);
        int off = (int)(sa[i] & 0xffffffff);
        const uint8_t *s = g_ms_seqs[si] + off;
        for (int d = 0; d < 150; d++) {
            if (s[d] >= 4) break;
            fputc("ACGT"[s[d]], fp);
        }
        fputc('\n', fp);
    }
    fclose(fp);
    fprintf(stderr, "[index-capture]   dumped %llu SA entries to %s\n",
            (unsigned long long)n_sa, path);
}

/* ── Step 5: Build BWT + OCC from concat text + SA ─────────────────── */

static void
capt_build_bwt_occ_ms(const uint64_t *sa, uint64_t n_sa,
                       uint8_t **bwt_out, uint64_t **occ_out)
{
    uint64_t n_blocks = n_sa / CAPT_OCC_INTV + 1;
    uint8_t  *bwt = calloc(n_sa, 1);
    uint64_t *occ[4];
    for (int c = 0; c < 4; c++)
        occ[c] = calloc(n_blocks, 8);

    uint64_t cnt[4] = {0, 0, 0, 0};

    for (uint64_t i = 0; i < n_sa; i++) {
        int si = (int)(sa[i] >> 32);
        int off = (int)(sa[i] & 0xffffffff);
        int base = (off == 0) ? 4 : (int)g_ms_seqs[si][off - 1];
        bwt[i] = (uint8_t)base;
        if (base < 4) cnt[base]++;
        if ((i + 1) % CAPT_OCC_INTV == 0) {
            uint64_t blk = (i + 1) / CAPT_OCC_INTV;
            for (int c = 0; c < 4; c++)
                occ[c][blk] = cnt[c];
        }
    }

    *bwt_out = bwt;
    for (int c = 0; c < 4; c++) occ_out[c] = occ[c];

    if (bwa_verbose >= 3) {
        uint64_t n_sentinel = n_sa - cnt[0] - cnt[1] - cnt[2] - cnt[3];
        fprintf(stderr, "[index-capture]   BWT: A=%llu C=%llu G=%llu T=%llu sentinel=%llu\n",
                (unsigned long long)cnt[0], (unsigned long long)cnt[1],
                (unsigned long long)cnt[2], (unsigned long long)cnt[3],
                (unsigned long long)n_sentinel);
    }
}

/* ── Dump BWT stats to file ────────────────────────────────────────── */

static void capt_dump_bwt(const uint8_t *bwt, uint64_t n_sa,
                           uint64_t *occ[4], const char *dir)
{
    char path[1024];

    /* --- BWT --- */
    snprintf(path, sizeof(path), "%s/bwt_dump.txt", dir);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "[index-capture] ERROR: cannot write '%s'\n", path);
        return;
    }
    uint64_t cnt[5] = {0, 0, 0, 0, 0};
    for (uint64_t i = 0; i < n_sa; i++)
        cnt[bwt[i] < 5 ? bwt[i] : 4]++;
    fprintf(fp, "# BWT length: %llu\n", (unsigned long long)n_sa);
    fprintf(fp, "# A=%llu C=%llu G=%llu T=%llu sentinel=%llu\n",
            (unsigned long long)cnt[0], (unsigned long long)cnt[1],
            (unsigned long long)cnt[2], (unsigned long long)cnt[3],
            (unsigned long long)cnt[4]);
    for (uint64_t i = 0; i < n_sa; i++)
        fprintf(fp, "%c\n", bwt[i] < 4 ? "ACGT"[bwt[i]] : 'N');
    fclose(fp);
    fprintf(stderr, "[index-capture]   dumped BWT to %s\n", path);

    /* --- OCC (complete table) --- */
    snprintf(path, sizeof(path), "%s/occ_dump.txt", dir);
    fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "[index-capture] ERROR: cannot write '%s'\n", path);
        return;
    }
    uint64_t n_blocks = n_sa / CAPT_OCC_INTV + 1;
    fprintf(fp, "# OCC blocks: %llu (interval=%d)\n",
            (unsigned long long)n_blocks, CAPT_OCC_INTV);
    fprintf(fp, "# blk\tA\tC\tG\tT\n");
    for (uint64_t b = 0; b < n_blocks; b++)
        fprintf(fp, "%llu\t%llu\t%llu\t%llu\t%llu\n",
                (unsigned long long)b,
                (unsigned long long)occ[0][b],
                (unsigned long long)occ[1][b],
                (unsigned long long)occ[2][b],
                (unsigned long long)occ[3][b]);
    fclose(fp);
    fprintf(stderr, "[index-capture]   dumped OCC to %s\n", path);
}

/* ── Step 6: Build LCP ─────────────────────────────────────────────── */

/*
 * Build LCP via Kasai.  Compares k-mers from genome PAC.
 *
 * Returns malloc'd uint32_t array, LCP[0]=0.
 */
static uint32_t *
capt_build_lcp_ms(const uint64_t *sa, uint64_t n_sa, int k)
{
    if (n_sa == 0) return NULL;
    uint32_t lcp_max = 0;
    uint32_t *lcp = calloc(n_sa, sizeof(uint32_t));
    lcp[0] = 0;

    if (bwa_verbose >= 3)
        fprintf(stderr, "[index-capture]   computing LCP for %llu suffixes...\n",
                (unsigned long long)n_sa);

    for (uint64_t i = 1; i < n_sa; i++) {
        int si_a = (int)(sa[i - 1] >> 32), off_a = (int)(sa[i - 1] & 0xffffffff);
        int si_b = (int)(sa[i] >> 32),      off_b = (int)(sa[i] & 0xffffffff);
        const uint8_t *A = g_ms_seqs[si_a] + off_a;
        const uint8_t *B = g_ms_seqs[si_b] + off_b;
        int rem_a = (int)g_ms_seq_lens[si_a] - off_a;
        int rem_b = (int)g_ms_seq_lens[si_b] - off_b;
        int maxd = rem_a < rem_b ? rem_a : rem_b;
        if (maxd > k) maxd = k;
        int h = 0;
        while (h < maxd && A[h] == B[h] && A[h] < 4)
            h++;
        lcp[i] = (uint32_t)h;
        if (h > (int)lcp_max) lcp_max = (uint32_t)h;
    }

    if (bwa_verbose >= 3) {
        uint64_t sum_lcp = 0;
        for (uint64_t j = 0; j < n_sa; j++) sum_lcp += lcp[j];
        fprintf(stderr, "[index-capture]   LCP: max=%u avg=%.1f\n",
                lcp_max, sum_lcp / (double)(n_sa > 1 ? n_sa - 1 : 1));
    }
    return lcp;
}

/* ── Dump LCP to file ──────────────────────────────────────────────── */

static void capt_dump_lcp(const uint32_t *lcp, uint64_t n_sa,
                           const char *dir)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/lcp_dump.txt", dir);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "[index-capture] ERROR: cannot write '%s'\n", path);
        return;
    }
    fprintf(fp, "# LCP length: %llu\n", (unsigned long long)n_sa);
    for (uint64_t i = 0; i < n_sa; i++)
        fprintf(fp, "%llu\t%u\n", (unsigned long long)i, lcp[i]);
    fclose(fp);
    fprintf(stderr, "[index-capture]   dumped LCP to %s\n", path);
}

/* ── Step 7: Build tail-RLE ────────────────────────────────────────── */

/*
 * Load genome BWT index from prefix (e.g. "chr7/chr7.fa" → .bwt, .sa).
 * Returns bwt_t with SA loaded, or NULL on error.
 */
static bwt_t *
capt_load_genome_bwt(const char *prefix)
{
    char buf[1024];
    bwt_t *bwt;
    snprintf(buf, sizeof(buf), "%s.bwt", prefix);
    bwt = bwt_restore_bwt(buf);
    if (!bwt) return NULL;
    snprintf(buf, sizeof(buf), "%s.sa", prefix);
    bwt_restore_sa(buf, bwt);
    return bwt;
}

/*
 * Compute tail-RLE vectors for each sub-SA suffix using genome BWT.
 *
 * For suffix i: extend k-mer on genome BWT, record genome-SA rank as
 * f(query_depth).  Store only runs for depths LCP[i]+1 .. k (the tail);
 * depths 1..LCP[i] are inherited from suffix i-1 at load time.
 *
 * Each RLE entry is 9 bytes: 8B val (genome-SA rank) + 1B len (run length).
 *
 * Returns: rle_data (uint64_t* alias for 9-byte-packed buffer),
 *          rle_offsets (per-suffix byte index), n_rle (total entries).
 */
static void
capt_build_rle_ms(const uint64_t *sa, uint64_t n_sa, int k,
                   const uint32_t *lcp,
                   const bwt_t *g_bwt,
                   uint64_t **rle_data, uint64_t **rle_offsets,
                   uint64_t *n_rle)
{
    size_t rle_cap = n_sa * 12 * 9;
    uint8_t  *rle_bytes = calloc(rle_cap, 1);
    uint64_t *off_bytes = calloc(n_sa, 8);
    uint64_t  byte_pos = 0;

    if (bwa_verbose >= 3)
        fprintf(stderr, "[index-capture]   computing RLE for %llu suffixes...\n",
                (unsigned long long)n_sa);

    for (uint64_t i = 0; i < n_sa; i++) {
        int si = (int)(sa[i] >> 32);
        int off = (int)(sa[i] & 0xffffffff);
        int share = (int)lcp[i];

        int b0 = (int)g_ms_seqs[si][off];
        off_bytes[i] = byte_pos;
        if (b0 > 3) continue;

        bwtintv_t ik;
        bwt_set_intv(g_bwt, b0, ik);

        uint64_t left_vals[256], right_vals[256];
        uint8_t  left_lens[256], right_lens[256];
        int n_left = 0, n_right = 0;

        uint64_t prev_left  = ik.x[0];
        uint64_t prev_right = ik.x[0] + ik.x[2] - 1;
        int left_start = 1, right_start = 1;

        for (int d = 1; d < k; d++) {
            int b = (int)g_ms_seqs[si][off + d];
            if (b > 3) break;

            bwtintv_t ok[4];
            bwt_extend(g_bwt, &ik, ok, 0);
            ik = ok[3 - b];

            int depth = d + 1;
            uint64_t cur_left  = ik.x[0];
            uint64_t cur_right = ik.x[0] + ik.x[2] - 1;

            if (cur_left != prev_left) {
                left_vals[n_left] = prev_left;
                left_lens[n_left] = (uint8_t)(depth - left_start);
                n_left++; prev_left = cur_left; left_start = depth;
            }
            if (cur_right != prev_right) {
                right_vals[n_right] = prev_right;
                right_lens[n_right] = (uint8_t)(depth - right_start);
                n_right++; prev_right = cur_right; right_start = depth;
            }
        }
        if (left_start <= k) {
            left_vals[n_left] = prev_left;
            left_lens[n_left] = (uint8_t)(k + 1 - left_start);
            n_left++;
        }
        if (right_start <= k) {
            right_vals[n_right] = prev_right;
            right_lens[n_right] = (uint8_t)(k + 1 - right_start);
            n_right++;
        }

        int consumed = 0;
        for (int j = 0; j < n_left; j++) {
            uint8_t rlen = left_lens[j];
            if (consumed + rlen <= share) { consumed += rlen; continue; }
            if (consumed < share) { rlen -= (share - consumed); consumed = share; }
            memcpy(rle_bytes + byte_pos, &left_vals[j], 8);
            rle_bytes[byte_pos + 8] = rlen;
            byte_pos += 9;
        }
        consumed = 0;
        for (int j = 0; j < n_right; j++) {
            uint8_t rlen = right_lens[j];
            if (consumed + rlen <= share) { consumed += rlen; continue; }
            if (consumed < share) { rlen -= (share - consumed); consumed = share; }
            memcpy(rle_bytes + byte_pos, &right_vals[j], 8);
            rle_bytes[byte_pos + 8] = rlen;
            byte_pos += 9;
        }

        if ((i + 1) % 1000000 == 0 && bwa_verbose >= 3)
            fprintf(stderr, "    %llu / %llu suffixes processed\n",
                    (unsigned long long)(i + 1), (unsigned long long)n_sa);
    }

    *n_rle = byte_pos / 9;
    *rle_data = (uint64_t *)rle_bytes;
    *rle_offsets = off_bytes;

    if (bwa_verbose >= 3)
        fprintf(stderr, "[index-capture]   RLE: %llu tail entries (%llu bytes)\n",
                (unsigned long long)*n_rle, (unsigned long long)byte_pos);
}

/* ── Dump RLE to file ──────────────────────────────────────────────── */

static void capt_dump_rle(const uint64_t *rle_data, uint64_t n_rle,
                           const uint64_t *rle_offsets, uint64_t n_sa,
                           const uint32_t *lcp, int k, const char *dir)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/rle_dump.txt", dir);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "[index-capture] ERROR: cannot write '%s'\n", path);
        return;
    }
    fprintf(fp, "# n_sa=%llu n_rle=%llu k=%d\n",
            (unsigned long long)n_sa, (unsigned long long)n_rle, k);
    const uint8_t *bytes = (const uint8_t *)rle_data;
    for (uint64_t i = 0; i < n_sa; i++) {
        uint64_t off = rle_offsets[i] / 9;  /* byte offset → entry index */
        uint64_t end = (i + 1 < n_sa) ? rle_offsets[i + 1] / 9 : n_rle;
        int share = (int)lcp[i];
        int tail = k - share;

        int printed = 0, cum = 0;
        uint64_t p = off;
        while (p < end && cum < tail) {
            uint64_t val;
            uint8_t len;
            memcpy(&val, bytes + p * 9, 8);
            len = bytes[p * 9 + 8];
            if (printed++) fputc(',', fp);
            fprintf(fp, "(%llu,%u)", (unsigned long long)val, len);
            cum += len;
            p++;
        }
        fputc('|', fp);
        /* right runs */
        printed = 0; cum = 0;
        while (p < end) {
            uint64_t val;
            uint8_t len;
            memcpy(&val, bytes + p * 9, 8);
            len = bytes[p * 9 + 8];
            if (printed++) fputc(',', fp);
            fprintf(fp, "(%llu,%u)", (unsigned long long)val, len);
            cum += len;
            p++;
        }
        fprintf(fp, "\t# %llu share=%d n_entries=%llu\n",
                (unsigned long long)i, share,
                (unsigned long long)(end - off));
    }
    fclose(fp);
    fprintf(stderr, "[index-capture]   dumped RLE to %s\n", path);
}

/* ── Step 8: Build position table ──────────────────────────────────── */

#define CAPT_UNIQ_THRESH 1000

/* Read base at PAC position p.  PAC stores forward half only;
 * RC half (p ≥ l_pac) is virtual: complement of remapped position. */
static inline int
capt_pac_base(const uint8_t *pac, int64_t p, int64_t l_pac)
{
    int64_t fp = (p < l_pac) ? p : ((l_pac << 1) - 1 - p);
    int b = (pac[fp >> 2] >> ((~fp & 3) << 1)) & 3;
    return (p < l_pac) ? b : 3 - b;
}

/*
 * Inverse of bwt_sa(): given a PAC position p, extend RIGHT (forward)
 * on the genome BWT until unique (or ≤THRESH), return genome-SA rank.
 */
static uint64_t
capt_pac2rank(const bwt_t *bwt, const uint8_t *pac, int64_t p, int64_t l_pac)
{
    int b = capt_pac_base(pac, p, l_pac);
    if (b > 3) return UINT64_MAX;
    int64_t pac_end = l_pac << 1;
    bwtintv_t ik;
    bwt_set_intv(bwt, b, ik);
    for (int64_t c = p + 1; c < pac_end && ik.x[2] > CAPT_UNIQ_THRESH && c - p <= 1000000; c++) {
        int cb = capt_pac_base(pac, c, l_pac);
        if (cb > 3) break;
        bwtintv_t ok[4];
        bwt_extend(bwt, &ik, ok, 0);
        ik = ok[3 - cb];
    }
    if (ik.x[2] == 1) return ik.x[0];
    if (ik.x[2] <= CAPT_UNIQ_THRESH) {
        for (bwtint_t s = 0; s < ik.x[2]; s++) {
            bwtint_t r = ik.x[0] + s;
            if ((int64_t)bwt_sa(bwt, r) == p) return r;
        }
    }
    return ik.x[0];
}

/* One PAC interval: [lo, hi) */
typedef struct { int64_t lo, hi; } capt_iv_t;
typedef kvec_t(capt_iv_t) capt_iv_v;

/* (rank, pac) pair */
typedef struct { uint64_t rank; int64_t pac; } capt_rp_t;

static int capt_rp_cmp(const void *a, const void *b) {
    uint64_t ra = ((const capt_rp_t *)a)->rank;
    uint64_t rb = ((const capt_rp_t *)b)->rank;
    return (ra > rb) ? 1 : (ra < rb) ? -1 : 0;
}

/*
 * Build position table: for each suffix position in concat text,
 * resolve genome-SA rank via bwt_extend, store (rank, pac) pair.
 */
static void
capt_build_pos_ms(const bwt_t *g_bwt,
                   const uint64_t *sa, uint64_t n_sa,
                   const int64_t *t2p,
                   uint64_t **pos_rank, int64_t **pos_pac,
                   uint64_t *n_pos)
{
    typedef kvec_t(capt_rp_t) capt_rp_v;
    capt_rp_v rp;
    kv_init(rp);

    for (uint64_t i = 0; i < n_sa; i++) {
        int si = (int)(sa[i] >> 32);
        int off = (int)(sa[i] & 0xffffffff);
        int b = (int)g_ms_seqs[si][off];
        if (b > 3) continue;

        bwtintv_t ik;
        bwt_set_intv(g_bwt, b, ik);
        /* extend until unique */
        for (int d = 1; ik.x[2] > 1; d++) {
            int nb = (int)g_ms_seqs[si][off + d];
            if (nb > 3) break;
            bwtintv_t ok[4];
            bwt_extend(g_bwt, &ik, ok, 0);
            ik = ok[3 - nb];
            if (d > 1000000) break;
        }
        int64_t pac = t2p[i];
        if (pac < 0) continue;
        kv_push(capt_rp_t, rp, ((capt_rp_t){ ik.x[0], pac }));
    }

    if (bwa_verbose >= 3)
        fprintf(stderr, "[index-capture]   %zu position-table entries\n", rp.n);

    qsort(rp.a, rp.n, sizeof(capt_rp_t), capt_rp_cmp);

    *n_pos = rp.n;
    *pos_rank = malloc(rp.n * 8);
    *pos_pac  = malloc(rp.n * 8);
    for (size_t i = 0; i < rp.n; i++) {
        (*pos_rank)[i] = rp.a[i].rank;
        (*pos_pac)[i]  = rp.a[i].pac;
    }
    free(rp.a);
}

/* ── Dump position table to file ───────────────────────────────────── */

static void capt_dump_pos(const uint64_t *pos_rank, const int64_t *pos_pac,
                           uint64_t n_pos, const char *dir)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/pos_dump.txt", dir);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "[index-capture] ERROR: cannot write '%s'\n", path);
        return;
    }
    fprintf(fp, "# n_pos=%llu\n", (unsigned long long)n_pos);
    fprintf(fp, "# rank\tpac\n");
    for (uint64_t i = 0; i < n_pos; i++)
        fprintf(fp, "%llu\t%lld\n",
                (unsigned long long)pos_rank[i], (long long)pos_pac[i]);
    fclose(fp);
    fprintf(stderr, "[index-capture]   dumped position table to %s\n", path);
}

/* ── Sanity check: position table ≡ inverse of bwt_sa ──────────────── */

static void capt_pos_sanity(const bwt_t *g_bwt, const uint8_t *pac,
                             int64_t l_pac,
                             const uint64_t *pos_rank,
                             const int64_t *pos_pac, uint64_t n_pos,
                             const char *dir)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/pos_sanity.txt", dir);
    FILE *fp = fopen(path, "w");
    if (!fp) return;

    fprintf(fp, "# sanity: pos_table[i] == bwt_sa(g_bwt, pos_rank[i])\n");
    fprintf(fp, "# checking all %llu entries...\n", (unsigned long long)n_pos);

    uint64_t ngood = 0, nbad = 0;

    for (uint64_t i = 0; i < n_pos; i++) {
        bwtint_t sa_val = bwt_sa(g_bwt, pos_rank[i]);
        int64_t expected = pos_pac[i];
        if ((int64_t)sa_val == expected) {
            ngood++;
        } else {
            nbad++;
            fprintf(fp, "MISMATCH  i=%llu rank=%llu  expected_pac=%lld  bwt_sa=%lld\n",
                    (unsigned long long)i,
                    (unsigned long long)pos_rank[i],
                    (long long)expected, (long long)sa_val);
        }
    }

    fprintf(fp, "# result: %llu OK  %llu MISMATCH\n",
            (unsigned long long)ngood, (unsigned long long)nbad);
    fclose(fp);
    fprintf(stderr, "[index-capture]   pos sanity: %llu OK, %llu MISMATCH → %s\n",
            (unsigned long long)ngood, (unsigned long long)nbad,
            nbad ? "FAIL" : "PASS");
}

/* ── Serialize .capt file ──────────────────────────────────────────── */

/*
 * Write the complete .capt binary file.
 *
 * Layout: 128-byte header → BWT → SA → OCC → LCP →
 *         RLE data → RLE offsets → pos ranks → pos coords
 */
static void
capt_dump(const char *prefix, const capt_t *capt)
{
    char fn[1024];
    snprintf(fn, sizeof(fn), "%s.capt", prefix);

    uint64_t n_blocks = capt->n_sa / CAPT_OCC_INTV + 1;

    /* section sizes */
    uint64_t sz_bwt    = capt->n_sa;
    uint64_t sz_sa     = capt->n_sa * 8;
    uint64_t sz_occ    = 4 * n_blocks * 8;
    uint64_t sz_lcp    = capt->n_sa * 4;
    uint64_t sz_rle    = capt->n_rle * 9;
    uint64_t sz_rleoff = capt->n_rle ? capt->n_sa * 8 : 0;
    uint64_t sz_posr   = capt->n_pos * 8;
    uint64_t sz_posp   = capt->n_pos * 8;

    /* offsets */
    uint64_t off = CAPT_HDR_SZ;
    uint64_t off_bwt    = off; off += sz_bwt;
    uint64_t off_sa     = off; off += sz_sa;
    uint64_t off_occ    = off; off += sz_occ;
    uint64_t off_lcp    = off; off += sz_lcp;
    uint64_t off_rle    = off; off += sz_rle;
    uint64_t off_rleof  = off; off += sz_rleoff;
    uint64_t off_posr   = off; off += sz_posr;
    uint64_t off_posp   = off; off += sz_posp;
    uint64_t off_endcnt = off; off += 32;  /* 4 × uint64_t */

    /* header */
    capt_hdr_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    memcpy(hdr.magic, "BWACAPT", 7);
    hdr.version       = CAPT_VERSION;
    hdr.k             = capt->k;
    hdr.padding       = capt->padding;
    hdr.n_sa          = capt->n_sa;
    hdr.n_rle         = capt->n_rle;
    hdr.n_pos         = capt->n_pos;
    hdr.occ_interval  = CAPT_OCC_INTV;
    hdr.off_bwt       = off_bwt;
    hdr.off_sa        = off_sa;
    hdr.off_occ       = off_occ;
    hdr.off_enc       = off_endcnt;  /* repurposed: end_cnt offset */
    hdr.off_lcp       = off_lcp;
    hdr.off_rle_data  = off_rle;
    hdr.off_rle_off   = off_rleof;
    hdr.off_pos_rank  = off_posr;
    hdr.off_pos_genome = off_posp;

    FILE *fp = fopen(fn, "wb");
    if (!fp) {
        fprintf(stderr, "[index-capture] ERROR: cannot write '%s'\n", fn);
        return;
    }
    fwrite(&hdr, CAPT_HDR_SZ, 1, fp);
    fwrite(capt->bwt, 1, sz_bwt, fp);
    fwrite(capt->sa, 8, capt->n_sa, fp);
    for (int c = 0; c < 4; c++)
        fwrite(capt->occ[c], 8, n_blocks, fp);
    fwrite(capt->lcp, 4, capt->n_sa, fp);
    if (capt->rle_data)    fwrite(capt->rle_data, 9, capt->n_rle, fp);
    if (capt->rle_offsets) fwrite(capt->rle_offsets, 8, capt->n_sa, fp);
    if (capt->pos_rank)    fwrite(capt->pos_rank, 8, capt->n_pos, fp);
    if (capt->pos_pac)     fwrite(capt->pos_pac, 8, capt->n_pos, fp);
    fwrite(capt->end_cnt, 8, 4, fp);
    fclose(fp);

    uint64_t total = off_endcnt + 32;
    fprintf(stderr, "[index-capture]   wrote %s (%.1f MB)\n",
            fn, (double)total / 1e6);
}

/* ── Top-level command ──────────────────────────────────────────────── */

static void usage(void)
{
    fprintf(stderr, "Usage: bwa index-capture <ref.fa> <target.bed> [options]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -p INT     padding around each region [0]\n");
    fprintf(stderr, "  -k INT     k-mer size, > max read length [200]\n");
    fprintf(stderr, "  -o STR     output prefix (writes prefix.capt)\n");
    fprintf(stderr, "  -t INT     number of threads [16]\n");
    fprintf(stderr, "  -d STR     dump merged regions to STR/region.txt\n");
}

int main_index_capture(int argc, char *argv[])
{
    int c, k = 200, padding = 0;
    const char *output = NULL, *dumpdir = NULL;

    while ((c = getopt(argc, argv, "p:k:o:t:d:")) >= 0) {
        switch (c) {
        case 'p': padding = atoi(optarg); break;
        case 'k': k = atoi(optarg); break;
        case 'o': output = optarg; break;
        case 't': /* n_threads */ break;
        case 'd': dumpdir = optarg; break;
        default:  usage(); return 1;
        }
    }

    if (optind + 2 > argc) {
        fprintf(stderr, "[index-capture] missing ref.fa and/or target.bed\n");
        usage();
        return 1;
    }
    if (!output) {
        fprintf(stderr, "[index-capture] -o <output_prefix> is required\n");
        usage();
        return 1;
    }

    const char *ref_path = argv[optind];
    const char *bed_path = argv[optind + 1];

    fprintf(stderr, "[index-capture] ref=%s bed=%s padding=%d k=%d output=%s.capt\n",
            ref_path, bed_path, padding, k, output);

    /* Allocate capture index */
    capt_t *capt = calloc(1, sizeof(capt_t));
    capt->k = k;
    capt->padding = padding;

    /* Load genome annotation */
    bntseq_t *bns = bns_restore(ref_path);
    if (!bns) {
        fprintf(stderr, "[index-capture] ERROR: cannot load reference '%s'\n", ref_path);
        free(capt); return 1;
    }

    /* 1. Parse BED */
    fprintf(stderr, "[index-capture] Step 1/8: Parsing BED...\n");
    capt_regions_v *regs = capt_bed_read(bed_path, bns);

    /* 2. Pad each region */
    fprintf(stderr, "[index-capture] Step 2/8: Padding regions (+%d bp)...\n", padding);
    capt_bed_pad(regs, padding);

    /* 3. Merge padded regions */
    fprintf(stderr, "[index-capture] Step 3/8: Merging padded regions...\n");
    capt_bed_merge(regs);

    /* Dump merged regions if requested */
    if (dumpdir) capt_dump_regions(regs, bns, dumpdir);

    /* Print merged intervals */
    //capt_regions_print(regs);

    /* Load genome PAC */
    char buf[1024];
    snprintf(buf, sizeof(buf), "%s.pac", ref_path);
    FILE *fp_pac = fopen(buf, "rb");
    if (!fp_pac) {
        fprintf(stderr, "[index-capture] ERROR: cannot open '%s'\n", buf);
        free(regs->a); free(regs); bns_destroy(bns); free(capt);
        return 1;
    }
    fseek(fp_pac, 0, SEEK_END);
    int64_t pac_size = ftell(fp_pac);
    uint8_t *pac = malloc(pac_size);
    fseek(fp_pac, 0, SEEK_SET);
    fread(pac, 1, pac_size, fp_pac);
    fclose(fp_pac);

    /* 4. Build multi-string SA */
    fprintf(stderr, "[index-capture] Step 4/8: Building multi-string SA...\n");
    uint64_t end_cnt[4] = {0, 0, 0, 0};
    int64_t *t2p = NULL;
    capt->sa = capt_build_sa_multistr(pac, bns->l_pac, regs, k,
                                       &capt->n_sa, end_cnt, &t2p);
    fprintf(stderr, "[index-capture]   SA: %llu entries\n",
            (unsigned long long)capt->n_sa);
    fprintf(stderr, "[index-capture]   end_cnt: A=%llu C=%llu G=%llu T=%llu\n",
            (unsigned long long)end_cnt[0], (unsigned long long)end_cnt[1],
            (unsigned long long)end_cnt[2], (unsigned long long)end_cnt[3]);
    memcpy(capt->end_cnt, end_cnt, sizeof(end_cnt));

    /* Dump SA if requested */
    if (dumpdir) capt_dump_sa_ms(capt->sa, capt->n_sa, dumpdir);

    /* Recompute g_ms_seq_lens for LCP (was set to NULL in build_sa) */
    {
        int n_seqs = (int)regs->n * 2;
        g_ms_seq_lens = calloc(n_seqs, sizeof(int64_t));
        for (int i = 0; i < n_seqs; i++) {
            if (!g_ms_seqs[i]) continue;
            int64_t slen = 0;
            while (g_ms_seqs[i][slen] < 4) slen++;
            g_ms_seq_lens[i] = slen;
        }
    }

    /* 5. Build BWT + OCC */
    fprintf(stderr, "[index-capture] Step 5/8: Building BWT+OCC...\n");
    capt_build_bwt_occ_ms(capt->sa, capt->n_sa, &capt->bwt, capt->occ);

    /* Dump BWT if requested */
    if (dumpdir) capt_dump_bwt(capt->bwt, capt->n_sa, capt->occ, dumpdir);

    /* 6. Build LCP */
    fprintf(stderr, "[index-capture] Step 6/8: Building LCP...\n");
    capt->lcp = capt_build_lcp_ms(capt->sa, capt->n_sa, k);

    /* Dump LCP if requested */
    if (dumpdir) capt_dump_lcp(capt->lcp, capt->n_sa, dumpdir);

    /* 7. Load genome BWT, build RLE */
    fprintf(stderr, "[index-capture] Step 7/8: Computing RLE...\n");
    bwt_t *g_bwt = capt_load_genome_bwt(ref_path);
    if (!g_bwt) {
        fprintf(stderr, "[index-capture] ERROR: cannot load genome BWT\n");
        free(t2p); capt_destroy(capt); free(pac);
        free(regs->a); free(regs); bns_destroy(bns);
        return 1;
    }
    capt_build_rle_ms(capt->sa, capt->n_sa, k, capt->lcp, g_bwt,
                       &capt->rle_data, &capt->rle_offsets, &capt->n_rle);

    /* Dump RLE if requested */
    if (dumpdir) capt_dump_rle(capt->rle_data, capt->n_rle, capt->rle_offsets,
                                capt->n_sa, capt->lcp, k, dumpdir);

    /* 8. Build position table */
    fprintf(stderr, "[index-capture] Step 8/8: Computing position table...\n");
    capt_build_pos_ms(g_bwt, capt->sa, capt->n_sa, t2p,
                       &capt->pos_rank, &capt->pos_pac, &capt->n_pos);

    /* Dump position table if requested */
    if (dumpdir) {
        capt_dump_pos(capt->pos_rank, capt->pos_pac, capt->n_pos, dumpdir);
    }

    free(t2p);
    bwt_destroy(g_bwt);

    /* Serialize .capt */
    capt_dump(output, capt);

    /* Load back and verify if requested */
    if (dumpdir) {
        char capt_fn[1024];
        snprintf(capt_fn, sizeof(capt_fn), "%s.capt", output);
        capt_t *loaded = capt_restore(capt_fn);
        if (loaded) {
            capt_load_verify(capt, loaded, dumpdir);
            capt_destroy(loaded);
        }
    }

    /* Final cleanup */
    free(g_ms_seq_lens);
    for (int i = 0; i < (int)regs->n * 2; i++) free(g_ms_seqs[i]);
    free(g_ms_seqs);
    capt_destroy(capt);
    free(pac);
    free(regs->a);
    free(regs);
    bns_destroy(bns);
    return 0;
}
