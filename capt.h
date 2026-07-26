/*
 * capt.h — Capture sub-index data structures, constants, and API.
 *
 * Shared by construction (bwa_index_capt.c) and runtime (capt.c).
 */

#ifndef CAPT_H
#define CAPT_H

#include <stdint.h>
#include "bwt.h"
#include "kvec.h"

/* forward-declare for bwamem.h (which includes us back after defining mem_alnreg_t) */
typedef struct capt_t capt_t;
#include "bwamem.h"

#define CAPT_MAGIC      "BWACAPT"
#define CAPT_VERSION    1
#define CAPT_HDR_SZ     128
#define CAPT_OCC_INTV   1

/* ── File header (packed, 128 bytes) ────────────────────────────────── */

#pragma pack(push, 1)
typedef struct {
    uint8_t  magic[8];       /* "BWACAPT\0" */
    uint32_t version;        /* 1 */
    int32_t  k;              /* k-mer size */
    int32_t  padding;        /* padding applied */
    uint64_t n_sa;           /* number of suffixes */
    uint64_t n_rle;          /* total RLE entries */
    uint64_t n_pos;          /* position table entries */
    uint32_t occ_interval;   /* 16 */
    uint32_t flags;          /* reserved */
    uint64_t off_bwt;        /* dense BWT (n_sa bytes) */
    uint64_t off_sa;         /* dense SA  (n_sa × 8) */
    uint64_t off_occ;        /* dense OCC (4 × n_blocks × 8) */
    uint64_t off_enc;        /* unused */
    uint64_t off_lcp;        /* LCP array (n_sa × 4) */
    uint64_t off_rle_data;   /* RLE entries (n_rle × 9) */
    uint64_t off_rle_off;    /* RLE per-suffix offsets (n_sa × 8) */
    uint64_t off_pos_rank;   /* position table ranks (n_pos × 8) */
    uint64_t off_pos_genome; /* position table coords (n_pos × 8) */
    int32_t  rlen_lo;        /* min read length */
    int32_t  rlen_hi;        /* max read length */
    uint64_t n_perfect;       /* number of precomputed perfect-match entries */
    uint64_t off_perfect;     /* offset to perfect-match table */
} capt_hdr_t;
#pragma pack(pop)

/* ── Precomputed alignments per (SA_idx, read_len) ──────────────────── */

typedef struct {            /* index entry: maps (SA_idx, read_len) → results */
    uint64_t offset;        /* offset into flat results array */
    uint32_t count;         /* number of mem_alnreg_t at this offset */
    uint32_t pad;
} capt_perfect_idx_t;

/* ── Data types for construction ───────────────────────────────────── */

/* One BED region: chrom is PAC-space offset on that chromosome */
typedef struct {
    int64_t  chrom_offset;
    int64_t  chrom_len;
    int64_t  start;
    int64_t  end;
} capt_region_t;

typedef kvec_t(capt_region_t) capt_regions_v;

/* ── Full in-memory capture index ──────────────────────────────────── */

struct capt_t {
    int       k;             /* k-mer size */
    int       padding;       /* padding applied */
    int       rlen_lo;       /* min read length */
    int       rlen_hi;       /* max read length */

    uint64_t  n_sa;          /* sub-SA size */
    uint64_t *sa;            /* dense suffix array (PAC positions) */
    uint8_t  *bwt;           /* dense BWT (1 byte/char) */
    uint64_t *occ[4];        /* dense occurrence table */
    uint64_t  l2[5];         /* cumulative counts (like bwt->L2) */
    uint32_t *lcp;           /* LCP array */

    uint64_t  n_rle;         /* total RLE entries */
    uint64_t *rle_data;      /* RLE entries (8B val + 1B len each) */
    uint64_t *rle_offsets;   /* per-suffix byte offset into rle_data */

    uint64_t  n_pos;         /* position table size */
    uint64_t *pos_rank;      /* genome SA ranks (sorted ascending) */
    int64_t  *pos_pac;       /* PAC coordinates (same order) */

    uint64_t  end_cnt[4];    /* bases at end of each region seq (fwd+RC) */

    capt_perfect_idx_t *perfect_idx; /* index: n_sa × n_lens entries */
    mem_alnreg_t       *perfect;     /* flat results: all alignments */
    uint64_t  n_perfect;             /* total results in perfect[] */
    uint64_t  n_perfect_idx;         /* = n_sa * n_lens */
};

/* O(1) lookup: returns (results, count) for (sa_idx, read_len) */
static inline mem_alnreg_t *
capt_perfect_get(const capt_t *capt, uint64_t sa_idx, int rlen, uint32_t *count) {
    int off = rlen - capt->rlen_lo;
    int n_lens = capt->rlen_hi - capt->rlen_lo + 1;
    if (off < 0 || off >= n_lens) { *count = 0; return NULL; }
    capt_perfect_idx_t *idx = &capt->perfect_idx[sa_idx * n_lens + off];
    if (idx->count == 0) { *count = 0; return NULL; }
    *count = idx->count;
    return &capt->perfect[idx->offset];
}

/* ── Capture-aware interval ────────────────────────────────────────── */

typedef struct {
    bwtint_t   x[3], info;  // genome-SA (cast-compatible with bwtintv_t)
    bwtintv_t  sub;         // sub-SA  (always valid)
    uint8_t     on_genome,depth;   // 1 = x[] valid, 0 = only sub valid
} capt_intv_t;

typedef kvec_t(capt_intv_t) capt_intv_v;

/* ── Runtime API (capt.c) ──────────────────────────────────────────── */

capt_t *capt_restore(const char *fn);
void    capt_destroy(capt_t *c);
uint64_t capt_occ(const capt_t *capt, int c, uint64_t k);
uint64_t capt_rle_left(const capt_t *capt, uint64_t idx, int depth);
uint64_t capt_rle_right(const capt_t *capt, uint64_t idx, int depth);
void    capt_load_verify(const capt_t *orig, const capt_t *loaded, const char *dir);

void    capt_set_intv(const capt_t *capt, int c, bwtintv_t *ik);
void    capt_extend(const capt_t *capt, const bwtintv_t *ik,
                    bwtintv_t ok[4], int is_back);
void    capt_smem1_init(capt_intv_t *ik, const capt_t *capt,
                        const bwt_t *g_bwt, int base);
void    capt_extend_bail(const capt_t *capt, const bwt_t *g_bwt,
                         const capt_intv_t *ik, int base, capt_intv_t *ok,
                         int is_back);
void    capt_int_print(const capt_intv_t *p, const char *label);

/* Push genome-SA bi-interval: RLE-translates sub→genome if on sub-index */
int64_t capt_lookup_pos(const capt_t *capt, uint64_t rank);
void    capt_translate_to_genome(const capt_t *capt, capt_intv_t *ik);
void    capt_push_curr(bwtintv_v *curr, const capt_intv_t *ik, int depth,
                       const capt_t *capt);

/* Capture-aware SMEM (forward+backward with sub-index fallback) */
int     capt_smem1(const capt_t *capt, const bwt_t *g_bwt,
                   int len, const uint8_t *q, int x, int min_intv,
                   bwtintv_v *mem, capt_intv_v *tmpvec[2]);

#endif
