/*
 * capt.c — Capture sub-index runtime: load, query, verify.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "capt.h"
#include "utils.h"

/* ── Load .capt file ───────────────────────────────────────────────── */

capt_t *capt_restore(const char *fn)
{
    FILE *fp = fopen(fn, "rb");
    if (!fp)
    {
        fprintf(stderr, "[E::%s] cannot open '%s'\n", __func__, fn);
        return NULL;
    }

    capt_hdr_t hdr;
    if (fread(&hdr, CAPT_HDR_SZ, 1, fp) != 1)
        goto fail;
    if (memcmp(hdr.magic, "BWACAPT", 7) != 0)
    {
        fprintf(stderr, "[E::%s] bad magic in '%s'\n", __func__, fn);
        goto fail;
    }
    if (hdr.version != CAPT_VERSION)
    {
        fprintf(stderr, "[E::%s] unsupported version %u\n", __func__, hdr.version);
        goto fail;
    }

    capt_t *c = calloc(1, sizeof(capt_t));
    c->k = hdr.k;
    c->padding = hdr.padding;
    c->n_sa = hdr.n_sa;
    c->n_pos = hdr.n_pos;

    fprintf(stderr, "[M::%s] loading %s: n_sa=%llu n_pos=%llu n_rle=%llu k=%d\n",
            __func__, fn, (unsigned long long)c->n_sa,
            (unsigned long long)c->n_pos,
            (unsigned long long)hdr.n_rle, c->k);

    /* Read SA */
    c->sa = malloc(c->n_sa * 8);
    fseek(fp, hdr.off_sa, SEEK_SET);
    fread(c->sa, 8, c->n_sa, fp);

    /* Read BWT */
    c->bwt = malloc(c->n_sa);
    fseek(fp, hdr.off_bwt, SEEK_SET);
    fread(c->bwt, 1, c->n_sa, fp);

    /* Read OCC */
    uint64_t n_blocks = c->n_sa / CAPT_OCC_INTV + 1;
    fseek(fp, hdr.off_occ, SEEK_SET);
    for (int i = 0; i < 4; i++)
    {
        c->occ[i] = malloc(n_blocks * 8);
        fread(c->occ[i], 8, n_blocks, fp);
    }
    /* compute L2 from BWT directly — count all bases including sentinels */
    uint64_t cnt[5] = {0, 0, 0, 0, 0};
    for (uint64_t i = 0; i < c->n_sa; i++)
    {
        uint8_t b = c->bwt[i];
        cnt[b < 5 ? b : 4]++;
    }
    c->l2[0] = 0;
    c->l2[1] = cnt[0];
    c->l2[2] = cnt[0] + cnt[1];
    c->l2[3] = cnt[0] + cnt[1] + cnt[2];
    c->l2[4] = cnt[0] + cnt[1] + cnt[2] + cnt[3]; /* non-sentinel total */

    /* Read LCP */
    c->lcp = malloc(c->n_sa * 4);
    fseek(fp, hdr.off_lcp, SEEK_SET);
    fread(c->lcp, 4, c->n_sa, fp);

    /* Read RLE data + offsets */
    c->n_rle = hdr.n_rle;
    c->rle_data = malloc(c->n_rle * 9);
    fseek(fp, hdr.off_rle_data, SEEK_SET);
    fread(c->rle_data, 9, c->n_rle, fp);

    c->rle_offsets = malloc(c->n_sa * 8);
    fseek(fp, hdr.off_rle_off, SEEK_SET);
    fread(c->rle_offsets, 8, c->n_sa, fp);

    fprintf(stderr, "[M::%s] RLE: %llu entries (%.1f MB)\n",
            __func__, (unsigned long long)c->n_rle,
            (double)(c->n_rle * 9) / 1e6);

    /* Read position table */
    c->pos_rank = malloc(c->n_pos * 8);
    c->pos_pac = malloc(c->n_pos * 8);
    fseek(fp, hdr.off_pos_rank, SEEK_SET);
    fread(c->pos_rank, 8, c->n_pos, fp);
    fseek(fp, hdr.off_pos_genome, SEEK_SET);
    fread(c->pos_pac, 8, c->n_pos, fp);

    /* Read end_cnt (stored at off_enc) */
    if (hdr.off_enc > 0)
    {
        fseek(fp, hdr.off_enc, SEEK_SET);
        fread(c->end_cnt, 8, 4, fp);
    }

    fclose(fp);
    return c;

fail:
    if (fp)
        fclose(fp);
    return NULL;
}

/* ── Cleanup ────────────────────────────────────────────────────────── */

void capt_destroy(capt_t *c)
{
    if (!c)
        return;
    free(c->sa);
    free(c->bwt);
    for (int i = 0; i < 4; i++)
        free(c->occ[i]);
    free(c->lcp);
    free(c->rle_data);
    free(c->rle_offsets);
    free(c->pos_rank);
    free(c->pos_pac);
    free(c);
}

/* ── OCC lookup ─────────────────────────────────────────────────────── */

uint64_t capt_occ(const capt_t *capt, int c, uint64_t k)
{
    uint64_t blk = k / CAPT_OCC_INTV;
    uint64_t base = capt->occ[c][blk];
    uint64_t start = blk * CAPT_OCC_INTV;
    for (uint64_t i = start; i < k; i++)
        if (capt->bwt[i] == (uint8_t)c)
            base++;
    return base;
}

/* ── RLE lookup ─────────────────────────────────────────────────────── */

/*
 * Walk LCP chain leftward from idx to find the suffix whose tail covers
 * depth `d`.  Returns that suffix's index.
 */
static uint64_t rle_find_ancestor(const capt_t *capt, uint64_t idx, int d)
{
    while (idx > 0 && d <= (int)capt->lcp[idx])
        idx--;
    return idx;
}

static uint64_t rle_lookup_tail(const capt_t *capt, uint64_t idx,
                                int d, int is_right)
{
    const uint8_t *rle = (const uint8_t *)capt->rle_data;
    uint64_t off = capt->rle_offsets[idx] / 9; /* byte offset → entry index */
    uint64_t end = (idx + 1 < capt->n_sa) ? capt->rle_offsets[idx + 1] / 9 : capt->n_rle;
    uint64_t pos = off;
    int cum = 0, share = (int)capt->lcp[idx];
    int target = d - share; /* depth beyond LCP that we need */

    if (pos >= end)
        return 0;

    if (is_right)
    {
        /* skip left runs: they cover depths share+1 .. k */
        int to_skip = capt->k - share;
        while (pos < end && cum < to_skip)
        {
            uint8_t len = rle[pos * 9 + 8];
            cum += len;
            if (cum >= to_skip)
            {
                pos++;
                break;
            }
            pos++;
        }
        cum = 0;
        while (pos < end)
        {
            uint8_t len = rle[pos * 9 + 8];
            cum += len;
            if (cum >= target)
            {
                uint64_t val;
                memcpy(&val, rle + pos * 9, 8);
                return val;
            }
            pos++;
        }
        return 0;
    }
    /* left: walk runs until depth d is covered */
    while (pos < end)
    {
        uint8_t len = rle[pos * 9 + 8];
        cum += len;
        if (cum >= target)
        {
            uint64_t val;
            memcpy(&val, rle + pos * 9, 8);
            return val;
        }
        pos++;
    }
    return 0;
}

uint64_t capt_rle_left(const capt_t *capt, uint64_t idx, int depth)
{
    uint64_t anc = rle_find_ancestor(capt, idx, depth);
    return rle_lookup_tail(capt, anc, depth, 0);
}

uint64_t capt_rle_right(const capt_t *capt, uint64_t idx, int depth)
{
    uint64_t anc = rle_find_ancestor(capt, idx, depth);
    return rle_lookup_tail(capt, anc, depth, 1);
}

/* ── Load-time sanity: compare loaded vs original capt_t ───────────── */

void capt_load_verify(const capt_t *orig, const capt_t *loaded, const char *dir)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/capt_verify.txt", dir);
    FILE *fp = fopen(path, "w");
    if (!fp)
        return;

    fprintf(fp, "# verify: loaded == original (field-by-field)\n");
    fprintf(fp, "# n_sa=%llu n_pos=%llu n_rle=%llu k=%d\n",
            (unsigned long long)orig->n_sa, (unsigned long long)orig->n_pos,
            (unsigned long long)orig->n_rle, orig->k);

    int nbad = 0;

/* compare scalars */
#define CHK(field, fmt)                                               \
    do                                                                \
    {                                                                 \
        if (orig->field != loaded->field)                             \
        {                                                             \
            nbad++;                                                   \
            fprintf(fp, "MISMATCH %s: orig=" fmt " loaded=" fmt "\n", \
                    #field, (unsigned long long)orig->field,          \
                    (unsigned long long)loaded->field);               \
        }                                                             \
    } while (0)

    CHK(k, "%lld");
    CHK(padding, "%lld");
    CHK(n_sa, "%llu");
    CHK(n_rle, "%llu");
    CHK(n_pos, "%llu");

/* compare arrays */
#define CHK_ARRAY(name, sz)                                      \
    do                                                           \
    {                                                            \
        if (orig->name && loaded->name &&                        \
            memcmp(orig->name, loaded->name, (size_t)(sz)) != 0) \
        {                                                        \
            nbad++;                                              \
            fprintf(fp, "MISMATCH %s: data differs\n", #name);   \
        }                                                        \
    } while (0)

    CHK_ARRAY(sa, orig->n_sa * 8);
    CHK_ARRAY(bwt, orig->n_sa);
    for (int c = 0; c < 4; c++)
    {
        if (orig->occ[c] && loaded->occ[c] &&
            memcmp(orig->occ[c], loaded->occ[c],
                   (orig->n_sa / CAPT_OCC_INTV + 1) * 8) != 0)
        {
            nbad++;
            fprintf(fp, "MISMATCH occ[%d]: data differs\n", c);
        }
    }
    CHK_ARRAY(lcp, orig->n_sa * 4);
    CHK_ARRAY(rle_data, orig->n_rle * 9);
    CHK_ARRAY(rle_offsets, orig->n_sa * 8);
    CHK_ARRAY(pos_rank, orig->n_pos * 8);
    CHK_ARRAY(pos_pac, orig->n_pos * 8);

#undef CHK
#undef CHK_ARRAY

    fprintf(fp, "# result: %d mismatch(es)\n", nbad);
    fclose(fp);
    fprintf(stderr, "[index-capture]   load verify: %d mismatch(es) → %s\n",
            nbad, nbad ? "FAIL" : "PASS");
}

/* ── 2-occurrence count on dense BWT ──────────────────────────────── */

static void capt_2occ(const capt_t *capt, uint64_t k, uint64_t l,
                      uint64_t tk[4], uint64_t tl[4])
{
    if (l > capt->n_sa) l = capt->n_sa;
    if (k > capt->n_sa) k = capt->n_sa;

    /* OCC at k: fetch block values + linear scan once for all 4 bases */
    uint64_t blk_k = k / CAPT_OCC_INTV, start_k = blk_k * CAPT_OCC_INTV;
    uint64_t oc_k[4];
    for (int c = 0; c < 4; c++) oc_k[c] = capt->occ[c][blk_k];
    for (uint64_t i = start_k; i < k; i++) {
        uint8_t b = capt->bwt[i];
        if (b < 4) oc_k[b]++;
    }

    /* OCC at l: same */
    uint64_t blk_l = l / CAPT_OCC_INTV, start_l = blk_l * CAPT_OCC_INTV;
    uint64_t oc_l[4];
    /* fast path: if same block as k, copy and extend */
    if (blk_l == blk_k) {
        for (int c = 0; c < 4; c++) oc_l[c] = oc_k[c];
        for (uint64_t i = k; i < l; i++) {
            uint8_t b = capt->bwt[i];
            if (b < 4) oc_l[b]++;
        }
    } else {
        for (int c = 0; c < 4; c++) oc_l[c] = capt->occ[c][blk_l];
        for (uint64_t i = start_l; i < l; i++) {
            uint8_t b = capt->bwt[i];
            if (b < 4) oc_l[b]++;
        }
    }

    for (int c = 0; c < 4; c++) {
        tk[c] = oc_l[c] - oc_k[c];
        tl[c] = oc_k[c];
    }
}

/* ── Sub-index interval init / extend (mirrors bwt_set_intv / bwt_extend) */

void capt_set_intv(const capt_t *capt, int c, bwtintv_t *ik)
{
    uint64_t ks_fwd = 0, ks_rc = 0;
    for (int b = 0; b < c; b++)
        ks_fwd += capt->end_cnt[b];
    for (int b = 0; b < 3 - c; b++)
        ks_rc += capt->end_cnt[b];
    ik->x[0] = capt->l2[c] + 1 + ks_fwd;
    ik->x[1] = capt->l2[3 - c] + 1 + ks_rc;
    ik->x[2] = capt->l2[c + 1] - capt->l2[c] + capt->end_cnt[c];
}

void capt_extend(const capt_t *capt, const bwtintv_t *ik,
                 bwtintv_t ok[4], int is_back)
{
    uint64_t tk[4], tl[4];
    uint64_t k = ik->x[!is_back] - 1;
    uint64_t l = k + ik->x[2];
    capt_2occ(capt, k, l, tk, tl);
    uint64_t ks[4];
    ks[0] = capt->end_cnt[0];
    ks[1] = capt->end_cnt[0] + capt->end_cnt[1];
    ks[2] = capt->end_cnt[0] + capt->end_cnt[1] + capt->end_cnt[2];
    ks[3] = capt->end_cnt[0] + capt->end_cnt[1] + capt->end_cnt[2] + capt->end_cnt[3];
    for (int i = 0; i < 4; i++) {
        ok[i].x[!is_back] = capt->l2[i] + 1 + tl[i] + ks[i];
        ok[i].x[2] = tk[i];
    }
    /* RC-strand rank: complement mapping A↔T(0↔3), C↔G(1↔2) */
    int sentinel_num = l - k - tk[0] - tk[1] - tk[2] - tk[3];
    ok[3].x[is_back] = ik->x[is_back] + sentinel_num;
    ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2];
    ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2];
    ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2];
}

/* ── Init capt_intv_t: set both genome and sub-index intervals ─────── */

void capt_smem1_init(capt_intv_t *ik, const capt_t *capt,
                     const bwt_t *g_bwt, int base)
{
    bwt_set_intv(g_bwt, base, *((bwtintv_t *)&ik->x));
    capt_set_intv(capt, base, &ik->sub);
    ik->on_genome = 0;
    ik->depth = 1;
    ik->info = 0;
}

/* ── Debug print ──────────────────────────────────────────────────── */

void capt_int_print(const capt_intv_t *p, const char *label)
{
    if (p->on_genome)
        fprintf(stderr, "  [%s] genome: k=%lld l=%lld s=%lld info=(hi=%u,lo=%u)\n",
                label, (long long)p->x[0], (long long)p->x[1], (long long)p->x[2],
                (unsigned)(p->info >> 32), (unsigned)(p->info & 0xffffffff));
    else
        fprintf(stderr, "  [%s] sub:    k=%lld l=%lld s=%lld info=(hi=%u,lo=%u)\n",
                label, (long long)p->sub.x[0], (long long)p->sub.x[1], (long long)p->sub.x[2],
                (unsigned)(p->info >> 32), (unsigned)(p->info & 0xffffffff));
}

/* ── capt_extend_bail: sub-index first, genome fallback via RLE ────── */

void capt_extend_bail(const capt_t *capt, const bwt_t *g_bwt,
                      const capt_intv_t *ik, int base, capt_intv_t *ok,
                      int is_back)
{
    ok->info = ik->info;
    ok->depth = ik->depth + 1;
    if (!ik->on_genome) {
        bwtintv_t saved = ik->sub;

        /* try sub-index extension */
        bwtintv_t ok_sub[4];
        capt_extend(capt, &ik->sub, ok_sub, is_back);

        if (ok_sub[base].x[2] > 0) {
            /* success: stay on sub-index */
            ok->x[0] = ik->x[0];
            ok->x[1] = ik->x[1];
            ok->x[2] = ik->x[2];
            ok->sub  = ok_sub[base];
            ok->on_genome = 0;
        } else if (saved.x[2] > 0) {
            /* fail: RLE translate → genome, extend on genome BWT */
            capt_intv_t tmp;
            tmp.on_genome = 0;
            tmp.depth = ik->depth;
            tmp.sub = saved;
            capt_translate_to_genome(capt, &tmp);
            if (tmp.x[2] > 0 && tmp.x[0] <= (uint64_t)(g_bwt->seq_len)) {
                bwtintv_t ok_gen[4];
                bwt_extend(g_bwt, (const bwtintv_t *)&tmp.x, ok_gen, is_back);
                ok->x[0] = ok_gen[base].x[0];
                ok->x[1] = ok_gen[base].x[1];
                ok->x[2] = ok_gen[base].x[2];
            } else {
                memset(ok, 0, sizeof(*ok));
            }
            ok->sub  = saved;
            ok->on_genome = 1;
        } else {
            memset(ok, 0, sizeof(*ok));
            ok->sub  = saved;
            ok->on_genome = 1;
        }
    } else {
        /* already on genome BWT */
        bwtintv_t ok_gen[4];
        bwt_extend(g_bwt, (const bwtintv_t *)&ik->x, ok_gen, is_back);
        ok->x[0] = ok_gen[base].x[0];
        ok->x[1] = ok_gen[base].x[1];
        ok->x[2] = ok_gen[base].x[2];
        ok->sub   = ik->sub;
        ok->on_genome = 1;
    }
}

/* ── Translate sub-SA interval to genome-SA in-place ───────────────── */

void capt_translate_to_genome(const capt_t *capt, capt_intv_t *ik)
{
    if (ik->on_genome || ik->sub.x[2] == 0) return;
    int d = ik->depth;
    /*fprintf(stderr, "[translate] depth=%d sub=(%lld,%lld,%lld)\n", d,
            (long long)ik->sub.x[0], (long long)ik->sub.x[1], (long long)ik->sub.x[2]);*/
    uint64_t gl = capt_rle_left(capt, ik->sub.x[0] - 1, d);
    uint64_t gr = capt_rle_right(capt,
        ik->sub.x[0] + ik->sub.x[2] - 2, d);
    uint64_t gl_inv = capt_rle_left(capt, ik->sub.x[1] - 1, d);
    /*fprintf(stderr, "[translate] gl=%llu gr=%llu gl_inv=%llu\n",
            (unsigned long long)gl, (unsigned long long)gr, (unsigned long long)gl_inv);*/
    ik->x[0] = gl;
    ik->x[1] = gl_inv;
    ik->x[2] = (gr >= gl) ? gr - gl + 1 : 0;
    ik->on_genome = 1;
}

/* ── Push genome-SA bi-interval ────────────────────────────────────── */

void capt_push_curr(bwtintv_v *curr, const capt_intv_t *ik, int depth,
                    const capt_t *capt)
{
    bwtintv_t t;
    if (!ik->on_genome) {
        capt_intv_t tmp = *ik;
        capt_translate_to_genome(capt, &tmp);
        t = *(bwtintv_t *)&tmp.x;
    } else {
        t = *(bwtintv_t *)&ik->x;
    }
    t.info = ik->info;
    kv_push(bwtintv_t, *curr, t);
}

/* ── Capture-aware SMEM ────────────────────────────────────────────── */

int cmp_capt_intv_t(const capt_intv_t* a,const capt_intv_t* b)
{
    if(!a->on_genome && b->on_genome){return 1;}
    if(a->on_genome && !b->on_genome){return -1;}
    if(a->on_genome){return (a->x[2])>(b->x[2])?1:((a->x[2])==(b->x[2])?0:-1);}
    else{return (a->sub.x[2])>(b->sub.x[2])?1:((a->sub.x[2])==(b->sub.x[2])?0:-1);}
}
int get_interval_len_capt(const capt_intv_t* a)
{
    return a->on_genome?a->x[2]:a->sub.x[2];
}

int capt_smem1(const capt_t *capt, const bwt_t *g_bwt,
               int len, const uint8_t *q, int x, int min_intv,
               bwtintv_v *mem, capt_intv_v *tmpvec[2])
{
    int max_intv=1;
    int i, j, c, ret;
	capt_intv_t ik, ok;
	capt_intv_v a[2], *prev, *curr, *swap;

	mem->n = 0;
	if (q[x] > 3) return x + 1;
	if (min_intv < 1) min_intv = 1; // the interval size should be at least 1
	kv_init(a[0]); kv_init(a[1]);
	prev = tmpvec && tmpvec[0]? tmpvec[0] : &a[0]; // use the temporary vector if provided
	curr = tmpvec && tmpvec[1]? tmpvec[1] : &a[1];
	capt_smem1_init(&ik, capt, g_bwt, q[x]); // the initial interval of a single base
	
    ik.info = x + 1;
    for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
		if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			capt_extend_bail(capt, g_bwt, &ik, c, &ok, 0);
			if (cmp_capt_intv_t(&ok,&ik)!=0) { // change of the interval size
				kv_push(capt_intv_t, *curr, ik);
                //capt_int_print(&ik, "ik");
				if (get_interval_len_capt(&ok) < min_intv) break; // the interval size is too small to be extended further
			}
			ik = ok; ik.info = i + 1;
		} else { // an ambiguous base
			kv_push(capt_intv_t, *curr, ik);
			break; // always terminate extension at an ambiguous base; in this case, i<len always stands
		}
	}
    
    if (i == len) kv_push(capt_intv_t, *curr, ik); // push the last interval if we reach the end
	//capt_int_print(&ik, "ik_last");
    /* reverse curr: smaller intervals (longer matches) first */
    for (int _r = 0; _r < (int)curr->n / 2; _r++) {
        capt_intv_t _t = curr->a[_r];
        curr->a[_r] = curr->a[curr->n - 1 - _r];
        curr->a[curr->n - 1 - _r] = _t;
    }
	ret = curr->a[0].info; // this will be the returned value
	swap = curr; curr = prev; prev = swap;
	/*fprintf(stderr, "[capt_smem1] prev (forward hits): n=%zu\n", prev->n);
	for (int _i = 0; _i < (int)prev->n; _i++)
		capt_int_print(&prev->a[_i], "prev");*/
	for (i = x - 1; i >= -1; --i) { // backward search for MEMs
		c = i < 0? -1 : q[i] < 4? q[i] : -1; // c==-1 if i<0 or q[i] is an ambiguous base
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			capt_intv_t *p = &prev->a[j];
			if (c >= 0 ) capt_extend_bail(capt, g_bwt, p, c, &ok, 1);//bwt_extend(bwt, p, ok, 1);
			if (c < 0 || get_interval_len_capt(&ok) < min_intv) { // keep the hit if reaching the beginning or an ambiguous base or the intv is small enough
				if (curr->n == 0) { // test curr->n>0 to make sure there are no longer matches
					if (mem->n == 0 || i + 1 < mem->a[mem->n-1].info>>32) { // skip contained matches
						ik = *p; ik.info |= (uint64_t)(i + 1)<<32;
						capt_push_curr(mem, &ik, (unsigned)(p->info & 0xffffffff)-i, capt);                        
					}
				} // otherwise the match is contained in another longer match
			} else if (curr->n == 0 || cmp_capt_intv_t(&ok, &(curr->a[curr->n-1]))!=0) 
            {
				ok.info = p->info;
				kv_push(capt_intv_t, *curr, ok);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}
	bwt_reverse_intvs(mem); // s.t. sorted by the start coordinate

	if (tmpvec == 0 || tmpvec[0] == 0) free(a[0].a);
	if (tmpvec == 0 || tmpvec[1] == 0) free(a[1].a);
	return ret;
}
/* ── Position table lookup (inverse of bwt_sa) ─────────────────────── */

int64_t capt_lookup_pos(const capt_t *capt, uint64_t rank)
{
    uint64_t lo = 0, hi = capt->n_pos;
    while (lo < hi) {
        uint64_t mid = (lo + hi) >> 1;
        if (capt->pos_rank[mid] < rank) lo = mid + 1;
        else                             hi = mid;
    }
    if (lo < capt->n_pos && capt->pos_rank[lo] == rank)
        return capt->pos_pac[lo];
    return -1;
}
