/*
 * bwamem_capt.c — Capture-aware alignment core.
 *
 * capt_mem_align1_core mirrors mem_align1_core but uses capt_smem1
 * for SMEM seeding instead of bwt_smem1.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwa.h"
#include "bwt.h"
#include "bntseq.h"
#include "bwamem.h"
#include "capt.h"

extern unsigned char nst_nt4_table[256];

mem_alnreg_v capt_mem_align1_core(const mem_opt_t *opt,
                                   const capt_t *capt,
                                   const bwt_t *g_bwt,
                                   const bntseq_t *bns,
                                   const uint8_t *pac,
                                   int l_seq, char *seq, void *buf)
{
    int i;
    mem_chain_v chn;
    mem_alnreg_v regs;

    for (i = 0; i < l_seq; ++i)
        seq[i] = seq[i] < 4 ? seq[i] : nst_nt4_table[(int)seq[i]];

    /* capture-aware SMEM */
    chn = mem_chain(opt, g_bwt, bns, l_seq, (uint8_t *)seq, buf);
    chn.n = mem_chain_flt(opt, chn.n, chn.a);
    mem_flt_chained_seeds(opt, bns, pac, l_seq, (uint8_t *)seq, chn.n, chn.a);
    if (bwa_verbose >= 4)
        mem_print_chain(bns, &chn);

    kv_init(regs);
    for (i = 0; i < chn.n; ++i) {
        mem_chain_t *p = &chn.a[i];
        if (bwa_verbose >= 4)
            err_printf("* ---> Processing chain(%d) <---\n", i);
        mem_chain2aln(opt, bns, pac, l_seq, (uint8_t *)seq, p, &regs);
        free(chn.a[i].seeds);
    }
    free(chn.a);
    regs.n = mem_sort_dedup_patch(opt, bns, pac, (uint8_t *)seq, regs.n, regs.a);
    if (bwa_verbose >= 4) {
        err_printf("* %ld chains remain after removing duplicated chains\n", regs.n);
        for (i = 0; i < regs.n; ++i) {
            mem_alnreg_t *p = &regs.a[i];
            printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
        }
    }
    for (i = 0; i < regs.n; ++i) {
        mem_alnreg_t *p = &regs.a[i];
        if (p->rid >= 0 && bns->anns[p->rid].is_alt)
            p->is_alt = 1;
    }
    return regs;
}