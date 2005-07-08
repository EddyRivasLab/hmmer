#include "config.h"
#include "squidconf.h"

#include "structs.h"
#include "funcs.h"

#if (defined __GNUC__) && (defined __APPLE__)
#include <ppc_intrinsics.h>
#endif

/* Echo a message about Altivec being used once: 
 *   (1) Confirm that the Altivec kernel is indeed being used
 *   (2) Trace bugs to the Altivec kernel instead of the vanilla code
 */
/*
 * Note:  This was originally in fast_algorithms.c, inside of an
 *        ifdef ALTIVEC preprocessor block, but it moved it here, with
 *        all of the other altivec-specific functions. - CRS 23 June 2005
 */
static int first_altivec = 1;
void
AltivecMessage()
{
  printf("Using Altivec/VMX code. Bug reports to Erik Lindahl <lindahl@sbc.su.se>.\n");
  first_altivec = 0;
}

/* Function: Viterbi()
 *
 * Note:     This was originally defined as P7Viterbi() in fast_algorithms.c,
 *           inside of an #ifdef ALTIVEC preprocessor block.  It is the 
 *           altivec specific impelmention of Viterbi, so I moved it here with
 *           the rest of the altivec specific (non-structure-related) functions.
 *             - CRS 23 June 2005
 * 
 * Purpose:  The Viterbi dynamic programming algorithm; Altivec implementation
 *           by Erik Lindahl, Stanford University, 2002.
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           mx     - DP matrix (may get grown here)
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
/* This first version of P7Viterbi has been accelerated with Altivec vectorization.
 * On Apple hardware, it is up to a factor 10 faster than the old non-altivec version.
 */
typedef union 
{
    vector signed int v;
    int i[4];
} ivector;

void 
printivec(vector signed int z)
{
    ivector q;
    q.v=z;
    printf("%d  %d  %d  %d\n",q.i[0],q.i[1],q.i[2],q.i[3]);
}


float
Viterbi(unsigned char *dsq, int L, struct plan7_s *hmm, struct dpmatrix_s *mx, struct p7trace_s **ret_tr) 
{
    struct p7trace_s  *tr;
    int **xmx;
    int **mmx;
    int **imx;
    int **dmx;
    int  *mmxi,*xmxi,*imxi,*dmxi;
    int  *lmmxi,*lxmxi,*limxi,*ldmxi;
    int  *p_tmm,*p_tim,*p_tdm,*p_bsc,*p_msc,*p_tmd,*p_tdd,*p_tmi,*p_tii,*p_isc,*p_esc;
    int   i,n,k;
    int   sc;
    /* gcc and motorola use different syntax for initializing vectors,
     * so we use a dummy variable instead to get an adress to load from...
     */   
    int t_lowscore = -INFTY;
    
    /* vector variables. We avoid the stupid gcc spill/fill code by
     * limiting ourselves to 32 generic variables and making the all registers.
     * (This reuse is the reason for the generic variable names).
     */
    vector signed int v_lowscore;
    vector signed int max_mmxesc;
    vector signed int v_xmb;
    vector unsigned int mask1;
    vector unsigned int mask2;
    vector unsigned int mask3;
    vector unsigned int mask4;
    vector signed int v_lmmx1;
    vector signed int v_lmmx2;
    vector signed int v_limx1;
    vector signed int v_limx2;
    vector signed int v_save_lmmx;
    vector signed int v_save_ldmx;
    vector signed int v_save_limx;
    vector signed int v_save_mmx;
    vector signed int v_save_dmx;
    vector signed int v1;
    vector signed int v2;
    vector signed int v3;
    vector signed int v4;
    vector signed int v5;
    vector signed int v6;
    vector signed int v7;
    vector signed int v8;
    vector signed int v9;
    vector signed int v10;
    vector signed int v11;
    vector signed int v12;
    vector signed int v13;
    vector signed int v14;
    vector signed int v15;

    if(first_altivec)
      AltivecMessage();

    /* load (-infinity) to all four elements in v_lowscore */
    v_lowscore      = vec_lde(0, &t_lowscore );
    mask1           = (vector unsigned int)vec_lvsl(0,&t_lowscore);
    v_lowscore      = vec_perm(v_lowscore,v_lowscore,(vector unsigned char)mask1);
    v_lowscore      = vec_splat(v_lowscore,0);
  
    v1 = vec_splat_s32(-1);
    v2 = vec_splat_s32(0);
    mask1 = (vector unsigned int)vec_sld(v1,v2,12); /* FF in first pos, rest. are 00 */
    mask2 = vec_sld(mask1,mask1,12);
    mask3 = vec_sld(mask1,mask1,8);
    mask4 = vec_sld(mask1,mask1,4);
    
    /* Make sure our DP matrix has 0..L rows, 0..M columns; grow it if needed. */
    ResizeDPMatrix(mx, L, hmm->M, &xmx, &mmx, &imx, &dmx);
    
    /* Initialization of the zero row. */
    xmx[0][XMN] = 0;		                     /* S->N, p=1            */
    xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
    xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  
    mmxi=mmx[0];
    imxi=imx[0];
    dmxi=dmx[0];
    xmxi=xmx[0];
    
    for (n = 0; n  < 5+hmm->M; n+=4) 
    {
        vec_st(v_lowscore, n*4, mmxi);
        vec_st(v_lowscore, n*4, imxi);
        vec_st(v_lowscore, n*4, dmxi);
    }

    /* Fill data beyound M with -INFTY, so we can take the maximum including
     * elements with k>M.
     */
    for(k=1+hmm->M;k<=(3+hmm->M);k++) 
    {
        hmm->esc[k]=-INFTY;
        hmm->bsc[k]=-INFTY;
        for(i=0;i<7;i++)
            hmm->tsc[i][k]=-INFTY;
        for(i=0;i<MAXCODE;i++)
        {
            hmm->msc[i][k]=-INFTY;
            hmm->isc[i][k]=-INFTY;
        }
    }
  
    /* Recursion. Done as a pull
     * Note some slightly wasteful boundary conditions:  
     * tsc[0] = -INFTY for all eight transitions (no node 0)
     * D_M and I_M are wastefully calculated (they don't exist)
     */

    for (i = 1; i <= L; i++) 
    {
        /* pointers to last (i-1) row */
        lmmxi=mmxi;
        limxi=imxi;
        ldmxi=dmxi;
        lxmxi=xmxi;
        
        /* get pointers to this row */
        mmxi=mmx[i];
        imxi=imx[i];
        dmxi=dmx[i];
        xmxi=xmx[i];
        
        /* Set everything that doesnt depend on k here */    
        
        /* load and splat (spread to all elements) XMX[i-1][XMB] */
        v13   = vec_lde(0,&(xmx[i-1][XMB]));
        v14   = (vector signed int)vec_lvsl(0,&(xmx[i-1][XMB]));
        v13   = vec_perm(v13,v13,(vector unsigned char)v14);
        v_xmb = vec_splat(v13,0);
        p_tmm = hmm->tsc[TMM];
        p_tim = hmm->tsc[TIM];
        p_tdm = hmm->tsc[TDM];
        p_bsc = hmm->bsc;
        k = dsq[i];
        p_msc = hmm->msc[k];
        p_isc = hmm->isc[k];
        p_tmd = hmm->tsc[TMD];
        p_tdd = hmm->tsc[TDD];
        p_tmi = hmm->tsc[TMI];
        p_tii = hmm->tsc[TII];
        p_esc = hmm->esc;
        max_mmxesc = v_lowscore;
        
        /* the 0 element of vectors are aligned 12 bytes up from the 16-byte boundary,
         * so we simply write the entire 16 bytes before the 16 byte boundary.
         */
        vec_st(v_lowscore,0,mmxi);
        vec_st(v_lowscore,0,imxi);
        vec_st(v_lowscore,0,dmxi);
        
        /* Load the first (i.e. 'previous') vector on last row for mmx,imx,dmx,
         * and load the first vector of dmx and mmx on this row.
         */
        v_save_lmmx = vec_ld(-12, lmmxi);
        v_save_limx = vec_ld(-12, limxi);
        v_save_ldmx = vec_ld(-12, ldmxi);
        v_save_mmx = vec_ld(-12, mmxi);
        v_save_dmx = vec_ld(-12, dmxi);
        
        /* we have allocated extra memory, so it is perfectly OK
         * to do the calculations for a couple of extra cells where
         * k>hmm->M. These cells just wont be used.
         */
        for (n = 4, k=1; k < (hmm->M-3) ; n+=32, k+=8) 
        {
            /* match state */
            
            /* 1: check which of mmx[i-1][k-1]+TMM and imx[i-1][k-1]+TIM is better,
             * but do it for 8 elements in parallel.
             * Since we are comparing with data on the previous row but one position
             * earlier, we have to shift stuff. Load two new vectors each round,
             * and use the saved from last round. 
             */
            /* load mmx data */
            v_lmmx1 = vec_ld(n,    lmmxi);
            v_lmmx2 = vec_ld(n+16, lmmxi);

            /* load imx data */
            v_limx1 = vec_ld(n,    limxi);
            v_limx2 = vec_ld(n+16, limxi);
            
            v5    = vec_ld(n,    ldmxi);	/* Load dmx data */
            v10   = vec_ld(n+16, ldmxi);	
            
            /* shift mmx, imx & dmx data */
            v1    = vec_sld(v_save_lmmx,v_lmmx1,12);
            v3    = vec_sld(v_save_limx,v_limx1,12);
            v9    = vec_sld(v_save_ldmx,v5,12);
            
            /* shift mmx, imx & dmx data */
            v2    = vec_sld(v_lmmx1,v_lmmx2,12);
            v4    = vec_sld(v_limx1,v_limx2,12);
            v_save_ldmx = v10;
            v10   = vec_sld(v5,v10,12);	
            
            v_save_lmmx = v_lmmx2;
            v_save_limx = v_limx2;
            
            /* v1,v2 now contains 8 element with mmx[i-1][k-1], 
             * v3,v4 contain 8 elements with imx[i-1][k-1],
             * and v9,v10 contain 8 elements with dmx[i-1][k-1].
             */
            /* load TMM, TIM & TDM entries from the HMM - these are aligned in memory */
            v5    = vec_ld(n-4, p_tmm);
            v6    = vec_ld(n+12, p_tmm);
            v7    = vec_ld(n-4, p_tim);
            v8    = vec_ld(n+12, p_tim);
            v11   = vec_ld(n-4, p_tdm);
            v12   = vec_ld(n+12, p_tdm);

            /* load bsc[k] */
            v14   = vec_ld(n, p_bsc);
            v15   = vec_ld(n+16, p_bsc);
            
            /* calc mmx+TMM, imx+TIM, dmx+TDM and XMX+bsc with saturated arithmetic, so 
             * we don't loop if we add the large negative numbers used for -infinity.
             */
            v1    = vec_adds(v1,v5);
            v2    = vec_adds(v2,v6);
            v3    = vec_adds(v3,v7);
            v4    = vec_adds(v4,v8);
            v9    = vec_adds(v9,v11);
            v10   = vec_adds(v10,v12);
            v14   = vec_adds(v14,v_xmb);
            v15   = vec_adds(v15,v_xmb);
            /* Select max of mmx+TMM and imx+TIM in each element */
            v1    = vec_max(v1,v3);
            v2    = vec_max(v2,v4);
            /* select max of dmx+TDM and XMX+bsc */
            v9    = vec_max(v9,v14);
            v10   = vec_max(v10,v15);	
            /* select max of the four alternatives */
            v1    = vec_max(v1,v9);
            v2    = vec_max(v2,v10);
            /* v1,v2 now contain the max values for the new mmx;
             * check if we should add msc.
             */
            
            v3    = vec_ld(n,    p_msc);
            v4    = vec_ld(n+16, p_msc);
            
            v5    = (vector signed int)vec_cmpgt(v3,v_lowscore);
            v6    = (vector signed int)vec_cmpgt(v4,v_lowscore);
            
            /* load esc[k] */
            v9    = vec_ld(n, p_esc);
            v10   = vec_ld(n+16, p_esc);
            
            v1    = vec_adds(v1,v3);
            v2    = vec_adds(v2,v4);
            v1    = vec_sel(v3,v1,(vector unsigned int)v5);
            v2    = vec_sel(v4,v2,(vector unsigned int)v6);      
            
            /* have final values for mmx on this row in v1,v2 - store it */
            vec_st(v1, n,    mmxi);
            vec_st(v2, n+16, mmxi);
            
            v9    = vec_adds(v9,v1);
            v10   = vec_adds(v10,v2);
            v9    = vec_max(v9,v10);
            max_mmxesc = vec_max(v9,max_mmxesc);
            
            /* OK. We finished the match state. Normally we could now first
             * do the delete and then the insertion. The problem is that the
             * delete is a pain in the ass to vectorize, since each element
             * depends on the previous one in the vector. This means I have
             * to write this relatively simple operation as eight independent
             * ones, just as we would have done in a non-vectorized code.
             * Since this is independent of the insertion state changes, I
             * try to hide the latencies by doing the delete and insert
             * calculations in parallel.
             * To make things easier I add 'del' in a comment on each
             * line for calculations that are on the delete state, and 'ins'
             * for the calculations for the insertion. Hang on...
             */
	    
            /* We already have the match data on this row from the previous
             * iteration in v_save_mmx. Rotate it so the element that used to be
             * is pos 4 last iteration is in position 1, and pos 2-4 contain
             * the the first three elements of mmx this iteration.
             * And do the same type of rotation for v1/v2...
             */

            v4 = vec_sld(v1,v2,12); /* del */
            v3 = vec_sld(v_save_mmx,v1,12); /* del */
            v_save_mmx = v2; /* Save for next iteration */
            
            /* Rotate last dmx data so we have the fourth element in pos 1. */
            v_save_dmx = vec_sld(v_save_dmx,v_save_dmx,12); /* del */
            
            /* load TMD & TDD data */
            v5   = vec_ld(n-4, p_tmd); /* del */
            v6   = vec_ld(n+12, p_tmd); /* del */
            v7   = vec_ld(n-4, p_tdd); /* del */
            v8   = vec_ld(n+12, p_tdd); /* del */
            
            /* calculate mmx+TMD */
            v3   = vec_adds(v3,v5); /* del */
            v4   = vec_adds(v4,v6); /* del */
            
            /* Start the ugly stuff. We have rotated last dmx data. Add TDD to
             * it and compare with v3/v4 (data from mmx+TMD alternative), but
             * we only compare & replace the first position, so we use a mask!
             */

            /* First position: Add TDD to v_save_dmx */
            v_save_dmx = vec_adds(v_save_dmx,v7); /* del */
            /* Select max of this and mmx+TMD and save it temporary to v_save_dmx */
            v_save_dmx = vec_max(v_save_dmx,v3); /* del */
            /* use a mask to select only the first element from v_save_dmx, rest from v3. */
            v3     = vec_sel(v3,v_save_dmx,mask1); /* del */
            
            /* Start doing the insertion calculations in parallel.
             * Load the TMI data.
             */
            v9    = vec_ld(n   , p_tmi); /* ins */
            v10   = vec_ld(n+16, p_tmi); /* ins */
            /* Deletion:
             * Now we have an accurate pos 1 in v3. continue with pos2 the same
             * way. Rotate to a temporary array, add to TDD, compare and use a mask
             * to write back to only position 2.
             */
            v_save_dmx = vec_sld(v3,v3,12); /* del */
            v_save_dmx = vec_adds(v_save_dmx,v7); /* del */
            v_save_dmx = vec_max(v_save_dmx,v3); /* del */
            v3     = vec_sel(v3,v_save_dmx,mask2); /* del */
            
            /* More insertion stuff - load TII data */
            v11   = vec_ld(n   , p_tii); /* ins */
            v12   = vec_ld(n+16, p_tii); /* ins */
            
            /* Deletion, position 3... */
            v_save_dmx = vec_sld(v3,v3,12); /* del */
            v_save_dmx = vec_adds(v_save_dmx,v7); /* del */
            v_save_dmx = vec_max(v_save_dmx,v3); /* del */
            v3     = vec_sel(v3,v_save_dmx,mask3); /* del */
            
            /* insertion stuff: calculate mmx+TMI */
            v9     = vec_adds(v_lmmx1,v9); /* ins */
            v10    = vec_adds(v_lmmx2,v10); /* ins */
            
            /* Deletion, position 4 */
            v_save_dmx = vec_sld(v3,v3,12); /* del */
            v_save_dmx = vec_adds(v_save_dmx,v7); /* del */
            v_save_dmx = vec_max(v_save_dmx,v3); /* del */
            v3     = vec_sel(v3,v_save_dmx,mask4); /* del */
            
            /* insertion stuff: calculate imx+TII */
            v11    = vec_adds(v_limx1,v11); /* ins */
            v12    = vec_adds(v_limx2,v12); /* ins */
            
            /* That was the first deletion vector, but we are unrolling.
             * The next step is position '5', i.e. the first in vector 2.
             * This one depends on the last position of vector 1, which we just finished.
             */
            v_save_dmx = vec_sld(v3,v3,12); /* del */
            v_save_dmx = vec_adds(v_save_dmx,v8); /* del */
            v_save_dmx = vec_max(v_save_dmx,v4); /* del */
            v4     = vec_sel(v4,v_save_dmx,mask1); /* del */
            
            /* insertion stuff: select max of mmx+TMI and imx+TII */
            v9     = vec_max(v9,v11); /* ins */
            v10    = vec_max(v10,v12); /* ins */
            /* insertion stuff: load data from hmm->isc[tmpidx] */
            v11    = vec_ld(n, p_isc); /* ins */
            v12    = vec_ld(n+16, p_isc); /* ins */
            
            /* position 6 (2 in vector 2) */
            v_save_dmx = vec_sld(v4,v4,12); /* del */
            v_save_dmx = vec_adds(v_save_dmx,v8); /* del */
            v_save_dmx = vec_max(v_save_dmx,v4); /* del */
            v4     = vec_sel(v4,v_save_dmx,mask2); /* del */
            
            /* insertion: compare max of mmx+TMI, imx+TII with hmm->isc */
            v13    = (vector signed int)vec_cmpgt(v11,v_lowscore);
            v14    = (vector signed int)vec_cmpgt(v12,v_lowscore);
            
            /* position 7 (3 in vector 2) */
            v_save_dmx = vec_sld(v4,v4,12); /* del */
            v_save_dmx = vec_adds(v_save_dmx,v8); /* del */
            v_save_dmx = vec_max(v_save_dmx,v4); /* del */
            v4     = vec_sel(v4,v_save_dmx,mask3); /* del */
            
            v9     = vec_adds(v9,v11);
            v10    = vec_adds(v10,v12);
            v9     = vec_sel(v11,v9,(vector unsigned int)v13);
            v10    = vec_sel(v12,v10,(vector unsigned int)v14);
            
            /* position 8 (4 in vector 2) */
            v_save_dmx = vec_sld(v4,v4,12); /* del */
            v_save_dmx = vec_adds(v_save_dmx,v8); /* del */
            v_save_dmx = vec_max(v_save_dmx,v4); /* del */
            v_save_dmx = vec_sel(v4,v_save_dmx,mask4); /* del */
            
            /* Puh! that was the deletion... v3/v_save_dmx now contain the updated dmx data;
             * save it to memory. (v_save_dmx will be used next iteration)
             */
            vec_st(v3,        n, dmxi); /* del */
            vec_st(v_save_dmx, n+16, dmxi); /* del */
            
            /* save insertion too */
            vec_st(v9, n, imxi);
            vec_st(v10,n+16,imxi);
        }
        /* odd loop */
        if(k< (1+hmm->M)) 
        {
            /* match state */
            /* 1: check which of mmx[i-1][k-1]+TMM and imx[i-1][k-1]+TIM is better,
            * but do it for 4 elements in parallel.
            * Since we are comparing with data on the previous row but one position
            * earlier, we have to shift stuff. Load two new vectors each round,
            * and use the saved from last round.
            */
            /* load mmx data */
            v_lmmx1 = vec_ld(n,    lmmxi);
            
            /* load imx data */
            v_limx1 = vec_ld(n,    limxi);
            
            v5    = vec_ld(n,    ldmxi);	/* Load dmx data */
            
            /* shift mmx, imx & dmx data */
            v1    = vec_sld(v_save_lmmx,v_lmmx1,12);
            v3    = vec_sld(v_save_limx,v_limx1,12);
            v9    = vec_sld(v_save_ldmx,v5,12);
            
            /* v1,v2 now contains 8 element with mmx[i-1][k-1],
             * v3,v4 contain 8 elements with imx[i-1][k-1],
             * and v9,v10 contain 8 elements with dmx[i-1][k-1].
             */
            /* load TMM, TIM & TDM entries from the HMM - these are aligned in memory */
            v5    = vec_ld(n-4, p_tmm);
            v7    = vec_ld(n-4, p_tim);
            v11   = vec_ld(n-4, p_tdm);
            /* load bsc[k] */
            v14   = vec_ld(n, p_bsc);
            
            /* calc mmx+TMM, imx+TIM, dmx+TDM and XMX+bsc with saturated arithmetic, so
             * we don't loop if we add the large negative numbers used for -infinity.
             */
            v1    = vec_adds(v1,v5);
            v3    = vec_adds(v3,v7);
            v9    = vec_adds(v9,v11);
            v14   = vec_adds(v14,v_xmb);
            /* Select max of mmx+TMM and imx+TIM in each element */
            v1    = vec_max(v1,v3);
            /* select max of dmx+TDM and XMX+bsc */
            v9    = vec_max(v9,v14);
            /* select max of the four alternatives */
            v1    = vec_max(v1,v9);
            /* v1,v2 now contain the max values for the new mmx;
             * check if we should add msc.
             */
            
            v3    = vec_ld(n,    p_msc);
            
            v5    = (vector signed int)vec_cmpgt(v3,v_lowscore);
        
            /* load esc[k] */
            v9    = vec_ld(n, p_esc);
            
            v1    = vec_adds(v1,v3);
            v1    = vec_sel(v3,v1,(vector unsigned int)v5);
            
            /* have final values for mmx on this row in v1,v2 - store it */
            vec_st(v1, n,    mmxi);
            
            v9    = vec_adds(v9,v1);
            max_mmxesc = vec_max(v9,max_mmxesc);
            
            /* OK. We finished the match state. Normally we could now first
             * do the delete and then the insertion. The problem is that the
             * delete is a pain in the ass to vectorize, since each element
             * depends on the previous one in the vector. This means I have
             * to write this relatively simple operation as eight independent
             * ones, just as we would have done in a non-vectorized code.
             * Since this is independent of the insertion state changes, I
             * try to hide the latencies by doing the delete and insert
             * calculations in parallel.
             * To make things easier I add 'del' in a comment on each
             * line for calculations that are on the delete state, and 'ins'
             * for the calculations for the insertion. Hang on...
             */
            
            /* We already have the match data on this row from the previous
             * iteration in v_save_mmx. Rotate it so the element that used to be
             * is pos 4 last iteration is in position 1, and pos 2-4 contain
             * the the first three elements of mmx this iteration.
             * And do the same type of rotation for v1/v2...
             */
        
            v3 = vec_sld(v_save_mmx,v1,12); /* del */
            
            /* Rotate last dmx data so we have the fourth element in pos 1. */
            v_save_dmx = vec_sld(v_save_dmx,v_save_dmx,12); /* del */
            
            /* load TMD & TDD data */
            v5   = vec_ld(n-4, p_tmd); /* del */
            v7   = vec_ld(n-4, p_tdd); /* del */
            
            /* calculate mmx+TMD */
            v3   = vec_adds(v3,v5); /* del */
            
            /* Start the ugly stuff. We have rotated last dmx data. Add TDD to
             * it and compare with v3/v4 (data from mmx+TMD alternative), but
             * we only compare & replace the first position, so we use a mask!
             */
            
            /* First position: Add TDD to v_save_dmx */
            v_save_dmx = vec_adds(v_save_dmx,v7); /* del */
            /* Select max of this and mmx+TMD and save it temporary to v_save_dmx */
            v_save_dmx = vec_max(v_save_dmx,v3); /* del */
            /* use a mask to select only the first element from v_save_dmx, rest from v3. */
            v3     = vec_sel(v3,v_save_dmx,mask1); /* del */
            
            /* Start doing the insertion calculations in parallel.
             * Load the TMI data.
             */
            v9    = vec_ld(n   , p_tmi); /* ins */
            /* Deletion:
             * Now we have an accurate pos 1 in v3. continue with pos2 the same
             * way. Rotate to a temporary array, add to TDD, compare and use a mask
             * to write back to only position 2.
             */
            v_save_dmx = vec_sld(v3,v3,12); /* del */
            v_save_dmx = vec_adds(v_save_dmx,v7); /* del */
            v_save_dmx = vec_max(v_save_dmx,v3); /* del */
            v3     = vec_sel(v3,v_save_dmx,mask2); /* del */
            
            /* More insertion stuff - load TII data */
            v11   = vec_ld(n   , p_tii); /* ins */
            
            /* Deletion, position 3... */
            v_save_dmx = vec_sld(v3,v3,12); /* del */
            v_save_dmx = vec_adds(v_save_dmx,v7); /* del */
            v_save_dmx = vec_max(v_save_dmx,v3); /* del */
            v3     = vec_sel(v3,v_save_dmx,mask3); /* del */
            
            /* insertion stuff: calculate mmx+TMI */
            v9     = vec_adds(v_lmmx1,v9); /* ins */
            
            /* Deletion, position 4 */
            v_save_dmx = vec_sld(v3,v3,12); /* del */
            v_save_dmx = vec_adds(v_save_dmx,v7); /* del */
            v_save_dmx = vec_max(v_save_dmx,v3); /* del */
            v3     = vec_sel(v3,v_save_dmx,mask4); /* del */
            
            /* insertion stuff: calculate imx+TII */
            v11    = vec_adds(v_limx1,v11); /* ins */
            
            /* insertion stuff: select max of mmx+TMI and imx+TII */
            v9     = vec_max(v9,v11); /* ins */
            /* insertion stuff: load data from hmm->isc[tmpidx] */
            v11    = vec_ld(n, p_isc); /* ins */
            
            /* insertion: compare max of mmx+TMI, imx+TII with hmm->isc */
            v13    = (vector signed int)vec_cmpgt(v11,v_lowscore);
            
            v9     = vec_adds(v9,v11);
            v9     = vec_sel(v11,v9,(vector unsigned int)v13);
            
            /* Puh! that was the deletion... v3/v_save_dmx now contain the updated dmx data;
             * save it to memory. (v_save_dmx will be used next iteration)
             */
            vec_st(v3,        n, dmxi); /* del */
            
            /* save insertion too */
            vec_st(v9, n, imxi);
        }
        /* end of k loops */
        
        /* Now the special states. Order is important here.
         * remember, C and J emissions are zero score by definition,
         */
        /* N state */
        xmx[i][XMN] = -INFTY;
        if ((sc = xmx[i-1][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
            xmx[i][XMN] = sc;
        
        /* E state */
        v2 = vec_sld(max_mmxesc,max_mmxesc,8);
        v2 = vec_max(v2,max_mmxesc);
        v1 = vec_sld(v2,v2,4);
        v1 = vec_max(v1,v2);
        vec_ste(v1,XME*4,xmxi);
        
        /* J state */
        xmxi[XMJ] = -INFTY;
        if ((sc = lxmxi[XMJ] + hmm->xsc[XTJ][LOOP]) > -INFTY)
            xmxi[XMJ] = sc;
        if ((sc = xmxi[XME]   + hmm->xsc[XTE][LOOP]) > xmxi[XMJ])
            xmxi[XMJ] = sc;
        
        /* B state */
        xmxi[XMB] = -INFTY;
        if ((sc = xmxi[XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
            xmxi[XMB] = sc;
        if ((sc = xmxi[XMJ] + hmm->xsc[XTJ][MOVE]) > xmxi[XMB])
            xmxi[XMB] = sc;
        
        /* C state */
        xmxi[XMC] = -INFTY;
        if ((sc = lxmxi[XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
            xmxi[XMC] = sc;
        if ((sc = xmxi[XME] + hmm->xsc[XTE][MOVE]) > xmxi[XMC])
            xmxi[XMC] = sc;
    }
    /* T state (not stored) */
    sc = xmx[L][XMC] + hmm->xsc[XTC][MOVE];
    
    if (ret_tr != NULL) 
    {
        P7ViterbiTrace(hmm, dsq, L, mx, &tr);
        *ret_tr = tr; 
    }
    
    /* Note (Lindahl): we do NOT free the dpmatrix here anymore - the code was 
     * spending 30% of the runtime allocating/freeing memory.
     * Provide a pointer to a dpmatrix_s structure to this routine,
     * and we try to reuse it. After the final call to P7Viterbi,
     * free it with FreeDPMatrix.
     */
    return Scorify(sc);		/* the total Viterbi score. */
}

/* Function: P7ViterbiNoTrace()
 *
 * Note:     This was originally defined in fast_algorithms.c, inside of an
 *           #ifdef ALTIVED preprocessor blcok, but I moved it here so that
 *           it is with the (non-structure-related) altivec-specific 
 *           functions.  - CRS 23 June 2005
 * 
 * Purpose:  The Viterbi dynamic programming algorithm, but a version
 *           that does not store the dynamic programming matrix, and thus
 *           doesn't return the trace - just the score.
 *           This boost the Altivec performance significantly, since 
 *           the tuned code otherwise would be memory bandwidth limited.
 *           Altivec implementation by Erik Lindahl, Stanford University, 2004.
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           mx     - DP matrix (may get grown here)
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
float
P7ViterbiNoTrace(unsigned char *dsq, int L, struct plan7_s *hmm, struct dpmatrix_s *mx) 
{
    struct p7trace_s *   tr;
    int *                workp;
    int                  i,n,nn,k,M;
    int                  sc;
    int *                tpmm;
    int *                tpim;
    int *                tpdm;
    int *                tpmd;
    int *                tpdd;
    int *                tpmi;
    int *                tpii;
    int *                bp;
    int *                ep;
    int *                ms;
    int *                is;
    int                  xmx_XMB;
    int                  xmx_XME;
    int                  xmx_XMJ;
    int                  xmx_XMC;
    int                  xmx_XMN;
    
    vector signed int m_infty = (vector signed int)(0x80000000,0x80000000,0x80000000,0x80000000);

    vector signed int v_xmb;
    vector signed int v_xme;
    vector signed int v_bp;
    vector signed int v_ep;
    vector signed int v_ms;
    vector signed int v_is;
    vector signed int v_mpp_1;
    vector signed int v_mpp_2;
    vector signed int v_mpp_rot;
    vector signed int v_dpp_1;
    vector signed int v_dpp_2;
    vector signed int v_dpp_rot;
    vector signed int v_ip_1;
    vector signed int v_ip_2;
    vector signed int v_ip_rot;
    vector signed int v_ic;
    vector signed int v_dc;
    vector signed int v_mc_1;
    vector signed int v_mc_2;
    vector signed int v_mc_rot;
    vector signed int v_tpdd;
    vector signed int v_tpmm;
    vector signed int v_tpim;
    vector signed int v_tpdm;
    vector signed int v_tpmi;
    vector signed int v_tpii;
    vector signed int v_tpmd;

    vector signed int v_work1;
    vector signed int v_work2;
    vector signed int v_work3;
    
    vector unsigned char perm;
    
    if(first_altivec)
      AltivecMessage();

    M = hmm->M;
            
    /* Initialization of the zero row. */
    xmx_XMN = 0;		                     /* S->N, p=1            */
    xmx_XMB = hmm->xsc[XTN][MOVE];           /* S->N->B, no N-tail   */
    xmx_XME = xmx_XMC = xmx_XMJ = -INFTY;    /* need seq to get here */
      
    workp = mx->workspace;
    
    for (n = 0; n  < M+5; n+=4) 
    {
        vec_st(m_infty,  0, workp);
        vec_st(m_infty, 16, workp);
        vec_st(m_infty, 32, workp);
        workp += 12;
    }

    tpmm  = hmm->tsc[TMM];
    tpim  = hmm->tsc[TIM];
    tpdm  = hmm->tsc[TDM];
    tpmd  = hmm->tsc[TMD];
    tpdd  = hmm->tsc[TDD];
    tpmi  = hmm->tsc[TMI];
    tpii  = hmm->tsc[TII];
    bp    = hmm->bsc;     
    ep    = hmm->esc;     
    
    /* Fill data beyound M with -INFTY, so we can take the maximum including
     * elements with k>M.
     */
    for(k=1+M;k<=(3+M);k++) {
        hmm->esc[k] = -INFTY;
        hmm->bsc[k] = -INFTY;
        
        tpmm[k] = -INFTY;
        tpim[k] = -INFTY;
        tpdm[k] = -INFTY;
        tpmd[k] = -INFTY;
        tpdd[k] = -INFTY;
        tpmi[k] = -INFTY;
        tpii[k] = -INFTY;  
        
        for(i=0;i<MAXCODE;i++) {
            hmm->msc[i][k]=-INFTY;
            hmm->isc[i][k]=-INFTY;
        }
    }

    tpmi += 1;
    tpii += 1;
    bp   += 1;
    ep   += 1;
    
    /* Recursion. Done as a pull
     * Note some slightly wasteful boundary conditions:  
     * tsc[0] = -INFTY for all eight transitions (no node 0)
     * D_M and I_M are wastefully calculated (they don't exist)
     */

    for (i = 1; i <= L; i++) 
    {
        workp = mx->workspace;

#if (defined __GNUC__) && (defined __APPLE__)
        __dcbtst(workp,256);
#endif
        
        ms    = hmm->msc[dsq[i]] + 1;
        is    = hmm->isc[dsq[i]] + 1;
        
        /* load and splat (spread to all elements) XMX[i-1][XMB] */
        v_xmb = vec_lde(0,&xmx_XMB);
        perm  = vec_lvsl(0,&xmx_XMB);
        v_xmb = vec_perm(v_xmb,v_xmb,perm);
        v_xmb = vec_splat(v_xmb,0);

        n  = 0;
        nn = 16;
        v_dc = m_infty;
        
        /* We will do the first iteration outside the loop since it is special,
         * but first we pre-prefetch the data we will need.
         */
        v_mpp_1 = vec_ld( 0, workp);
        v_ip_1  = vec_ld(16, workp);
        v_dpp_1 = vec_ld(32, workp);
        
        /* Load transition probabilities we are going to use for the match state. */   
        v_tpmm  = vec_ld( n, tpmm);
        v_tpim  = vec_ld( n, tpim);
        v_tpdm  = vec_ld( n, tpdm);
        
        v_bp    = vec_ld( n, bp);
        v_ep    = vec_ld( n, ep);
        v_ms    = vec_ld( n, ms);
        
        /* Load transition probabilities we are going to use for the insert state. */   
        v_tpmi  = vec_ld( n, tpmi);
        v_tpii  = vec_ld( n, tpii);
        v_is    = vec_ld( n, is);
        
        /* Rotate it into place with the first -infty elements,
         * so we have values from the previous column in the previous
         * row for all three of them.
         */   
        v_mpp_rot   = vec_sld(m_infty, v_mpp_1, 12);
        v_dpp_rot   = vec_sld(m_infty, v_dpp_1, 12);
        v_ip_rot    = vec_sld(m_infty, v_ip_1,  12);
        
        /* Prefetch the values we will use in the next iteration */
        v_mpp_2  = vec_ld( 48, workp);
        v_ip_2   = vec_ld( 64, workp);
        v_dpp_2  = vec_ld( 80, workp);        
        
        /* Match state */
        v_mc_1   = vec_adds(v_mpp_rot,v_tpmm);
        v_work1  = vec_adds(v_ip_rot, v_tpim);
        v_mc_1   = vec_max(v_mc_1,v_work1);
        v_work2  = vec_adds(v_dpp_rot, v_tpdm);
        v_work3  = vec_adds(v_xmb, v_bp);
        v_work2  = vec_max(v_work2,v_work3);
        v_mc_1   = vec_max(v_mc_1,v_work2);
        v_mc_1   = vec_adds(v_mc_1,v_ms);

        v_xme    = vec_adds(v_mc_1,v_ep);
        
        /* Finished match state, preload values for next iteration */
        v_tpmm  = vec_ld(nn, tpmm);
        v_tpim  = vec_ld(nn, tpim);
        v_tpdm  = vec_ld(nn, tpdm);
        
        v_bp    = vec_ld(nn, bp);
        v_ep    = vec_ld(nn, ep);
        v_ms    = vec_ld(nn, ms);
        
        /* Save the new match values to memory */
        vec_st(v_mc_1, 0, workp);
        
        v_mc_rot = vec_sld(m_infty,v_mc_1,12);
        
        /* Insertion state */
        v_ic    = vec_adds(v_mpp_1,v_tpmi);
        v_work1 = vec_adds(v_ip_1, v_tpii);
        v_ic    = vec_max(v_ic,v_work1);
        v_ic    = vec_adds(v_ic,v_is);
        
        /* Prefetch values for next iteration */
        v_tpmi  = vec_ld(nn, tpmi);
        v_tpii  = vec_ld(nn, tpii);
        v_is    = vec_ld(nn, is);        

        /* Store the new insertion state values to memory */
        vec_st(v_ic, 16, workp);
        
        /* Prefetch values for delete state */
        v_tpdd  = vec_ld(n, tpdd);
        v_tpmd  = vec_ld(n, tpmd);
        
        workp += 12;
        n     += 16;
        nn    += 16;

        /* rotate in place for next iteration */
        v_mpp_rot   = vec_sld(v_mpp_1, v_mpp_2, 12);
        v_ip_rot    = vec_sld(v_ip_1,  v_ip_2,  12);
        v_dpp_rot   = vec_sld(v_dpp_1, v_dpp_2, 12);

        for ( k=5 ; k < M-3 ; k+=8 ) 
        {

#if (defined __GNUC__) && (defined __APPLE__)
            __dcbtst(workp,256);
#endif

            /* Calculate delete state for LAST iteration */
            /* dc[k] = mc[k-1] + tpmd[k-1] */
            v_work1 = vec_adds(v_mc_rot,v_tpmd);
            
            /* Do the first propagation step of delete state for LAST iteration */
            /* calculate dc[k-1] + tpdd[k-1], and assign if larger than dc[k] */
            v_work2 = vec_sld(v_dc,     v_work1,12);
            v_work2 = vec_adds(v_work2, v_tpdd);
            v_dc    = vec_max(v_work1,  v_work2);
            
            /* Second propagation step of delete state for LAST iteration */
            v_work1 = vec_sld(m_infty,v_dc,12);
            v_work1 = vec_adds(v_work1,v_tpdd);
            v_dc    = vec_max(v_dc,v_work1);
            
            /* Preload previous row of DP matrix for next iteration.
             * Note that this version of the Viterbi routine only stores one row.
             */
            v_mpp_1  = vec_ld( 48, workp);
            v_ip_1   = vec_ld( 64, workp);
            v_dpp_1  = vec_ld( 80, workp);        
            
            /* Calculate MATCH_ALT1: mpp[k-1] + tpmm[k-1] */
            v_mc_2   = vec_adds(v_mpp_rot,v_tpmm);
            /* Calculate MATCH_ALT2: ip[k-1] + tpim[k-1] */
            v_work1  = vec_adds(v_ip_rot, v_tpim);
            /* mc[k] = max(MATCH_ALT1,MATCH_ALT2) */
            v_mc_2   = vec_max(v_mc_2,v_work1);
            /* Calculate MATCH_ALT3: dpp[k-1] + tpdm[k-1] */
            v_work2  = vec_adds(v_dpp_rot, v_tpdm);
            /* Calculate MATCH_ALT4: xmb + bp[k]  */
            v_work3  = vec_adds(v_xmb, v_bp);
            /* max of MATCH_ALT3, MATCH_ALT4 */
            v_work2  = vec_max(v_work2,v_work3);
            /* mc[k] = max(ALT1,ALT2,ALT3,ALT4) */
            v_mc_2   = vec_max(v_mc_2,v_work2);
            
            /* mc[k] += ms[k] */
            v_mc_2 = vec_adds(v_mc_2,v_ms);
            
            /* Calculation for E state that originally was done in a separate loop.
             * xme should be max of mc[k] + ep[k] in each row.
             */
            /* mc[k] + ep[k] */
            v_work1  = vec_adds(v_mc_2,v_ep);
            /* xme = max(xme, mc[k] + ep[k]) */
            v_xme    = vec_max(v_xme,v_work1);
            
            /* Finished match state, preload values for next iteration */
            v_tpmm  = vec_ld(nn, tpmm);
            v_tpim  = vec_ld(nn, tpim);
            v_tpdm  = vec_ld(nn, tpdm);
            
            /* Third propagation step of delete state for LAST iteration */
            v_work1 = vec_sld(m_infty,v_dc,12);
            v_work1 = vec_adds(v_work1,v_tpdd);
            v_dc    = vec_max(v_dc,v_work1);
            
            v_bp    = vec_ld(nn, bp);
            v_ep    = vec_ld(nn, ep);
            v_ms    = vec_ld(nn, ms);
            
            /* Save the new match values to memory */
            vec_st(v_mc_2, 0, workp);
            
            v_mc_rot = vec_sld(v_mc_1,v_mc_2,12);
            
            
            /* Prefetch for delete state */
            v_tpmd  = vec_ld( n, tpmd);

            /* Insertion state. We have taken care so all transition probabilities
             * and other vectors are filled with -2e31-1, so we don't need to treat
             * k==M as a special case.
             *
             * For the insertion state we need mpp[k] and ip[k], instead of the
             * rotated alternatives mpp[k-1] and ip[k-1]. This is already available
             * in the variables v_mpp_old and v_ip_old.   
             */  
            
            /* INSERT_ALT1: mpp[k] + tpmi[k] */
            v_ic    = vec_adds(v_mpp_2,v_tpmi);
            /* INSERT_ALT2: ip[k] + tpii[k] */
            v_work1 = vec_adds(v_ip_2, v_tpii);
            /* ic[k] = max(ALT1,ALT2) */
            v_ic    = vec_max(v_ic,v_work1);
            /* ic[k] += is[k] */
            v_ic    = vec_adds(v_ic,v_is);


            /* Prefetch values for next iteration */
            v_tpmi  = vec_ld(nn, tpmi);
            v_tpii  = vec_ld(nn, tpii);
            v_is    = vec_ld(nn, is);
            
            /* Fourth propagation step of delete state for LAST iteration */
            v_work1 = vec_sld(m_infty,v_dc,12);
            v_work1 = vec_adds(v_work1,v_tpdd);
            v_dc    = vec_max(v_dc,v_work1);
            
            /* Prefetch values for delete state */
            v_tpdd  = vec_ld( n, tpdd);
            
            
            /* Store the delete state for last iteration */
            vec_st(v_dc, -16, workp);
            
            /* rotate in place for next iteration */
            v_mpp_rot   = vec_sld(v_mpp_2, v_mpp_1, 12);
            v_ip_rot    = vec_sld(v_ip_2,  v_ip_1,  12);
            v_dpp_rot   = vec_sld(v_dpp_2, v_dpp_1, 12);
            
            /* Store the new insertion state values to memory */
            vec_st(v_ic, 16, workp);
            
            workp += 12;        
            n     += 16;
            nn    += 16;
            
            /* Second unrolled iteration */
            
            /* Calculate delete state for LAST iteration */
            /* dc[k] = mc[k-1] + tpmd[k-1] */
            v_work1 = vec_adds(v_mc_rot,v_tpmd);
            
            /* Do the first propagation step of delete state for LAST iteration */
            /* calculate dc[k-1] + tpdd[k-1], and assign if larger than dc[k] */
            v_work2 = vec_sld(v_dc,v_work1,12);
            v_work2 = vec_adds(v_work2,v_tpdd);
            v_dc    = vec_max(v_work1,v_work2);
            
            /* Second propagation step of delete state for LAST iteration */
            v_work1 = vec_sld(m_infty,v_dc,12);
            v_work1 = vec_adds(v_work1,v_tpdd);
            v_dc    = vec_max(v_dc,v_work1);
            
            /* Preload previous row of DP matrix for next iteration.
             * Note that this version of the Viterbi routine only stores one row.
             */
            v_mpp_2  = vec_ld( 48, workp);
            v_ip_2   = vec_ld( 64, workp);
            v_dpp_2  = vec_ld( 80, workp);        
            
            /* Prefetch the values we will use in the next iteration */
            
            /* Calculate MATCH_ALT1: mpp[k-1] + tpmm[k-1] */
            v_mc_1   = vec_adds(v_mpp_rot,v_tpmm);
            /* Calculate MATCH_ALT2: ip[k-1] + tpim[k-1] */
            v_work1  = vec_adds(v_ip_rot, v_tpim);
            /* mc[k] = max(MATCH_ALT1,MATCH_ALT2) */
            v_mc_1   = vec_max(v_mc_1,v_work1);
            /* Calculate MATCH_ALT3: dpp[k-1] + tpdm[k-1] */
            v_work2  = vec_adds(v_dpp_rot, v_tpdm);
            /* Calculate MATCH_ALT4: xmb + bp[k]  */
            v_work3  = vec_adds(v_xmb, v_bp);
            /* max of MATCH_ALT3, MATCH_ALT4 */
            v_work2  = vec_max(v_work2,v_work3);
            /* mc[k] = max(ALT1,ALT2,ALT3,ALT4) */
            v_mc_1   = vec_max(v_mc_1,v_work2);
            
            /* mc[k] += ms[k] */
            v_mc_1 = vec_adds(v_mc_1,v_ms);
            
            /* Calculation for E state that originally was done in a separate loop.
             * xme should be max of mc[k] + ep[k] in each row.
             */
            /* mc[k] + ep[k] */
            v_work1  = vec_adds(v_mc_1,v_ep);
            /* xme = max(xme, mc[k] + ep[k]) */
            v_xme    = vec_max(v_xme,v_work1);
            
            
            /* Finished match state, preload values for next iteration */
            v_tpmm  = vec_ld(nn, tpmm);
            v_tpim  = vec_ld(nn, tpim);
            v_tpdm  = vec_ld(nn, tpdm);
            
            /* Third propagation step of delete state for LAST iteration */
            v_work1 = vec_sld(m_infty,v_dc,12);
            v_work1 = vec_adds(v_work1,v_tpdd);
            v_dc    = vec_max(v_dc,v_work1);
            
            v_bp    = vec_ld(nn, bp);
            v_ep    = vec_ld(nn, ep);
            v_ms    = vec_ld(nn, ms);
            
            /* Save the new match values to memory */
            vec_st(v_mc_1, 0, workp);
            
            v_mc_rot = vec_sld(v_mc_2,v_mc_1,12);
            
            /* Prefetch for delete state */
            v_tpmd  = vec_ld( n, tpmd);
            
            /* Insertion state. We have taken care so all transition probabilities
             * and other vectors are filled with -2e31-1, so we don't need to treat
             * k==M as a special case.
             *
             * For the insertion state we need mpp[k] and ip[k], instead of the
             * rotated alternatives mpp[k-1] and ip[k-1]. This is already available
             * in the variables v_mpp_old and v_ip_old.   
             */  
            
            /* INSERT_ALT1: mpp[k] + tpmi[k] */
            v_ic    = vec_adds(v_mpp_1,v_tpmi);
            /* INSERT_ALT2: ip[k] + tpii[k] */
            v_work1 = vec_adds(v_ip_1, v_tpii);
            /* ic[k] = max(ALT1,ALT2) */
            v_ic    = vec_max(v_ic,v_work1);
            /* ic[k] += is[k] */
            v_ic    = vec_adds(v_ic,v_is);
            
            /* Prefetch values for next iteration */
            v_tpmi  = vec_ld(nn, tpmi);
            v_tpii  = vec_ld(nn, tpii);
            v_is    = vec_ld(nn, is);
            
            /* Fourth propagation step of delete state for LAST iteration */
            v_work1 = vec_sld(m_infty,v_dc,12);
            v_work1 = vec_adds(v_work1,v_tpdd);
            v_dc    = vec_max(v_dc,v_work1);
            
            /* Prefetch values for delete state */
            v_tpdd  = vec_ld( n, tpdd);
            
            /* Store the new insertion state values to memory */
            vec_st(v_ic, 16, workp);
            
            /* rotate in place for next iteration */
            v_mpp_rot   = vec_sld(v_mpp_1, v_mpp_2, 12);
            v_ip_rot    = vec_sld(v_ip_1,  v_ip_2,  12);
            v_dpp_rot   = vec_sld(v_dpp_1, v_dpp_2, 12);
            
            /* Store the delete state for last iteration */
            vec_st(v_dc, -16, workp);
            
            workp += 12;        
            n     += 16;
            nn    += 16;
        }
        
        /* Odd loop */
        if(k <= M)
        {
            /* Calculate delete state for LAST iteration */
            /* dc[k] = mc[k-1] + tpmd[k-1] */
            v_work1 = vec_adds(v_mc_rot,v_tpmd);
            
            /* Do the first propagation step of delete state for LAST iteration */
            /* calculate dc[k-1] + tpdd[k-1], and assign if larger than dc[k] */
            v_work2 = vec_sld(v_dc,v_work1,12);
            v_work2 = vec_adds(v_work2,v_tpdd);
            v_dc    = vec_max(v_work1,v_work2);
            
            /* Second propagation step of delete state for LAST iteration */
            v_work1 = vec_sld(m_infty,v_dc,12);
            v_work1 = vec_adds(v_work1,v_tpdd);
            v_dc    = vec_max(v_dc,v_work1);
            
            /* Prefetch the values we will use in the next iteration */
            
            /* Calculate MATCH_ALT1: mpp[k-1] + tpmm[k-1] */
            v_mc_2   = vec_adds(v_mpp_rot,v_tpmm);
            /* Calculate MATCH_ALT2: ip[k-1] + tpim[k-1] */
            v_work1  = vec_adds(v_ip_rot, v_tpim);
            /* mc[k] = max(MATCH_ALT1,MATCH_ALT2) */
            v_mc_2   = vec_max(v_mc_2,v_work1);
            /* Calculate MATCH_ALT3: dpp[k-1] + tpdm[k-1] */
            v_work2  = vec_adds(v_dpp_rot, v_tpdm);
            /* Calculate MATCH_ALT4: xmb + bp[k]  */
            v_work3  = vec_adds(v_xmb, v_bp);
            /* max of MATCH_ALT3, MATCH_ALT4 */
            v_work2  = vec_max(v_work2,v_work3);
            /* mc[k] = max(ALT1,ALT2,ALT3,ALT4) */
            v_mc_2   = vec_max(v_mc_2,v_work2);
            
            /* mc[k] += ms[k] */
            v_mc_2   = vec_adds(v_mc_2,v_ms);
            
            /* Calculation for E state that originally was done in a separate loop.
             * xme should be max of mc[k] + ep[k] in each row.
             */
            /* mc[k] + ep[k] */
            v_work1  = vec_adds(v_mc_2,v_ep);
            /* xme = max(xme, mc[k] + ep[k]) */
            v_xme    = vec_max(v_xme,v_work1);
            
            /* Third propagation step of delete state for LAST iteration */
            v_work1 = vec_sld(m_infty,v_dc,12);
            v_work1 = vec_adds(v_work1,v_tpdd);
            v_dc    = vec_max(v_dc,v_work1);
            
            /* Save the new match values to memory */
            vec_st(v_mc_2, 0, workp);
            
            
            v_mc_rot = vec_sld(v_mc_1,v_mc_2,12);
            
            
            /* Insertion state. We have taken care so all transition probabilities
             * and other vectors are filled with -2e31-1, so we don't need to treat
             * k==M as a special case.
             *
             * For the insertion state we need mpp[k] and ip[k], instead of the
             * rotated alternatives mpp[k-1] and ip[k-1]. This is already available
             * in the variables v_mpp_old and v_ip_old.   
             */  
            
            /* INSERT_ALT1: mpp[k] + tpmi[k] */
            v_ic    = vec_adds(v_mpp_2,v_tpmi);
            /* INSERT_ALT2: ip[k] + tpii[k] */
            v_work1 = vec_adds(v_ip_2, v_tpii);
            /* ic[k] = max(ALT1,ALT2) */
            v_ic    = vec_max(v_ic,v_work1);
            /* ic[k] += is[k] */
            v_ic    = vec_adds(v_ic,v_is);
            
            
            /* Fourth propagation step of delete state for LAST iteration */
            v_work1 = vec_sld(m_infty,v_dc,12);
            v_work1 = vec_adds(v_work1,v_tpdd);
            v_dc    = vec_max(v_dc,v_work1);

            /* Prefetch values for delete state */
            v_tpdd  = vec_ld( n, tpdd);
            v_tpmd  = vec_ld( n, tpmd);

            /* Store the new insertion state values to memory */
            vec_st(v_ic, 16, workp);
            
            /* Store the delete state for last iteration */
            vec_st(v_dc, -16, workp);
            
            workp += 12;   
        }
        /* Do the epilogue */
            
        /* Calculate delete state for LAST iteration */
        /* dc[k] = mc[k-1] + tpmd[k-1] */
        v_work1 = vec_adds(v_mc_rot,v_tpmd);
        
        /* Do the first propagation step of delete state for LAST iteration */
        /* calculate dc[k-1] + tpdd[k-1], and assign if larger than dc[k] */
        v_work2 = vec_sld(v_dc,v_work1,12);
        v_work2 = vec_adds(v_work2,v_tpdd);
        v_dc    = vec_max(v_work1,v_work2);
        
        /* Second propagation step of delete state for LAST iteration */
        v_work1 = vec_sld(m_infty,v_dc,12);
        v_work1 = vec_adds(v_work1,v_tpdd);
        v_dc    = vec_max(v_dc,v_work1);
        
        /* Third propagation step of delete state for LAST iteration */
        v_work1 = vec_sld(m_infty,v_dc,12);
        v_work1 = vec_adds(v_work1,v_tpdd);
        v_dc    = vec_max(v_dc,v_work1);
        
        /* Fourth propagation step of delete state for LAST iteration */
        v_work1 = vec_sld(m_infty,v_dc,12);
        v_work1 = vec_adds(v_work1,v_tpdd);
        v_dc    = vec_max(v_dc,v_work1);
        
        /* Store the delete state for last iteration */
        vec_st(v_dc, -16, workp);
        
        
        
        /* end of k loop */

        
        /* Now the special states. Order is important here.
         * remember, C and J emissions are zero score by definition,
         */
        /* N state */
        if ((sc = xmx_XMN + hmm->xsc[XTN][LOOP]) > -INFTY)
            xmx_XMN = sc;
        else
            xmx_XMN = -INFTY;
        
        
        /* E state */

        v_work1 = vec_sld(v_xme,v_xme,8);
        v_work1 = vec_max(v_work1,v_xme);
        v_work2 = vec_sld(v_work1,v_work1,4);
        v_work1 = vec_max(v_work1,v_work2);
        
        vec_ste(v_work1,0,&xmx_XME);

        if(xmx_XME < -INFTY)
            xmx_XME = -INFTY;
        
        /* J state */
        if ((sc = xmx_XMJ + hmm->xsc[XTJ][LOOP]) > -INFTY)
            xmx_XMJ = sc;
        else
            xmx_XMJ = -INFTY;
        
        if ((sc = xmx_XME   + hmm->xsc[XTE][LOOP]) > xmx_XMJ)
            xmx_XMJ = sc;
        
        /* B state */
        if ((sc = xmx_XMN + hmm->xsc[XTN][MOVE]) > -INFTY)
            xmx_XMB = sc;
        else
            xmx_XMB = -INFTY;
        
        if ((sc = xmx_XMJ + hmm->xsc[XTJ][MOVE]) > xmx_XMB)
            xmx_XMB = sc;
        
        /* C state */
        if ((sc = xmx_XMC + hmm->xsc[XTC][LOOP]) > -INFTY)
            xmx_XMC = sc;
        else
            xmx_XMC = -INFTY;
        
        if ((sc = xmx_XME + hmm->xsc[XTE][MOVE]) > xmx_XMC)
            xmx_XMC = sc;
        
    }
    /* T state (not stored) */
    sc = xmx_XMC + hmm->xsc[XTC][MOVE];
        
    return Scorify(sc);
}
