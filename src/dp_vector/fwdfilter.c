/* Forwards/Backwards filters.
 *
 * This is code in the acceleration pipeline:
 * SSVFilter -> MSVFilter -> VitFilter -> ForwardFilter -> BackwardFilter
 *                                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 *                                                (you are here)
 *
 * Striped SIMD vector implementation [Farrar07], using checkpointing 
 * to guarantee O(M sqrt L) memory [Grice97,TarnasHughey98,Newberg08],
 * in probability space using sparse rescaling [Eddy11].
 *
 * Called in two pieces. p7_ForwardFilter() returns the Forward score
 * in nats. Caller chooses whether or not to proceed with Backward and
 * posterior decoding. If caller proceeds, p7_BackwardFilter() does
 * Backward and posterior decoding. Based on the posterior decoding
 * probabilities on each row i, it determines which cells are to be
 * added to the sparse DP mask for subsequent local/glocal
 * reprocessing.  The sparse DP mask is returned in a <P7_SPARSEMASK>
 * structure.
 * 
 * Any cell (i,k) with total posterior probability (i.e., summed over
 * M,D,I) >= sm_thresh is marked and included in the sparse mask.
 * Cells needed for glocal entry/exit delete paths (G->DDDD->Mk,
 * Mk->DDDD->E) are not marked, because the filter only uses local
 * alignment, not glocal. The default for sm_thresh is
 * <p7_SPARSIFY_THRESH>, in <p7_config.h.in>; currently set
 * to 0.01.
 *
 * ForwardFilter() and BackwardFilter() are a dependent pair, sharing
 * the same DP matrix object. They must be called sequentially,
 * because BackwardFilter() is doing posterior decoding on the fly,
 * and for this it needs both forward and backward scores; it uses
 * checkpointed Forward rows that have been stored by the preceding
 * ForwardFilter() call.
 * 
 * Prob-space implementations using sparse rescaling require multihit
 * local alignment mode for numeric range reasons. Unihit or glocal
 * will result in errors due to underflow.
 *
 * Contents:
 *    1. ForwardFilter() and BackwardFilter() API.
 *    2. Internal functions: inlined recursions.
 *    3. Debugging tools.
 *    4. Stats driver (memory requirements)
 *    5. Benchmark driver.
 *    6. Unit tests.
 *    7. Test driver.
 *    8. Example.
 *    9. Notes
 *       a. On debugging and testing methods.
 *       b. On running time, in theory and in practice.
 *       c. Copyright and license information.
 */
#include "p7_config.h"

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"

#ifdef p7_build_AVX2
#include <immintrin.h>
#include "esl_avx.h"
#endif
#ifdef p7_build_AVX512
#include <immintrin.h>
#include "esl_avx_512.h"
#endif
#include "esl_vectorops.h"

#include "dp_reference/p7_refmx.h"
#include "dp_sparse/p7_sparsemx.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/fwdfilter.h"


/* vectorized DP recursions are tediously lengthy, so for some
 * semblance of clarity, they're broken out into one-page-ish
 * chunks, using static inlined functions.
 */
#ifdef p7_build_SSE
static inline float forward_row      (ESL_DSQ xi, const P7_OPROFILE *om, const __m128 *dpp, __m128 *dpc, int Q);
static inline void  backward_row_main(ESL_DSQ xi, const P7_OPROFILE *om,       __m128 *dpp, __m128 *dpc, int Q, float scalefactor);
static inline void  backward_row_L   (            const P7_OPROFILE *om,                    __m128 *dpc, int Q, float scalefactor);
static inline void  backward_row_finish(          const P7_OPROFILE *om,                    __m128 *dpc, int Q, __m128 dcv);
static inline void  backward_row_rescale(float *xc, __m128 *dpc, int Q, float scalefactor);
static inline int   posterior_decode_row(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc);
#endif

#ifdef p7_build_AVX2
static inline float
forward_row_AVX(ESL_DSQ xi, const P7_OPROFILE *om, const __m256 *dpp, __m256 *dpc, int Q);
static inline void  backward_row_main_AVX(ESL_DSQ xi, const P7_OPROFILE *om,       __m256 *dpp, __m256 *dpc, int Q, float scalefactor);
static inline void  backward_row_L_AVX   (            const P7_OPROFILE *om,                    __m256 *dpc, int Q, float scalefactor);
static inline void  backward_row_finish_AVX(          const P7_OPROFILE *om,                    __m256 *dpc, int Q, __m256 dcv);
static inline void  backward_row_rescale_AVX(float *xc, __m256 *dpc, int Q, float scalefactor);
static inline int   posterior_decode_row_AVX(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc);
#endif

#ifdef p7_build_AVX512
static inline float
forward_row_AVX_512(ESL_DSQ xi, const P7_OPROFILE *om, const __m512 *dpp, __m512 *dpc, int Q);
static inline void  backward_row_main_AVX_512(ESL_DSQ xi, const P7_OPROFILE *om,       __m512 *dpp, __m512 *dpc, int Q, float scalefactor);
static inline void  backward_row_L_AVX_512   (            const P7_OPROFILE *om,                    __m512 *dpc, int Q, float scalefactor);
static inline void  backward_row_finish_AVX_512(          const P7_OPROFILE *om,                    __m512 *dpc, int Q, __m512 dcv);
static inline void  backward_row_rescale_AVX_512(float *xc, __m512 *dpc, int Q, float scalefactor);
static inline int   posterior_decode_row_AVX_512(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc);
#endif

#ifdef p7_DEBUGGING
static inline float backward_row_zero(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox);
static        void  save_debug_row_pp(P7_CHECKPTMX *ox,               __m128 *dpc, int i);
static        void  save_debug_row_fb(P7_CHECKPTMX *ox, P7_REFMX *gx, __m128 *dpc, int i, float totscale);

#endif

/*****************************************************************
 * 1. Forward and Backward API calls
 *****************************************************************/

/* Function:  p7_ForwardFilter()
 * Synopsis:  Checkpointed striped vector Forward calculation, producing score.
 *
 * Purpose:   Calculate the Forward algorithm for target sequence <dsq>
 *            of <L> residues aligned to query profile <om>, using the
 *            checkpointed DP matrix <ox> provided by the caller. Upon
 *            successful return, <ox> contains the filled Forward
 *            matrix, and <*opt_sc> optionally contains the raw Forward
 *            score in nats.
 *            
 *            <ox> will be reallocated, if needed, for the MxL problem;
 *            caller does not need to call <p7_checkptmx_GrowTo()> itself.
 *  
 * Args:      dsq    - digital target sequence, 1..L
 *            L      - length of dsq, residues
 *            om     - optimized profile (multihit local)
 *            ox     - checkpointed DP matrix
 *            opt_sc - optRETURN: raw Forward score (nats)
 *
 * Returns:   <eslOK> on success, <ox> contains the checkpointed
 *            Forward matrix calculation, and <*opt_sc> optionally
 *            has the raw Forward score in nats.
 *
 * Throws:    <eslEMEM> on reallocation error.
 * 
 * Xref:      For layout of checkpointed <ox> see exegesis in p7_checkptmx.h.
 */
int
p7_ForwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc)
{
#ifdef p7_build_SSE
  int           Q     = P7_NVF(om->M);                   /* segment length; # of MDI vectors on each row      */
  __m128       *dpp   = NULL;                            /* dpp=prev row. start on dpf[2]; rows 0,1=Backwards */
  __m128       *dpc   = NULL;		                 /* dpc points at current row         */
  const __m128  zerov = _mm_setzero_ps();     
  float        *xc    = NULL;                            /* specials E,N,JJ,J,B,CC,C,SCALE    */
  float         totsc = 0.0f;	                         /* accumulates Forward score in nats */
  int     q;      /* counter over vectors 0..Q-1                        */
   ox->Qf = Q;
#endif
#ifdef p7_build_AVX2
  int           Q_AVX     = P7_NVF_AVX(om->M);                   /* segment length; # of MDI vectors on each row      */
  __m256       *dpp_AVX   = NULL;                            /* dpp=prev row. start on dpf[2]; rows 0,1=Backwards */
  __m256       *dpc_AVX   = NULL;                    /* dpc points at current row         */
  const __m256  zerov_AVX = _mm256_setzero_ps();     
  float        *xc_AVX    = NULL;                            /* specials E,N,JJ,J,B,CC,C,SCALE    */
  float         totsc_AVX = 0.0f;                          /* accumulates Forward score in nats */
  int     q_AVX;      /* counter over vectors 0..Q-1                        */
   ox->Qf_AVX = Q_AVX;
#endif
#ifdef p7_build_AVX512
  int           Q_AVX_512     = P7_NVF_AVX_512(om->M);                   /* segment length; # of MDI vectors on each row      */
  __m512       *dpp_AVX_512   = NULL;                            /* dpp=prev row. start on dpf[2]; rows 0,1=Backwards */
  __m512       *dpc_AVX_512   = NULL;                    /* dpc points at current row         */
  const __m512  zerov_AVX_512 = _mm512_setzero_ps();     
  float        *xc_AVX_512    = NULL;                            /* specials E,N,JJ,J,B,CC,C,SCALE    */
  float         totsc_AVX_512 = 0.0f;                          /* accumulates Forward score in nats */
  int     q_AVX_512;      /* counter over vectors 0..Q-1                        */
   ox->Qf_AVX_512 = Q_AVX_512;
#endif


  int     i;			/* counter over residues/rows 1..L                    */
  int     b;			/* counter down through checkpointed blocks, Rb+Rc..1 */
  int     w;			/* counter down through rows in a checkpointed block  */

  /* Contract checks */
  //  ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  //  ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  //  ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
                                          /*  ... which you can disable, if you're playing w/ config     */
  // SRE: disabled 21 Jul 14 because unit tests in sparse_asc_fwdback need to create
  // sparse mask w/ checkpointed F/B w/ unilocal L=0 models

  /* Make sure <ox> is allocated big enough.
   * DO NOT set any ptrs into the matrix until after this potential reallocation!
   */
  p7_checkptmx_GrowTo(ox, om->M, L);
 #ifdef p7_build_SSE
  dpp =  (__m128 *) ox->dpf[ox->R0-1];    /* dpp=prev row. start on dpf[2]; rows 0,1=Backwards */
  xc  =  (float *) (dpp + Q*p7C_NSCELLS); /* specials E,N,JJ,J,B,CC,C,SCALE    */
#endif
#ifdef p7_build_AVX2
  dpp_AVX =  (__m256 *) ox->dpf_AVX[ox->R0-1];    /* dpp=prev row. start on dpf[2]; rows 0,1=Backwards */
  xc_AVX  =  (float *) (dpp_AVX + Q_AVX*p7C_NSCELLS); /* specials E,N,JJ,J,B,CC,C,SCALE    */
#endif
#ifdef p7_build_AVX512
  dpp_AVX_512 =  (__m512 *) ox->dpf_AVX_512[ox->R0-1];    /* dpp=prev row. start on dpf[2]; rows 0,1=Backwards */
  xc_AVX_512  =  (float *) (dpp_AVX_512 + Q_AVX_512*p7C_NSCELLS); /* specials E,N,JJ,J,B,CC,C,SCALE    */
#endif
  /* Set the size of the problem in <ox> now, not later
   * Debugging dumps need this information, for example
   * (Ditto for any debugging copy of the fwd mx)
   */
  ox->M  = om->M;	
  ox->L  = L;
 
#ifdef p7_DEBUGGING
  ox->dump_flags |= p7_SHOW_LOG;                     /* also sets for Backward dumps, since <ox> shared */
  if (ox->do_dumping) p7_checkptmx_DumpFBHeader(ox);
  if (ox->fwd) { ox->fwd->M = om->M; ox->fwd->L = L; ox->fwd->type = p7R_FORWARD; }
#endif

  /* Initialization of the zero row, including specials */
#ifdef p7_build_SSE
  for (q = 0; q < p7C_NSCELLS*Q; q++) dpp[q] = zerov;
  xc[p7C_N]     = 1.;
  xc[p7C_B]     = om->xf[p7O_N][p7O_MOVE]; 
  xc[p7C_E]     = xc[p7C_JJ] = xc[p7C_J]  = xc[p7C_CC] = xc[p7C_C]  = 0.;			
  xc[p7C_SCALE] = 1.;
#endif
#ifdef p7_build_AVX2
  for (q_AVX = 0; q_AVX < p7C_NSCELLS*Q_AVX; q_AVX++) dpp_AVX[q_AVX] = zerov_AVX;
  xc_AVX[p7C_N]     = 1.;
  xc_AVX[p7C_B]     = om->xf[p7O_N][p7O_MOVE]; 
  xc_AVX[p7C_E]     = xc_AVX[p7C_JJ] = xc_AVX[p7C_J]  = xc_AVX[p7C_CC] = xc_AVX[p7C_C]  = 0.;     
  xc_AVX[p7C_SCALE] = 1.;
#endif           
#ifdef p7_build_AVX512
  for (q_AVX_512 = 0; q_AVX_512 < p7C_NSCELLS*Q_AVX_512; q_AVX_512++) dpp_AVX_512[q_AVX_512] = zerov_AVX_512;
  xc_AVX_512[p7C_N]     = 1.;
  xc_AVX_512[p7C_B]     = om->xf[p7O_N][p7O_MOVE]; 
  xc_AVX_512[p7C_E]     = xc_AVX_512[p7C_JJ] = xc_AVX_512[p7C_J]  = xc_AVX_512[p7C_CC] = xc_AVX_512[p7C_C]  = 0.;     
  xc_AVX_512[p7C_SCALE] = 1.;
#endif  

#ifdef p7_DEBUGGING
  if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, 0, dpp, "f1 O"); 
  if (ox->fwd)        save_debug_row_fb(ox, ox->fwd, dpp, 0, totsc); 
#endif

  /* Phase one: the "a" region: all rows in this region are saved */
  for (i = 1; i <= ox->La; i++)
    {
#ifdef p7_build_SSE      
      dpc = (__m128 *) ox->dpf[ox->R0+ox->R]; ox->R++;    /* idiomatic for "get next save/checkpoint row" */

      totsc += forward_row(dsq[i], om, dpp, dpc, Q);
      dpp = dpc;	 /* current row becomes prev row */
#endif
#ifdef p7_build_AVX2      
      dpc_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX]; ox->R_AVX++;    /* idiomatic for "get next save/checkpoint row" */

      totsc_AVX += forward_row_AVX(dsq[i], om, dpp_AVX, dpc_AVX, Q_AVX);
      dpp_AVX = dpc_AVX;   /* current row becomes prev row */
#endif 

#ifdef p7_build_AVX512      
      dpc_AVX_512 = (__m512 *) ox->dpf_AVX_512[ox->R0+ox->R_AVX_512]; ox->R_AVX_512++;    /* idiomatic for "get next save/checkpoint row" */

      totsc_AVX_512 += forward_row_AVX_512(dsq[i], om, dpp_AVX_512, dpc_AVX_512, Q_AVX_512);
      dpp_AVX_512 = dpc_AVX_512;   /* current row becomes prev row */
#endif       
	                          
#ifdef p7_build_check_AVX2

    float *unstriped, *unstriped_AVX;
    union { __m128 v; float x[4]; } tmp_check;
    union { __m256 v; float x[8]; } tmp_check_AVX;

    unstriped = malloc( sizeof(float) * ((Q*4)+1));  // Yes, these allocates are slow,but this is check code that won't be
    unstriped_AVX = malloc( sizeof(float) * ((Q_AVX*8)+1));  // compiled in production
    int q_temp;
    int z_temp;
    unstriped[0] = 0.;
    unstriped_AVX[0] = 0.;

    /* Line 1. M cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_MQ(dpc, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
  for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_MQ(dpc_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }

  // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("forward filter M miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", i, q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }

  /* Line 2: I cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_IQ(dpc, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_IQ(dpc_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("forward filter I miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", i, q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }

 /* Line 3: D cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_DQ(dpc, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_DQ(dpc_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("forward filter D miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", i, q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }
  free(unstriped);
  free(unstriped_AVX); 

#endif      

#ifdef p7_build_check_AVX512

    float *unstriped, *unstriped_AVX_512;
    union { __m128 v; float x[4]; } tmp_check;
    union { __m512 v; float x[16]; } tmp_check_AVX_512;

    unstriped = malloc( sizeof(float) * ((Q*4)+1));  // Yes, these allocates are slow,but this is check code that won't be
    unstriped_AVX_512 = malloc( sizeof(float) * ((Q_AVX_512*16)+1));  // compiled in production
    int q_temp;
    int z_temp;
    unstriped[0] = 0.;
    unstriped_AVX_512[0] = 0.;

    /* Line 1. M cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_MQ(dpc, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
  for (q_temp = 0; q_temp < Q_AVX_512; q_temp++) {
    tmp_check_AVX_512.v = P7C_MQ(dpc_AVX_512, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX_512; z_temp++) unstriped_AVX_512[q_temp+Q_AVX_512*z_temp+1] = tmp_check_AVX_512.x[z_temp];
  }

  // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX_512[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("forward filter M miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX-512)\n", i, q_temp, unstriped[q_temp], unstriped_AVX_512[q_temp]);
    }
  }

  /* Line 2: I cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_IQ(dpc, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX_512; q_temp++) {
    tmp_check_AVX_512.v = P7C_IQ(dpc_AVX_512, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX_512; z_temp++) unstriped_AVX_512[q_temp+Q_AVX_512*z_temp+1] = tmp_check_AVX_512.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX_512[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("forward filter I miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX-512)\n", i, q_temp, unstriped[q_temp], unstriped_AVX_512[q_temp]);
    }
  }

 /* Line 3: D cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_DQ(dpc, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX_512; q_temp++) {
    tmp_check_AVX_512.v = P7C_DQ(dpc_AVX_512, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX_512; z_temp++) unstriped_AVX_512[q_temp+Q_AVX_512*z_temp+1] = tmp_check_AVX_512.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX_512[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("forward filter D miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX-512)\n", i, q_temp, unstriped[q_temp], unstriped_AVX_512[q_temp]);
    }
  }
  free(unstriped);
  free(unstriped_AVX_512); 

#endif      


#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, dpc, "f1 O"); 
      if (ox->fwd)        save_debug_row_fb(ox, ox->fwd, dpc, i, totsc); 
#endif
    }

  /* Phase two: "b" and "c" regions: partially and fully checkpointed */
  for (b = ox->Rb + ox->Rc, w = (ox->Rb ? ox->Lb : ox->Rc+1); i <= L; i++)
    {
      /* this section, deciding whether to get a checkpointed vs. tmp_check
       * row, is why we set <fwd>, <dpp> here, rather than having the
       * inlined forward_row() set them.
       */
      if (! (--w)) { 		                   /* we're on the last row in segment: this row is saved    */
#ifdef p7_build_SSE       
	dpc = (__m128 *) ox->dpf[ox->R0+ox->R]; ox->R++;  /* idiomatic for "get next save/checkpoint row"    */
#endif    
#ifdef p7_build_AVX2      
  dpc_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX]; ox->R_AVX++;  /* idiomatic for "get next save/checkpoint row"    */
#endif  
#ifdef p7_build_AVX512      
  dpc_AVX_512 = (__m512 *) ox->dpf_AVX_512[ox->R0+ox->R_AVX_512]; ox->R_AVX_512++;  /* idiomatic for "get next save/checkpoint row"    */
#endif    

	w = b;  			           /* next segment has this many rows, ending in a saved row */
	b--;					   /* decrement segment number counter; last segment is r=1  */
      } else{
#ifdef p7_build_SSE        
       dpc = (__m128 *) ox->dpf[i%2];        /* idiomatic for "next tmp row", 0/1; i%2 makes sure dpp != dpc */
#endif
#ifdef p7_build_AVX2        
       dpc_AVX = (__m256 *) ox->dpf_AVX[i%2];        /* idiomatic for "next tmp row", 0/1; i%2 makes sure dpp != dpc */
#endif
#ifdef p7_build_AVX512        
       dpc_AVX_512 = (__m512 *) ox->dpf_AVX_512[i%2];        /* idiomatic for "next tmp row", 0/1; i%2 makes sure dpp != dpc */
#endif      
      }

#ifdef p7_build_SSE      
      totsc += forward_row(dsq[i], om, dpp, dpc, Q);
      dpp = dpc;
#endif
#ifdef p7_build_AVX2      
      totsc_AVX += forward_row_AVX(dsq[i], om, dpp_AVX, dpc_AVX, Q_AVX);
      dpp_AVX = dpc_AVX;
#endif 
#ifdef p7_build_AVX512      
      totsc_AVX_512 += forward_row_AVX_512(dsq[i], om, dpp_AVX_512, dpc_AVX_512, Q_AVX_512);
      dpp_AVX_512 = dpc_AVX_512;
#endif 
#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, dpc, w ? "f1 X" : "f1 O"); 
      if (ox->fwd)        save_debug_row_fb(ox, ox->fwd, dpc, i, totsc); 
#endif
#ifdef p7_build_check_AVX2

    float *unstriped, *unstriped_AVX;
    union { __m128 v; float x[4]; } tmp_check;
    union { __m256 v; float x[8]; } tmp_check_AVX;

    unstriped = malloc( sizeof(float) * ((Q*4)+1));  // Yes, these allocates are slow,but this is check code that won't be
    unstriped_AVX = malloc( sizeof(float) * ((Q_AVX*8)+1));  // compiled in production
    int q_temp;
    int z_temp;
    unstriped[0] = 0.;
    unstriped_AVX[0] = 0.;

    /* Line 1. M cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_MQ(dpc, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
  for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_MQ(dpc_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }

  // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("forward filter M miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", i, q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }

  /* Line 2: I cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_IQ(dpc, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_IQ(dpc_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("forward filter I miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", i, q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }

 /* Line 3: D cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_DQ(dpc, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_DQ(dpc_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("forward filter D miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", i, q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }
  free(unstriped);
  free(unstriped_AVX); 

#endif      
  #ifdef p7_build_check_AVX512

    float *unstriped, *unstriped_AVX_512;
    union { __m128 v; float x[4]; } tmp_check;
    union { __m512 v; float x[16]; } tmp_check_AVX_512;

    unstriped = malloc( sizeof(float) * ((Q*4)+1));  // Yes, these allocates are slow,but this is check code that won't be
    unstriped_AVX_512 = malloc( sizeof(float) * ((Q_AVX_512*16)+1));  // compiled in production
    int q_temp;
    int z_temp;
    unstriped[0] = 0.;
    unstriped_AVX_512[0] = 0.;

    /* Line 1. M cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_MQ(dpc, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
  for (q_temp = 0; q_temp < Q_AVX_512; q_temp++) {
    tmp_check_AVX_512.v = P7C_MQ(dpc_AVX_512, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX_512; z_temp++) unstriped_AVX_512[q_temp+Q_AVX_512*z_temp+1] = tmp_check_AVX_512.x[z_temp];
  }

  // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX_512[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("forward filter M miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX-512)\n", i, q_temp, unstriped[q_temp], unstriped_AVX_512[q_temp]);
    }
  }

  /* Line 2: I cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_IQ(dpc, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX_512; q_temp++) {
    tmp_check_AVX_512.v = P7C_IQ(dpc_AVX_512, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX_512; z_temp++) unstriped_AVX_512[q_temp+Q_AVX_512*z_temp+1] = tmp_check_AVX_512.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX_512[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("forward filter I miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX-512)\n", i, q_temp, unstriped[q_temp], unstriped_AVX_512[q_temp]);
    }
  }

 /* Line 3: D cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_DQ(dpc, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX_512; q_temp++) {
    tmp_check_AVX_512.v = P7C_DQ(dpc_AVX_512, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX_512; z_temp++) unstriped_AVX_512[q_temp+Q_AVX_512*z_temp+1] = tmp_check_AVX_512.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX_512[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("forward filter D miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX-512)\n", i, q_temp, unstriped[q_temp], unstriped_AVX_512[q_temp]);
    }
  }
  free(unstriped);
  free(unstriped_AVX_512); 

#endif      
    }

#ifdef p7_build_SSE
  xc     = (float *) (dpc + Q*p7C_NSCELLS);

  ESL_DASSERT1( (ox->R == ox->Ra+ox->Rb+ox->Rc) );
  ESL_DASSERT1( ( (! isnan(xc[p7C_C])) && (! isinf(xc[p7C_C]))) );

  if (opt_sc) *opt_sc = totsc + logf(xc[p7C_C] * om->xf[p7O_C][p7O_MOVE]);
#endif
#ifdef p7_build_AVX2
  xc_AVX     = (float *) (dpc_AVX + Q_AVX*p7C_NSCELLS);

  ESL_DASSERT1( (ox->R_AVX == ox->Ra+ox->Rb+ox->Rc) );
  ESL_DASSERT1( ( (! isnan(xc_AVX[p7C_C])) && (! isinf(xc_AVX[p7C_C]))) );

  if (opt_sc) *opt_sc = totsc_AVX + logf(xc_AVX[p7C_C] * om->xf[p7O_C][p7O_MOVE]);
#endif
#ifdef p7_build_AVX512
  xc_AVX_512     = (float *) (dpc_AVX_512 + Q_AVX_512*p7C_NSCELLS);

  ESL_DASSERT1( (ox->R_AVX_512 == ox->Ra+ox->Rb+ox->Rc) );
  ESL_DASSERT1( ( (! isnan(xc_AVX_512[p7C_C])) && (! isinf(xc_AVX_512[p7C_C]))) );

  if (opt_sc) *opt_sc = totsc_AVX_512 + logf(xc_AVX_512[p7C_C] * om->xf[p7O_C][p7O_MOVE]);
#endif

#ifdef p7_build_check_AVX2
  float check = totsc + logf(xc[p7C_C] * om->xf[p7O_C][p7O_MOVE]);
  float check_AVX = totsc_AVX + logf(xc_AVX[p7C_C] * om->xf[p7O_C][p7O_MOVE]);
  if(fabs(check_AVX - check) > fabs(check/1000)){
    printf("SSE and AVX results disagree in forward filter %f vs %f\n", check, check_AVX);
  }
/*  else{
    printf("SSE and AVX results agree in forward filter %f and %f", totsc, totsc_AVX);
  } */
#endif  
#ifdef p7_build_check_AVX512
  float check = totsc + logf(xc[p7C_C] * om->xf[p7O_C][p7O_MOVE]);
  float check_AVX_512 = totsc_AVX_512 + logf(xc_AVX_512[p7C_C] * om->xf[p7O_C][p7O_MOVE]);
  if(fabs(check_AVX_512 - check) > fabs(check/1000)){
    printf("SSE and AVX-512 results disagree in forward filter %f vs %f\n", check, check_AVX_512);
  }
/*  else{
    printf("SSE and AVX results agree in forward filter %f and %f", totsc, totsc_AVX);
  } */
#endif
  return eslOK;
}


/* Function:  p7_BackwardFilter()
 * Synopsis:  Checkpointed striped vector Backward calculation, producing sparse mask.
 *
 * Purpose:   Given a target sequence <dsq> of length <L>, a query model
 *            <om>, and a DP matrix <ox> resulting from a successful
 *            call to <p7_ForwardFilter()>; calculate the Backward and
 *            posterior decoding algorithms. On each row <i=1..L>, use
 *            posterior decoding to determine which <i,k> cells pass a
 *            significance threshold for inclusion in the sparse DP
 *            mask. Store those sparse cells in <sm>, which the caller
 *            allocates (or reuses) and provides. <sm> will be resized
 *            here as needed; caller does not need to call 
 *            <p7_sparsemask_Reinit()> on it.
 *            
 *            
 * Args:      dsq       - digital target sequence, 1..L
 *            L         - length of dsq, residues
 *            om        - optimized profile (multihit local)
 *            ox        - checkpointed DP matrix, ForwardFilter already run
 *            sm        - allocated P7_SPARSEMASK structure to hold sparse DP mask
 *            sm_thresh - Threshold for determining 'significant' posterior
 *                        alignment probability, passed in to posterior_decode_row()
 *
 * Throws:    <eslEINVAL> if something's awry with a data structure's internals.
 *            <eslEMEM> on allocation failure.
 */
int
p7_BackwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh)
{
//printf("Calling Backward Filter\n");
#ifdef p7_build_SSE
    float  *xf;
  int Q = ox->Qf;
  __m128 *fwd;
  __m128 *bck;
  __m128 *dpp;
  int i;
#endif

#ifdef p7_build_AVX2
    float  *xf_AVX;
  int Q_AVX = ox->Qf_AVX;
  __m256 *fwd_AVX;
  __m256 *bck_AVX;
  __m256 *dpp_AVX;
  int i_AVX;
#endif

#ifdef p7_build_AVX512
    float  *xf_AVX_512;
  int Q_AVX_512 = ox->Qf_AVX_512;
  __m512 *fwd_AVX_512;
  __m512 *bck_AVX_512;
  __m512 *dpp_AVX_512;
  int i_AVX_512;
#endif


  float   Tvalue, Tvalue_AVX, Tvalue_AVX_512;
  int     b, w, i2;
  int     status;

  p7_sparsemask_Reinit(sm, om->M, L);

  /* Contract checks */
  //ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  //ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  //ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
                                          /*  ... which you can disable, if you're playing w/ config     */
  // SRE: disabled 21 Jul 14 because unit tests in sparse_asc_fwdback need to create
  // sparse mask w/ checkpointed F/B w/ unilocal L=0 models

#ifdef p7_DEBUGGING
  /* Debugging instrumentations. */
   if (ox->bck) { ox->bck->M = om->M; ox->bck->L = L; ox->bck->type = p7R_BACKWARD; }
   if (ox->pp)  { ox->pp->M  = om->M; ox->pp->L  = L; ox->pp->type  = p7R_DECODING; }
#endif

  /* Row L is a special case for Backwards; so in checkpointed Backward,
   * we have to special-case the last block, rows L-1 and L
   */  
 
 

#ifdef p7_build_SSE
  i = L;
  ox->R--;
  fwd = (__m128 *) ox->dpf[ox->R0 + ox->R];      /* pop row for fwd[L] off the checkpointed stack */
  xf  = (float *) (fwd + Q*p7C_NSCELLS);
  Tvalue = xf[p7C_C] * om->xf[p7O_C][p7O_MOVE];  /* i.e. scaled fwd[L] val at T state = scaled overall score */
  bck = (__m128 *) ox->dpf[i%2];	         /* get tmp space for bck[L]                                 */
  backward_row_L(om, bck, Q, xf[p7C_SCALE]);     /* calculate bck[L] row                                     */

#ifdef p7_DEBUGGING
  ox->bcksc = logf(xf[p7C_SCALE]);
  if (ox->do_dumping) { p7_checkptmx_DumpFBRow(ox, L, fwd, "f2 O"); if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, L, bck, "bck");  }
  if (ox->bck)          save_debug_row_fb(ox, ox->bck, bck, L, ox->bcksc); 

#endif
  if ( (status = posterior_decode_row(ox, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
  i--;
  dpp = bck;
#endif

#ifdef p7_build_AVX2
   i_AVX = L;
  ox->R_AVX--;
  fwd_AVX = (__m256 *) ox->dpf_AVX[ox->R0 + ox->R_AVX];      /* pop row for fwd[L] off the checkpointed stack */
  xf_AVX  = (float *) (fwd_AVX + Q_AVX*p7C_NSCELLS);
  Tvalue_AVX = xf_AVX[p7C_C] * om->xf[p7O_C][p7O_MOVE];  /* i.e. scaled fwd[L] val at T state = scaled overall score */
  bck_AVX = (__m256 *) ox->dpf_AVX[i_AVX%2];           /* get tmp space for bck[L]                                 */
  backward_row_L_AVX(om, bck_AVX, Q_AVX, xf_AVX[p7C_SCALE]);     /* calculate bck[L] row                                     */

  if ( (status = posterior_decode_row_AVX(ox, i_AVX, sm, sm_thresh, Tvalue_AVX)) != eslOK) return status;
  i_AVX--;
  dpp_AVX = bck_AVX;
#endif

#ifdef p7_build_check_AVX2  // Check that the computed rows match

    float *unstriped, *unstriped_AVX;
    union { __m128 v; float x[4]; } tmp_check;
    union { __m256 v; float x[8]; } tmp_check_AVX;

    unstriped = malloc( sizeof(float) * ((Q*4)+1));  // Yes, these allocates are slow,but this is check code that won't be
    unstriped_AVX = malloc( sizeof(float) * ((Q_AVX*8)+1));  // compiled in production
    int q_temp;
    int z_temp;
    unstriped[0] = 0.;
    unstriped_AVX[0] = 0.;

    /* Line 1. M cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_MQ(bck, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
  for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_MQ(bck_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }

  // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("backward filter L-type M miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX), R = %d, R_AVX = %d\n", (i+1), q_temp, unstriped[q_temp], unstriped_AVX[q_temp], ox->R, ox->R_AVX);
    }
  }

  /* Line 2: I cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_IQ(bck, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_IQ(bck_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("backward filter L-type I miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", (i+1), q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }

 /* Line 3: D cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_DQ(bck, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_DQ(bck_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("backward filter L-type D miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", (i+1), q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }
  free(unstriped);
  free(unstriped_AVX); 

#endif      

  /* If there's any checkpointing, there's an L-1 row to fill now. */
  if (ox->Rb+ox->Rc > 0)
    {
#ifdef p7_build_SSE      
      /* Compute fwd[L-1] from last checkpoint, which we know is fwd[L-2] */
      dpp = (__m128 *) ox->dpf[ox->R0+ox->R-1];  /* fwd[L-2] values, already known        */
      fwd = (__m128 *) ox->dpf[ox->R0+ox->R];    /* get free row memory from top of stack */
      forward_row(dsq[i], om, dpp, fwd, Q);      /* calculate fwd[L-1]                    */
#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, fwd, "f2 X");
#endif

      /* Compute bck[L-1] from bck[L]. */
      xf  = (float *) (fwd + Q*p7C_NSCELLS);
      dpp = (__m128 *) ox->dpf[(i+1)%2]; 
      bck = (__m128 *) ox->dpf[i%2];             /* get space for bck[L-1]                */
      backward_row_main(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#ifdef p7_DEBUGGING
      ox->bcksc += logf(xf[p7C_SCALE]);
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, bck, "bck");
      if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, i, ox->bcksc); 
#endif
      /* And decode. */
      if ( (status = posterior_decode_row(ox, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
      dpp = bck;
      i--;			/* i is now L-2 if there's checkpointing; else it's L-1 */
#endif


#ifdef p7_build_AVX2     
      /* Compute fwd[L-1] from last checkpoint, which we know is fwd[L-2] */
      dpp_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX-1];  /* fwd[L-2] values, already known        */
      fwd_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX];    /* get free row memory from top of stack */
      forward_row_AVX(dsq[i_AVX], om, dpp_AVX, fwd_AVX, Q_AVX);      /* calculate fwd[L-1]                    */


      /* Compute bck[L-1] from bck[L]. */
      xf_AVX  = (float *) (fwd_AVX + Q_AVX*p7C_NSCELLS);
      dpp_AVX = (__m256 *) ox->dpf_AVX[(i_AVX+1)%2]; 
      bck_AVX = (__m256 *) ox->dpf_AVX[i_AVX%2];             /* get space for bck[L-1]                */
      backward_row_main_AVX(dsq[i_AVX+1], om, dpp_AVX, bck_AVX, Q_AVX, xf_AVX[p7C_SCALE]);

      /* And decode. */
      if ( (status = posterior_decode_row_AVX(ox, i_AVX, sm, sm_thresh, Tvalue_AVX)) != eslOK) return status;
      dpp_AVX = bck_AVX;
      i_AVX--;      /* i is now L-2 if there's checkpointing; else it's L-1 */
#endif

      #ifdef p7_build_check_AVX2  // Check that the computed rows match

    float *unstriped, *unstriped_AVX;
    union { __m128 v; float x[4]; } tmp_check;
    union { __m256 v; float x[8]; } tmp_check_AVX;

    unstriped = malloc( sizeof(float) * ((Q*4)+1));  // Yes, these allocates are slow,but this is check code that won't be
    unstriped_AVX = malloc( sizeof(float) * ((Q_AVX*8)+1));  // compiled in production
    int q_temp;
    int z_temp;
    unstriped[0] = 0.;
    unstriped_AVX[0] = 0.;

    /* Line 1. M cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_MQ(bck, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
  for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_MQ(bck_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }

  // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("backward filter main1-type M miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", (i+1), q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }

  /* Line 2: I cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_IQ(bck, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_IQ(bck_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("backward filter main1-type I miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", (i+1), q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }

 /* Line 3: D cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_DQ(bck, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_DQ(bck_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("backward filter main1-type D miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", (i+1), q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }
  free(unstriped);
  free(unstriped_AVX); 

#endif      
    }

  /* Main loop for checkpointed regions (b,c) */
   for (b = 2; b <= ox->Rb+ox->Rc; b++)
    {				/* i=L-2 as we enter here, and <dpp> is on bck[L-1] */
      w = (b <= ox->Rc ? b+1 : ox->Lb);

#ifdef p7_build_SSE
      /* We know current row i (r=R0+R-1) ends a block and is checkpointed in fwd. */
      ox->R--;
      fwd = (__m128 *) ox->dpf[ox->R0+ox->R];      /* pop checkpointed forward row off "stack" */
      xf  = (float *) (fwd + Q*p7C_NSCELLS);
      Tvalue = xf[p7C_C] * om->xf[p7O_C][p7O_MOVE];  /* i.e. scaled fwd[L] val at T state = scaled overall score */
      /* Calculate bck[i]; <dpp> is already bck[i+1] */
      bck = (__m128 *) ox->dpf[i%2];	    /* get available tmp memory for row     */
      backward_row_main(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#ifdef p7_DEBUGGING
      ox->bcksc += logf(xf[p7C_SCALE]);
      if (ox->do_dumping) { p7_checkptmx_DumpFBRow(ox, i, fwd, "f2 O");	p7_checkptmx_DumpFBRow(ox, i, bck, "bck"); }
      if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, i, ox->bcksc); 
#endif
      /* And decode checkpointed row i. */
      if ( (status = posterior_decode_row(ox, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
      
      /* The rest of the rows in the block weren't checkpointed.
       * Compute Forwards from last checkpoint ...
       */
      dpp = (__m128 *) ox->dpf[ox->R0+ox->R-1];       /* get last Fwd checkpoint. */
      for (i2 = i-w+1; i2 <= i-1; i2++)
	{
	  fwd = (__m128 *) ox->dpf[ox->R0+ox->R]; ox->R++;  /* push new forward row on "stack"     */
	  forward_row(dsq[i2], om, dpp, fwd, Q);
	  dpp = fwd;	  
	}

      /* ... and compute Backwards over the block we just calculated, while decoding. */
      dpp = bck;
      for (i2 = i-1; i2 >= i-w+1; i2--)
	{
	  ox->R--;
	  fwd = (__m128 *) ox->dpf[ox->R0+ox->R]; /* pop just-calculated forward row i2 off "stack" */
	  xf  = (float *) (fwd + Q*p7C_NSCELLS);
	  bck = (__m128 *) ox->dpf[i2%2];	  /* get available for calculating bck[i2]          */
	  backward_row_main(dsq[i2+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#ifdef p7_DEBUGGING
	  ox->bcksc += logf(xf[p7C_SCALE]);
	  if (ox->do_dumping) { p7_checkptmx_DumpFBRow(ox, i2, fwd, "f2 X"); p7_checkptmx_DumpFBRow(ox, i2, bck, "bck"); }
	  if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, i2, ox->bcksc); 
#endif

	  if ((status = posterior_decode_row(ox, i2, sm, sm_thresh, Tvalue)) != eslOK) return status;
	  dpp = bck;
	}
      i -= w;
#endif

#ifdef p7_build_AVX2
      /* We know current row i (r=R0+R-1) ends a block and is checkpointed in fwd. */
      ox->R_AVX--;
      fwd_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX];      /* pop checkpointed forward row off "stack" */
      xf_AVX  = (float *) (fwd_AVX + Q_AVX*p7C_NSCELLS);

      /* Calculate bck[i]; <dpp> is already bck[i+1] */
      bck_AVX = (__m256 *) ox->dpf_AVX[i_AVX%2];      /* get available tmp memory for row     */
      backward_row_main_AVX(dsq[i_AVX+1], om, dpp_AVX, bck_AVX, Q_AVX, xf_AVX[p7C_SCALE]);

      /* And decode checkpointed row i. */
      if ( (status = posterior_decode_row_AVX(ox, i_AVX, sm, sm_thresh, Tvalue_AVX)) != eslOK) return status;
      
      /* The rest of the rows in the block weren't checkpointed.
       * Compute Forwards from last checkpoint ...
       */
      dpp_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX-1];       /* get last Fwd checkpoint. */
      for (i2 = i_AVX-w+1; i2 <= i_AVX-1; i2++)
  {
    fwd_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX]; ox->R_AVX++;  /* push new forward row on "stack"     */
    forward_row_AVX(dsq[i2], om, dpp_AVX, fwd_AVX, Q_AVX);
    dpp_AVX = fwd_AVX;    
  }

      /* ... and compute Backwards over the block we just calculated, while decoding. */
      dpp_AVX = bck_AVX;
      for (i2 = i_AVX-1; i2 >= i_AVX-w+1; i2--)
  {
    ox->R_AVX--;
    fwd_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX]; /* pop just-calculated forward row i2 off "stack" */
    xf_AVX  = (float *) (fwd_AVX + Q_AVX*p7C_NSCELLS);
    bck_AVX = (__m256 *) ox->dpf_AVX[i2%2];   /* get available for calculating bck[i2]          */
    backward_row_main_AVX(dsq[i2+1], om, dpp_AVX, bck_AVX, Q_AVX, xf_AVX[p7C_SCALE]);

    if ((status = posterior_decode_row_AVX(ox, i2, sm, sm_thresh, Tvalue_AVX)) != eslOK) return status;
    dpp_AVX = bck_AVX;
  }
      i_AVX -= w;
#endif

#ifdef p7_build_check_AVX2  // Check that the computed rows match

    float *unstriped, *unstriped_AVX;
    union { __m128 v; float x[4]; } tmp_check;
    union { __m256 v; float x[8]; } tmp_check_AVX;

    unstriped = malloc( sizeof(float) * ((Q*4)+1));  // Yes, these allocates are slow,but this is check code that won't be
    unstriped_AVX = malloc( sizeof(float) * ((Q_AVX*8)+1));  // compiled in production
    int q_temp;
    int z_temp;
    unstriped[0] = 0.;
    unstriped_AVX[0] = 0.;

    /* Line 1. M cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_MQ(bck, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
  for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_MQ(bck_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }

  // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("backward filter main2-type M miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", i, q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }

  /* Line 2: I cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_IQ(bck, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_IQ(bck_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("backward filter main2-type I miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", i, q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }

 /* Line 3: D cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_DQ(bck, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_DQ(bck_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("backward filter main2-type D miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", i, q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }
  free(unstriped);
  free(unstriped_AVX); 

#endif      
    }
   /* now i=La as we leave the checkpointed regions; or i=L-1 if there was no checkpointing */
 #ifndef p7_build_SSE // need to set up the i variable for next loop
    #ifdef p7_build_AVX2
      int i = i_AVX;
    #else
      int i = i_AVX_512;
    #endif
#endif


   /* The uncheckpointed "a" region */
   for (; i >= 1; i--)
     {
#ifdef p7_build_SSE
       ox->R--; 
       fwd = (__m128 *) ox->dpf[ox->R0+ox->R]; /* pop off calculated row fwd[i]           */
       xf  = (float *) (fwd + Q*p7C_NSCELLS);
       bck = (__m128 *) ox->dpf[i%2];	       /* get open space for bck[i]               */
       backward_row_main(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#ifdef p7_DEBUGGING
       ox->bcksc += logf(xf[p7C_SCALE]);
       if (ox->do_dumping) { p7_checkptmx_DumpFBRow(ox, i, fwd, "f2 O"); p7_checkptmx_DumpFBRow(ox, i, bck, "bck"); }
       if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, i, ox->bcksc); 
#endif
       if ((status = posterior_decode_row(ox, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
       dpp = bck;
#endif
#ifdef p7_build_AVX2
     // Use regular i, not i_AVX here because it's the loop iterator variable
       ox->R_AVX--; 
       fwd_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX]; /* pop off calculated row fwd[i]           */
       xf_AVX  = (float *) (fwd_AVX + Q_AVX*p7C_NSCELLS);
       bck_AVX = (__m256 *) ox->dpf_AVX[i%2];        /* get open space for bck[i]               */
       backward_row_main_AVX(dsq[i+1], om, dpp_AVX, bck_AVX, Q_AVX, xf_AVX[p7C_SCALE]);
 
       if ((status = posterior_decode_row_AVX(ox, i, sm, sm_thresh, Tvalue_AVX)) != eslOK) return status;
       dpp_AVX = bck_AVX;
#endif

       #ifdef p7_build_check_AVX2  // Check that the computed rows match

    float *unstriped, *unstriped_AVX;
    union { __m128 v; float x[4]; } tmp_check;
    union { __m256 v; float x[8]; } tmp_check_AVX;

    unstriped = malloc( sizeof(float) * ((Q*4)+1));  // Yes, these allocates are slow,but this is check code that won't be
    unstriped_AVX = malloc( sizeof(float) * ((Q_AVX*8)+1));  // compiled in production
    int q_temp;
    int z_temp;
    unstriped[0] = 0.;
    unstriped_AVX[0] = 0.;

    /* Line 1. M cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_MQ(bck, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
  for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_MQ(bck_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }

  // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("backward filter main3-type M miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", i, q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }

  /* Line 2: I cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_IQ(bck, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_IQ(bck_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("backward filter main3-type I miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", i, q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }

 /* Line 3: D cells: unpack, unstripe, print */
  for (q_temp = 0; q_temp < Q; q_temp++) {
    tmp_check.v = P7C_DQ(bck, q_temp);
    for (z_temp = 0; z_temp < p7_VNF; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.x[z_temp];
  }
 for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
    tmp_check_AVX.v = P7C_DQ(bck_AVX, q_temp);
    for (z_temp = 0; z_temp < p7_VNF_AVX; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.x[z_temp];
  }
 // now, compare
  for(q_temp = 0; q_temp < Q * p7_VNF; q_temp++){
    if(fabs(unstriped[q_temp] - unstriped_AVX[q_temp]) > fabs(unstriped[q_temp]/ 100)){
      printf("backward filter main3-type D miss-match at row %d, position %d, %.8f (SSE) vs. %.8f (AVX)\n", i, q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
    }
  }
  free(unstriped);
  free(unstriped_AVX); 

#endif      
     }
#ifdef p7_DEBUGGING
   /* To get the backward score, we need to complete row 0 too, where
    * only the N and B states are reachable: B from Mk, and N from B.
    * Seeing the backward score match the forward score is useful in
    * debugging. But to do posterior decoding, which is all we need in 
    * production code, we can stop at row 1.
    */
   float xN;
   ox->R--;
   fwd = (__m128 *) ox->dpf[ox->R0];
   bck = (__m128 *) ox->dpf[i%2];	       
   xN = backward_row_zero(dsq[1], om, ox); 
   if (ox->do_dumping) { p7_checkptmx_DumpFBRow(ox, 0, fwd, "f2 O"); p7_checkptmx_DumpFBRow(ox, 0, bck, "bck"); }
   if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, 0, ox->bcksc); 
   if ((status = posterior_decode_row(ox, 0, sm, sm_thresh, Tvalue)) != eslOK) return status;
   ox->bcksc += xN;
#endif

   p7_sparsemask_Finish(sm);
   return eslOK;
}
/*----------- end forward/backward API calls --------------------*/



/*****************************************************************
 * 2. Internal functions: inlined recursions
 *****************************************************************/

/* forward_row()
 * 
 * <xi>   dsq[i]; residue on this row.
 * <om>   query model
 * <dpp>  ptr to previous row; may be a tmp row or a checkpointed row.
 * <dpc>  ptr to current row; ditto.
 * Q      number of vectors in <dpp>, <dpc>; ox->Qf.
 * 
 * Upon return, <dpc> contains scaled Forward values.
 * Returns log(scalevalue), for the caller to accumulate
 * as it accumulates the total Forward score.
 */
#ifdef p7_build_SSE 
static inline float
forward_row(ESL_DSQ xi, const P7_OPROFILE *om, const __m128 *dpp, __m128 *dpc, int Q)
{
  const    __m128 *rp   = om->rfv[xi];
  const    __m128 zerov = _mm_setzero_ps();
  const    __m128 *tp   = om->tfv;
  const    float  *xp   = (float *) (dpp + Q * p7C_NSCELLS);
  float           *xc   = (float *) (dpc + Q * p7C_NSCELLS);
  __m128          dcv   = _mm_setzero_ps();
  __m128          xEv   = _mm_setzero_ps();
  __m128          xBv   = _mm_set1_ps(xp[p7C_B]);
  __m128 mpv, dpv, ipv;
  __m128 sv;
  int    q;
  int    j;

  mpv = esl_sse_rightshift_ps(P7C_MQ(dpp, Q-1), zerov); 
  ipv = esl_sse_rightshift_ps(P7C_IQ(dpp, Q-1), zerov); 
  dpv = esl_sse_rightshift_ps(P7C_DQ(dpp, Q-1), zerov); 

  /* DP recursion for main states, all but the D->D path */
  for (q = 0; q < Q; q++)
    {
      /* Calculate M(i,q); hold it in tmp var <sv> */
      sv     =                _mm_mul_ps(xBv, *tp);  tp++; /* B->Mk    */
      sv     = _mm_add_ps(sv, _mm_mul_ps(mpv, *tp)); tp++; /* Mk-1->Mk */
      sv     = _mm_add_ps(sv, _mm_mul_ps(ipv, *tp)); tp++; /* Ik-1->Mk */
      sv     = _mm_add_ps(sv, _mm_mul_ps(dpv, *tp)); tp++; /* Dk-1->Dk */
      sv     = _mm_mul_ps(sv, *rp);                  rp++; /* e_Mk(x_i)*/
      xEv    = _mm_add_ps(xEv, sv);			   /* Mk->E    */

      /* Advance on previous row, picking up M,D,I values. */
      mpv = *dpp++;
      dpv = *dpp++;
      ipv = *dpp++;

      /* Delayed store of M,D */
      P7C_MQ(dpc, q) = sv;		
      P7C_DQ(dpc, q) = dcv;

      /* Partial calculation of *next* D(i,q+1); M->D only; delay storage, hold in dcv */
      dcv    = _mm_mul_ps(sv, *tp); tp++;

      /* Calculate and store I(i,q) */
      sv             =                _mm_mul_ps(mpv, *tp);  tp++;
      P7C_IQ(dpc, q) = _mm_add_ps(sv, _mm_mul_ps(ipv, *tp)); tp++;
    }

  /* Now the DD paths. We would rather not serialize them but 
   * in an accurate Forward calculation, we have few options.
   * dcv has carried through from end of q loop above; store it 
   * in first pass, we add M->D and D->D path into DMX.
   */ 
  /* We're almost certainly're obligated to do at least one complete 
   * DD path to be sure: 
   */
  dcv            = esl_sse_rightshift_ps(dcv, zerov);
  P7C_DQ(dpc, 0) = zerov;
  tp             = om->tfv + 7*Q;	/* set tp to start of the DD's */
  for (q = 0; q < Q; q++) 
    {
      P7C_DQ(dpc,q) = _mm_add_ps(dcv, P7C_DQ(dpc,q));	
      dcv           = _mm_mul_ps(P7C_DQ(dpc,q), *tp); tp++; /* extend DMO(q), so we include M->D and D->D paths */
    }
  
  /* now. on small models, it seems best (empirically) to just go
   * ahead and serialize. on large models, we can do a bit better,
   * by testing for when dcv (DD path) accrued to DMO(q) is below
   * machine epsilon for all q, in which case we know DMO(q) are all
   * at their final values. The tradeoff point is (empirically) somewhere around M=100,
   * at least on my desktop. We don't worry about the conditional here;
   * it's outside any inner loops.
   */
  if (om->M < 100)
    {			/* Fully serialized version */
      for (j = 1; j < 4; j++)
	{
	  dcv = esl_sse_rightshift_ps(dcv, zerov);
	  tp  = om->tfv + 7*Q;	/* reset tp to start of the DD's */
	  for (q = 0; q < Q; q++) 
	    { /* note, extend dcv, not DMO(q); only adding DD paths now */
	      P7C_DQ(dpc,q) = _mm_add_ps(dcv, P7C_DQ(dpc,q));	
	      dcv           = _mm_mul_ps(dcv, *tp);   tp++; 
	    }	    
	}
    } 
  else
    {			/* Slightly parallelized version, but which incurs some overhead */
      for (j = 1; j < 4; j++)
	{
	  register __m128 cv = zerov;	/* keeps track of whether any DD's change DMO(q) */

	  dcv = esl_sse_rightshift_ps(dcv, zerov);
	  tp  = om->tfv + 7*Q;	/* set tp to start of the DD's */
	  for (q = 0; q < Q; q++) 
	    { /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
	      sv            = _mm_add_ps(dcv, P7C_DQ(dpc,q));	
	      cv            = _mm_or_ps(cv, _mm_cmpgt_ps(sv, P7C_DQ(dpc,q))); 
	      P7C_DQ(dpc,q) = sv;	                                    /* store new DMO(q) */
	      dcv           = _mm_mul_ps(dcv, *tp);   tp++;                 /* note, extend dcv, not DMO(q) */
	    }	    
	  if (! _mm_movemask_ps(cv)) break; /* DD's didn't change any DMO(q)? Then done, break out. */
	}
    }
  
  /* Add Dk's to xEv */
  for (q = 0; q < Q; q++) xEv = _mm_add_ps(P7C_DQ(dpc,q), xEv);

  /* Specials, in order: E N JJ J B CC C */
  esl_sse_hsum_ps(xEv, &xc[p7C_E]);  
  xc[p7C_N]  =                                       xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_JJ] =                                       xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_J]  = xc[p7C_JJ]                          + xc[p7C_E] * om->xf[p7O_E][p7O_LOOP];
  xc[p7C_B]  = xc[p7C_N] * om->xf[p7O_N][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_J][p7O_MOVE];
  xc[p7C_CC] =                                       xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  xc[p7C_C]  = xc[p7C_CC]                          + xc[p7C_E] * om->xf[p7O_E][p7O_MOVE];

  /* Sparse rescaling. xE above threshold? Then trigger a rescaling event.            */
  if (xc[p7C_E] > 1.0e4)	/* that's a little less than e^10, ~10% of our dynamic range */
    {
      xc[p7C_N]  /= xc[p7C_E];
      xc[p7C_JJ] /= xc[p7C_E];
      xc[p7C_J]  /= xc[p7C_E];
      xc[p7C_B]  /= xc[p7C_E];
      xc[p7C_CC] /= xc[p7C_E];
      xc[p7C_C]  /= xc[p7C_E];
      xEv = _mm_set1_ps(1.0 / xc[p7C_E]);

      for (q = 0; q < Q; q++)
	{
	  P7C_MQ(dpc,q) = _mm_mul_ps(P7C_MQ(dpc,q), xEv);
	  P7C_DQ(dpc,q) = _mm_mul_ps(P7C_DQ(dpc,q), xEv);
	  P7C_IQ(dpc,q) = _mm_mul_ps(P7C_IQ(dpc,q), xEv);
	}

      xc[p7C_SCALE] = xc[p7C_E];
      xc[p7C_E]     = 1.0f;
    }
  else xc[p7C_SCALE] = 1.0f;

  return logf(xc[p7C_SCALE]);
}
#endif

#ifdef p7_build_AVX2
static inline float
forward_row_AVX(ESL_DSQ xi, const P7_OPROFILE *om, const __m256 *dpp, __m256 *dpc, int Q)
{
  const    __m256 *rp   = om->rfv_AVX[xi];
  const    __m256 zerov = _mm256_setzero_ps();
  const    __m256 *tp   = om->tfv_AVX;
  const    float  *xp   = (float *) (dpp + Q * p7C_NSCELLS);
  float           *xc   = (float *) (dpc + Q * p7C_NSCELLS);
  __m256          dcv   = _mm256_setzero_ps();
  __m256          xEv   = _mm256_setzero_ps();
  __m256          xBv   = _mm256_set1_ps(xp[p7C_B]);
  __m256 mpv, dpv, ipv;
  __m256 sv;
  int    q;
  int    j;

  mpv = esl_avx_leftshift_ps(P7C_MQ(dpp, Q-1)); 
  ipv = esl_avx_leftshift_ps(P7C_IQ(dpp, Q-1)); 
  dpv = esl_avx_leftshift_ps(P7C_DQ(dpp, Q-1)); 

  /* DP recursion for main states, all but the D->D path */
  for (q = 0; q < Q; q++)
    {
      /* Calculate M(i,q); hold it in tmp var <sv> */
      sv     =                _mm256_mul_ps(xBv, *tp);  tp++; /* B->Mk    */
      sv     = _mm256_add_ps(sv, _mm256_mul_ps(mpv, *tp)); tp++; /* Mk-1->Mk */
      sv     = _mm256_add_ps(sv, _mm256_mul_ps(ipv, *tp)); tp++; /* Ik-1->Mk */
      sv     = _mm256_add_ps(sv, _mm256_mul_ps(dpv, *tp)); tp++; /* Dk-1->Dk */
      sv     = _mm256_mul_ps(sv, *rp);                  rp++; /* e_Mk(x_i)*/
      xEv    = _mm256_add_ps(xEv, sv);        /* Mk->E    */

      /* Advance on previous row, picking up M,D,I values. */
      mpv = *dpp++;
      dpv = *dpp++;
      ipv = *dpp++;

      /* Delayed store of M,D */
      P7C_MQ(dpc, q) = sv;    
      P7C_DQ(dpc, q) = dcv;

      /* Partial calculation of *next* D(i,q+1); M->D only; delay storage, hold in dcv */
      dcv    = _mm256_mul_ps(sv, *tp); tp++;

      /* Calculate and store I(i,q) */
      sv             =                _mm256_mul_ps(mpv, *tp);  tp++;
      P7C_IQ(dpc, q) = _mm256_add_ps(sv, _mm256_mul_ps(ipv, *tp)); tp++;
    }

  /* Now the DD paths. We would rather not serialize them but 
   * in an accurate Forward calculation, we have few options.
   * dcv has carried through from end of q loop above; store it 
   * in first pass, we add M->D and D->D path into DMX.
   */ 
  /* We're almost certainly're obligated to do at least one complete 
   * DD path to be sure: 
   */
  dcv            = esl_avx_leftshift_ps(dcv);  // Function naming issue: The esl_sse_(right/left)shift_ps functions 
  // describe the direction of logical shift of the vectors.  The esl_avx_leftshift functions that I wrote for
  // other filters talk about the direction of the shift instruction, which is influenced by the little-endian 
  // nature of x86.  For consistency with myself, I'm sticking to function names that match the direction of the 
  // shift instruction
  P7C_DQ(dpc, 0) = zerov;
  tp             = om->tfv_AVX + 7*Q; /* set tp to start of the DD's */
  for (q = 0; q < Q; q++) 
    {
      P7C_DQ(dpc,q) = _mm256_add_ps(dcv, P7C_DQ(dpc,q)); 
      dcv           = _mm256_mul_ps(P7C_DQ(dpc,q), *tp); tp++; /* extend DMO(q), so we include M->D and D->D paths */
    }
  
  /* now. on small models, it seems best (empirically) to just go
   * ahead and serialize. on large models, we can do a bit better,
   * by testing for when dcv (DD path) accrued to DMO(q) is below
   * machine epsilon for all q, in which case we know DMO(q) are all
   * at their final values. The tradeoff point is (empirically) somewhere around M=100,
   * at least on my desktop. We don't worry about the conditional here;
   * it's outside any inner loops.
   */
  if (om->M < 100)
    {     /* Fully serialized version */
      for (j = 1; j < 8 /*was 4 */; j++)  /* What is this number?  Should it be vector length? Guess yes because of the 
      rightshift inside the loop*/
  { 
    dcv = esl_avx_leftshift_ps(dcv);
    tp  = om->tfv_AVX + 7*Q;  /* reset tp to start of the DD's */
    for (q = 0; q < Q; q++) 
      { /* note, extend dcv, not DMO(q); only adding DD paths now */
        P7C_DQ(dpc,q) = _mm256_add_ps(dcv, P7C_DQ(dpc,q)); 
        dcv           = _mm256_mul_ps(dcv, *tp);   tp++; 
      }     
  }
    } 
  else
    {     /* Slightly parallelized version, but which incurs some overhead */
      for (j = 1; j < 8; j++) /* ditto for this number, see previous j loop */
  {
    register __m256 cv = zerov; /* keeps track of whether any DD's change DMO(q) */

    dcv = esl_avx_leftshift_ps(dcv);
    tp  = om->tfv_AVX + 7*Q;  /* set tp to start of the DD's */
    for (q = 0; q < Q; q++) 
      { /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
        sv            = _mm256_add_ps(dcv, P7C_DQ(dpc,q)); 
        cv            = _mm256_or_ps(cv, _mm256_cmp_ps(sv, P7C_DQ(dpc,q), 14));   // 14 = code for compare greater than
        P7C_DQ(dpc,q) = sv;                                     /* store new DMO(q) */
        dcv           = _mm256_mul_ps(dcv, *tp);   tp++;                 /* note, extend dcv, not DMO(q) */
      }     
    if (! _mm256_movemask_ps(cv)) break; /* DD's didn't change any DMO(q)? Then done, break out. */
  }
    }
  
  /* Add Dk's to xEv */
  for (q = 0; q < Q; q++) xEv = _mm256_add_ps(P7C_DQ(dpc,q), xEv);

  /* Specials, in order: E N JJ J B CC C */
  /* xc, xp are pointers into vectors provided as inputs, so should be no problem with running multiple
  ISAs at same time */

  esl_avx_hsum_ps(xEv, &xc[p7C_E]);  
  xc[p7C_N]  =                                       xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_JJ] =                                       xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_J]  = xc[p7C_JJ]                          + xc[p7C_E] * om->xf[p7O_E][p7O_LOOP];
  xc[p7C_B]  = xc[p7C_N] * om->xf[p7O_N][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_J][p7O_MOVE];
  xc[p7C_CC] =                                       xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  xc[p7C_C]  = xc[p7C_CC]                          + xc[p7C_E] * om->xf[p7O_E][p7O_MOVE];

  /* Sparse rescaling. xE above threshold? Then trigger a rescaling event.            */
  if (xc[p7C_E] > 1.0e4)  /* that's a little less than e^10, ~10% of our dynamic range */
    {
      xc[p7C_N]  /= xc[p7C_E];
      xc[p7C_JJ] /= xc[p7C_E];
      xc[p7C_J]  /= xc[p7C_E];
      xc[p7C_B]  /= xc[p7C_E];
      xc[p7C_CC] /= xc[p7C_E];
      xc[p7C_C]  /= xc[p7C_E];
      xEv = _mm256_set1_ps(1.0 / xc[p7C_E]);

      for (q = 0; q < Q; q++)
  {
    P7C_MQ(dpc,q) = _mm256_mul_ps(P7C_MQ(dpc,q), xEv);
    P7C_DQ(dpc,q) = _mm256_mul_ps(P7C_DQ(dpc,q), xEv);
    P7C_IQ(dpc,q) = _mm256_mul_ps(P7C_IQ(dpc,q), xEv);
  }

      xc[p7C_SCALE] = xc[p7C_E];
      xc[p7C_E]     = 1.0f;
    }
  else xc[p7C_SCALE] = 1.0f;

  return logf(xc[p7C_SCALE]);
}
#endif
#ifdef p7_build_AVX512
static inline float
forward_row_AVX_512(ESL_DSQ xi, const P7_OPROFILE *om, const __m512 *dpp, __m512 *dpc, int Q)
{
  const    __m512 *rp   = om->rfv_AVX_512[xi];
  const    __m512 zerov = _mm512_setzero_ps();
  const    __m512 *tp   = om->tfv_AVX_512;
  const    float  *xp   = (float *) (dpp + Q * p7C_NSCELLS);
  float           *xc   = (float *) (dpc + Q * p7C_NSCELLS);
  __m512          dcv   = _mm512_setzero_ps();
  __m512          xEv   = _mm512_setzero_ps();
  __m512          xBv   = _mm512_set1_ps(xp[p7C_B]);
  __m512 mpv, dpv, ipv;
  __m512 sv;
  int    q;
  int    j;

  mpv = esl_avx_512_leftshift_ps(P7C_MQ(dpp, Q-1)); 
  ipv = esl_avx_512_leftshift_ps(P7C_IQ(dpp, Q-1)); 
  dpv = esl_avx_512_leftshift_ps(P7C_DQ(dpp, Q-1)); 

  /* DP recursion for main states, all but the D->D path */
  for (q = 0; q < Q; q++)
    {
      /* Calculate M(i,q); hold it in tmp var <sv> */
      sv     =                _mm512_mul_ps(xBv, *tp);  tp++; /* B->Mk    */
      sv     = _mm512_add_ps(sv, _mm512_mul_ps(mpv, *tp)); tp++; /* Mk-1->Mk */
      sv     = _mm512_add_ps(sv, _mm512_mul_ps(ipv, *tp)); tp++; /* Ik-1->Mk */
      sv     = _mm512_add_ps(sv, _mm512_mul_ps(dpv, *tp)); tp++; /* Dk-1->Dk */
      sv     = _mm512_mul_ps(sv, *rp);                  rp++; /* e_Mk(x_i)*/
      xEv    = _mm512_add_ps(xEv, sv);        /* Mk->E    */

      /* Advance on previous row, picking up M,D,I values. */
      mpv = *dpp++;
      dpv = *dpp++;
      ipv = *dpp++;

      /* Delayed store of M,D */
      P7C_MQ(dpc, q) = sv;    
      P7C_DQ(dpc, q) = dcv;

      /* Partial calculation of *next* D(i,q+1); M->D only; delay storage, hold in dcv */
      dcv    = _mm512_mul_ps(sv, *tp); tp++;

      /* Calculate and store I(i,q) */
      sv             =                _mm512_mul_ps(mpv, *tp);  tp++;
      P7C_IQ(dpc, q) = _mm512_add_ps(sv, _mm512_mul_ps(ipv, *tp)); tp++;
    }

  /* Now the DD paths. We would rather not serialize them but 
   * in an accurate Forward calculation, we have few options.
   * dcv has carried through from end of q loop above; store it 
   * in first pass, we add M->D and D->D path into DMX.
   */ 
  /* We're almost certainly're obligated to do at least one complete 
   * DD path to be sure: 
   */
  dcv            = esl_avx_512_leftshift_ps(dcv);  // Function naming issue: The esl_sse_(right/left)shift_ps functions 
  // describe the direction of logical shift of the vectors.  The esl_avx_leftshift functions that I wrote for
  // other filters talk about the direction of the shift instruction, which is influenced by the little-endian 
  // nature of x86.  For consistency with myself, I'm sticking to function names that match the direction of the 
  // shift instruction
  P7C_DQ(dpc, 0) = zerov;
  tp             = om->tfv_AVX_512 + 7*Q; /* set tp to start of the DD's */
  for (q = 0; q < Q; q++) 
    {
      P7C_DQ(dpc,q) = _mm512_add_ps(dcv, P7C_DQ(dpc,q)); 
      dcv           = _mm512_mul_ps(P7C_DQ(dpc,q), *tp); tp++; /* extend DMO(q), so we include M->D and D->D paths */
    }
  
  /* now. on small models, it seems best (empirically) to just go
   * ahead and serialize. on large models, we can do a bit better,
   * by testing for when dcv (DD path) accrued to DMO(q) is below
   * machine epsilon for all q, in which case we know DMO(q) are all
   * at their final values. The tradeoff point is (empirically) somewhere around M=100,
   * at least on my desktop. We don't worry about the conditional here;
   * it's outside any inner loops.
   */
  if (om->M < 100)
    {     /* Fully serialized version */
      for (j = 1; j < 16; j++)  
  { 
    dcv = esl_avx_512_leftshift_ps(dcv);
    tp  = om->tfv_AVX_512 + 7*Q;  /* reset tp to start of the DD's */
    for (q = 0; q < Q; q++) 
      { /* note, extend dcv, not DMO(q); only adding DD paths now */
        P7C_DQ(dpc,q) = _mm512_add_ps(dcv, P7C_DQ(dpc,q)); 
        dcv           = _mm512_mul_ps(dcv, *tp);   tp++; 
      }     
  }
    } 
  else
    {     /* Slightly parallelized version, but which incurs some overhead */
      for (j = 1; j < 16; j++) 
  {
    register __mmask16 cv = 0; /* keeps track of whether any DD's change DMO(q) */

    dcv = esl_avx_512_leftshift_ps(dcv);
    tp  = om->tfv_AVX_512 + 7*Q;  /* set tp to start of the DD's */
    for (q = 0; q < Q; q++) 
      { /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
        sv            = _mm512_add_ps(dcv, P7C_DQ(dpc,q)); 
        cv            = _mm512_kor(cv,_mm512_cmp_ps_mask(sv, P7C_DQ(dpc,q), 14));   // 14 = code for compare greater than
        P7C_DQ(dpc,q) = sv;                                     /* store new DMO(q) */
        dcv           = _mm512_mul_ps(dcv, *tp);   tp++;                 /* note, extend dcv, not DMO(q) */
      }     
    if (! cv) break; /* DD's didn't change any DMO(q)? Then done, break out. */
  }
    }
  
  /* Add Dk's to xEv */
  for (q = 0; q < Q; q++) xEv = _mm512_add_ps(P7C_DQ(dpc,q), xEv);

  /* Specials, in order: E N JJ J B CC C */
  /* xc, xp are pointers into vectors provided as inputs, so should be no problem with running multiple
  ISAs at same time */

  esl_avx_512_hsum_ps(xEv, &xc[p7C_E]);  
  xc[p7C_N]  =                                       xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_JJ] =                                       xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_J]  = xc[p7C_JJ]                          + xc[p7C_E] * om->xf[p7O_E][p7O_LOOP];
  xc[p7C_B]  = xc[p7C_N] * om->xf[p7O_N][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_J][p7O_MOVE];
  xc[p7C_CC] =                                       xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  xc[p7C_C]  = xc[p7C_CC]                          + xc[p7C_E] * om->xf[p7O_E][p7O_MOVE];

  /* Sparse rescaling. xE above threshold? Then trigger a rescaling event.            */
  if (xc[p7C_E] > 1.0e4)  /* that's a little less than e^10, ~10% of our dynamic range */
    {
      xc[p7C_N]  /= xc[p7C_E];
      xc[p7C_JJ] /= xc[p7C_E];
      xc[p7C_J]  /= xc[p7C_E];
      xc[p7C_B]  /= xc[p7C_E];
      xc[p7C_CC] /= xc[p7C_E];
      xc[p7C_C]  /= xc[p7C_E];
      xEv = _mm512_set1_ps(1.0 / xc[p7C_E]);

      for (q = 0; q < Q; q++)
  {
    P7C_MQ(dpc,q) = _mm512_mul_ps(P7C_MQ(dpc,q), xEv);
    P7C_DQ(dpc,q) = _mm512_mul_ps(P7C_DQ(dpc,q), xEv);
    P7C_IQ(dpc,q) = _mm512_mul_ps(P7C_IQ(dpc,q), xEv);
  }

      xc[p7C_SCALE] = xc[p7C_E];
      xc[p7C_E]     = 1.0f;
    }
  else xc[p7C_SCALE] = 1.0f;

  return logf(xc[p7C_SCALE]);
}
#endif

/* backward_row_main()
 * 
 * Backward calculation for rows 1..L-1. Rows 0, L are special cases.
 * 
 * <xi>        = residue dsq[i+1] on NEXT row
 * <om>        = query model
 * <dpp>       = 'previous' row i+1
 * <dpc>       = current row i
 * Q           = number of vectors in row; ox->Qf
 * scalefactor = scalefactor that Forward set for this row.
 *
 * Upon return, <dpc> contains backwards values. 
 *
 * Also, <dpp> is modified as a side effect, and is no longer
 * a valid backward row; caller should consider it destroyed. Its M
 * values now include the emission score contribution. This is
 * done for efficiency; better to add them in once and be done with
 * it.
 */
#ifdef p7_build_SSE 
static inline void
backward_row_main(ESL_DSQ xi, const P7_OPROFILE *om, __m128 *dpp, __m128 *dpc, int Q, float scalefactor)
{
  const __m128 *rp       = om->rfv[xi];			      /* emission scores on row i+1, for bck; xi = dsq[i+1]  */
  float       * const xc = (float *) (dpc + Q * p7C_NSCELLS); /* E N JJ J B CC C SCALE */
  const float * const xp = (float *) (dpp + Q * p7C_NSCELLS);
  const __m128 *tp, *tpdd;
  const __m128  zerov = _mm_set1_ps(0.0f);
  __m128        xBv;
  __m128       *dp;
  __m128        xEv;
  __m128        dcv, mcv, ipv, mpv;
  int           q;
  __m128 tmmv, timv, tdmv;				   /* copies of transition prob quads; a leftshift is needed as boundary cond */
 
  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. 
   * This invalidates <dpp> as a backwards row; its values are now
   * intermediates in the calculation of the current <dpc> row.
   */
  dp  = dpp;
  xBv = zerov;
  tp  = om->tfv;		/* on first transition vector */
  for (q = 0; q < Q; q++)
    {
      *dp = _mm_mul_ps(*dp, *rp); rp++;
      xBv = _mm_add_ps(xBv, _mm_mul_ps(*dp, *tp)); dp+= p7C_NSCELLS; tp += 7;
    }

  /* Specials. Dependencies dictate partial order C,CC,B < N,J,JJ < E */
  xc[p7C_C] = xc[p7C_CC] = xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  esl_sse_hsum_ps(xBv, &(xc[p7C_B]));
  xc[p7C_J] = xc[p7C_JJ] = xc[p7C_B] * om->xf[p7O_J][p7O_MOVE] + xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_N]              = xc[p7C_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_E]              = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_E][p7O_LOOP];

  /* Initialize for the row calculation */
  mpv  = esl_sse_leftshift_ps(*dpp,       zerov); /* [1 5 9 13] -> [5 9 13 x], M(i+1,k+1) * e(M_k+1, x_{i+1}) */
  tmmv = esl_sse_leftshift_ps(om->tfv[1], zerov);
  timv = esl_sse_leftshift_ps(om->tfv[2], zerov);
  tdmv = esl_sse_leftshift_ps(om->tfv[3], zerov);
  xEv  = _mm_set1_ps(xc[p7C_E]);
  tp   = om->tfv + 7*Q - 1;
  tpdd = tp + Q;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--)
    {
      ipv                 = P7C_IQ(dpp, q);
      P7C_IQ(dpc,q)       = _mm_add_ps( _mm_mul_ps(ipv, *tp),   _mm_mul_ps(mpv, timv)); tp--;   /* II,IM; I is done         */
      mcv                 = _mm_add_ps( _mm_mul_ps(ipv, *tp),   _mm_mul_ps(mpv, tmmv)); tp-=2;  /* MI,MM; ME,MD remain      */
      dcv                 = _mm_add_ps( _mm_mul_ps(dcv, *tpdd), _mm_mul_ps(mpv, tdmv)); tpdd--; /* DM and one segment of DD */

      P7C_DQ(dpc,q) = dcv = _mm_add_ps( xEv, dcv);
      P7C_MQ(dpc,q)       = _mm_add_ps( xEv, mcv);

      mpv  = P7C_MQ(dpp, q);
      tdmv = *tp; tp--;
      timv = *tp; tp--;
      tmmv = *tp; tp-=2;
    }
  backward_row_finish(om, dpc, Q, dcv);
  backward_row_rescale(xc, dpc, Q, scalefactor);
}
#endif

#ifdef p7_build_AVX2
static inline void
backward_row_main_AVX(ESL_DSQ xi, const P7_OPROFILE *om, __m256 *dpp, __m256 *dpc, int Q, float scalefactor)
{
  const __m256 *rp       = om->rfv_AVX[xi];           /* emission scores on row i+1, for bck; xi = dsq[i+1]  */
  float       * const xc = (float *) (dpc + Q * p7C_NSCELLS); /* E N JJ J B CC C SCALE */
  const float * const xp = (float *) (dpp + Q * p7C_NSCELLS);
  const __m256 *tp, *tpdd;
  const __m256  zerov = _mm256_setzero_ps();
  __m256        xBv;
  __m256       *dp;
  __m256        xEv;
  __m256        dcv, mcv, ipv, mpv;
  int           q;
  __m256 tmmv, timv, tdmv;           /* copies of transition prob quads; a leftshift is needed as boundary cond */
 
  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. 
   * This invalidates <dpp> as a backwards row; its values are now
   * intermediates in the calculation of the current <dpc> row.
   */
  dp  = dpp;
  xBv = zerov;
  tp  = om->tfv_AVX;    /* on first transition vector */
  for (q = 0; q < Q; q++)
    {
      *dp = _mm256_mul_ps(*dp, *rp); rp++;
      xBv = _mm256_add_ps(xBv, _mm256_mul_ps(*dp, *tp)); dp+= p7C_NSCELLS; tp += 7;
    }

  /* Specials. Dependencies dictate partial order C,CC,B < N,J,JJ < E */
  xc[p7C_C] = xc[p7C_CC] = xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  esl_avx_hsum_ps(xBv, &(xc[p7C_B]));
  xc[p7C_J] = xc[p7C_JJ] = xc[p7C_B] * om->xf[p7O_J][p7O_MOVE] + xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_N]              = xc[p7C_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_E]              = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_E][p7O_LOOP];

  /* Initialize for the row calculation */
  mpv  = esl_avx_rightshift_ps(*dpp      ); /* [1 5 9 13] -> [5 9 13 x], M(i+1,k+1) * e(M_k+1, x_{i+1}) */
  tmmv = esl_avx_rightshift_ps(om->tfv_AVX[1]);
  timv = esl_avx_rightshift_ps(om->tfv_AVX[2]);
  tdmv = esl_avx_rightshift_ps(om->tfv_AVX[3]);
  xEv  = _mm256_set1_ps(xc[p7C_E]);
  tp   = om->tfv_AVX + 7*Q - 1;
  tpdd = tp + Q;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--)
    {
      ipv                 = P7C_IQ(dpp, q);
      P7C_IQ(dpc,q)       = _mm256_add_ps( _mm256_mul_ps(ipv, *tp),   _mm256_mul_ps(mpv, timv)); tp--;   /* II,IM; I is done         */
      mcv                 = _mm256_add_ps( _mm256_mul_ps(ipv, *tp),   _mm256_mul_ps(mpv, tmmv)); tp-=2;  /* MI,MM; ME,MD remain      */
      dcv                 = _mm256_add_ps( _mm256_mul_ps(dcv, *tpdd), _mm256_mul_ps(mpv, tdmv)); tpdd--; /* DM and one segment of DD */

      P7C_DQ(dpc,q) = dcv = _mm256_add_ps( xEv, dcv);
      P7C_MQ(dpc,q)       = _mm256_add_ps( xEv, mcv);

      mpv  = P7C_MQ(dpp, q);
      tdmv = *tp; tp--;
      timv = *tp; tp--;
      tmmv = *tp; tp-=2;
    }
  backward_row_finish_AVX(om, dpc, Q, dcv);
  backward_row_rescale_AVX(xc, dpc, Q, scalefactor);
}
#endif

/* backward_row_L()
 * 
 * Backward calculation for row L; 
 * a special case because the matrix <ox> has no 'previous' row L+1.
 * Otherwise identical to backward_row_main().
 */
 #ifdef p7_build_SSE
static inline void
backward_row_L(const P7_OPROFILE *om,  __m128 *dpc, int Q, float scalefactor)
{
  const __m128  zerov = _mm_setzero_ps();
  float        *xc    = (float *) (dpc + Q * p7C_NSCELLS);
  const __m128 *tpdd;
  __m128       *dp;
  __m128       xEv, dcv;
  int          q;

  /* Backwards from T <- C,CC <- E;  all other specials unreachable, impossible on row L.
   * specials are stored in order E N JJ J B CC C.  
   */
  xc[p7C_C] = xc[p7C_CC] = om->xf[p7O_C][p7O_MOVE];
  xc[p7C_B] = xc[p7C_J] = xc[p7C_JJ] = xc[p7C_N] = 0.0;
  xc[p7C_E] = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE];

  xEv  = _mm_set1_ps(xc[p7C_E]);
  dp   = dpc + Q*p7C_NSCELLS - 1;
  tpdd = om->tfv + 8*Q - 1;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--) 
    {
      *dp--       = zerov;		                              /* I */
      *dp-- = dcv = _mm_add_ps(xEv, _mm_mul_ps(dcv, *tpdd)); tpdd--;  /* D */
      *dp--       = xEv; 	                                      /* M */
    }
  backward_row_finish(om, dpc, Q, dcv);
  backward_row_rescale(xc, dpc, Q, scalefactor);
}
#endif
#ifdef p7_build_AVX2
static inline void
backward_row_L_AVX(const P7_OPROFILE *om,  __m256 *dpc, int Q, float scalefactor)
{
  const __m256  zerov = _mm256_setzero_ps();
  float        *xc    = (float *) (dpc + Q * p7C_NSCELLS);
  const __m256 *tpdd;
  __m256       *dp;
  __m256       xEv, dcv;
  int          q;

  /* Backwards from T <- C,CC <- E;  all other specials unreachable, impossible on row L.
   * specials are stored in order E N JJ J B CC C.  
   */
  xc[p7C_C] = xc[p7C_CC] = om->xf[p7O_C][p7O_MOVE];
  xc[p7C_B] = xc[p7C_J] = xc[p7C_JJ] = xc[p7C_N] = 0.0;
  xc[p7C_E] = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE];

  xEv  = _mm256_set1_ps(xc[p7C_E]);
  dp   = dpc + Q*p7C_NSCELLS - 1;
  tpdd = om->tfv_AVX + 8*Q - 1;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--) 
    {
      *dp--       = zerov;                                  /* I */
      *dp-- = dcv = _mm256_add_ps(xEv, _mm256_mul_ps(dcv, *tpdd)); tpdd--;  /* D */
      *dp--       = xEv;                                        /* M */
    }
  backward_row_finish_AVX(om, dpc, Q, dcv);
  backward_row_rescale_AVX(xc, dpc, Q, scalefactor);
}
#endif


/* backward_row_finish()
 * 
 * Finishing DD and MD path contributions along the row.
 * Both L and main (1..L-1) recursions share this.
 * First segment of DD contributions has been calculated,
 * and up to three remain; then once D values are known,
 * we can add the M->D contribution to M.
 * 
 * Caller gets <dcv> from doing the first DD pass. 
 * <dcv> contains D scores for first vector, [1 5 9 13] in 
 * the example. We need it as an initialization condition;
 * we leftshift it to [5 9 13 *] and continue propagating
 * DD path.
 * 
 * The DD path calculation is the most annoying part of 
 * the otherwise gorgeous Farrar algorithm; it has to be
 * at least partially serialized, if not fully.
 * 
 * <om>  = query model
 * <dpc> = current row, partially calculated
 * Q     = width of rows in # vectors; ox->Qf
 * dcv   = the first D vector [1 5 9 13] from caller's earlier calculation
 */
#ifdef p7_build_SSE 
static inline void
backward_row_finish(const P7_OPROFILE *om, __m128 *dpc, int Q, __m128 dcv)
{
  const __m128 zerov = _mm_setzero_ps();
  const __m128 *tp;
  __m128       *dp;
  int           j,q;
  
  /* See notes on forward calculation: 
   * we have two options, either full serialization
   * or on long models, it becomes worthwhile to
   * check that all propagating DD contributions have
   * become negligble.
   */
  if (om->M < 100)
    { /* Full serialization */
      for (j = 1; j < 4; j++)
	{
	  dcv = esl_sse_leftshift_ps(dcv, zerov); /* [1 5 9 13] => [5 9 13 *]          */
	  tp  = om->tfv + 8*Q - 1;	          /* <*tp> now the [4 8 12 x] TDD quad */
	  dp  = dpc + Q*p7C_NSCELLS - 2;          /* init to point at D(i,q) vector    */
	  for (q = Q-1; q >= 0; q--)
	    {
	      dcv = _mm_mul_ps(dcv, *tp); tp--;
	      *dp = _mm_add_ps(*dp, dcv); dp -= p7C_NSCELLS;
	    }
	}
    }
  else
    { /* With check for early convergence */
      __m128 sv;
      __m128 cv;  /* keeps track of whether any DD addition changes DQ(q) value */
      for (j = 1; j < 4; j++)
	{
	  dcv = esl_sse_leftshift_ps(dcv, zerov);
	  tp  = om->tfv + 8*Q - 1;	
	  dp  = dpc + Q*p7C_NSCELLS - 2;
	  cv  = zerov;
	  for (q = Q-1; q >= 0; q--)
	    { /* using cmpgt below tests if DD changed any DMO(q) without conditional branch (i.e. no if) */
	      dcv  = _mm_mul_ps(dcv, *tp); tp--;
	      sv   = _mm_add_ps(*dp, dcv);
	      cv   = _mm_or_ps(cv, _mm_cmpgt_ps(sv, *dp)); /* if DD path changed DQ(dpc,q), cv bits know it now */
	      *dp  = sv; 
	      dp  -= p7C_NSCELLS;
	    }
	  if (! _mm_movemask_ps(cv)) break; /* if no DD path changed DQ(q) in this segment, then done, no more segments needed */
	}
    }

  /* Finally, M->D path contribution
   * these couldn't be added to M until we'd finished calculating D values on row.
   */
  dcv = esl_sse_leftshift_ps(P7C_DQ(dpc, 0), zerov);
  tp  = om->tfv + 7*Q - 3;	 
  dp  = dpc + (Q-1)*p7C_NSCELLS; 
  for (q = Q-1; q >= 0; q--)
    {
      *dp  = _mm_add_ps(*dp, _mm_mul_ps(dcv, *tp)); tp -= 7; 
      dcv  = *(dp+1);                               dp -= p7C_NSCELLS;
    }
}
#endif
#ifdef p7_build_AVX2
static inline void
backward_row_finish_AVX(const P7_OPROFILE *om, __m256 *dpc, int Q, __m256 dcv)
{
  const __m256 zerov = _mm256_setzero_ps();
  const __m256 *tp;
  __m256       *dp;
  int           j,q;
  
  /* See notes on forward calculation: 
   * we have two options, either full serialization
   * or on long models, it becomes worthwhile to
   * check that all propagating DD contributions have
   * become negligble.
   */
  if (om->M < 100)
    { /* Full serialization */
      for (j = 1; j < 8; j++)
  {
    dcv = esl_avx_rightshift_ps(dcv); /* [1 5 9 13] => [5 9 13 *]          */
    tp  = om->tfv_AVX + 8*Q - 1;            /* <*tp> now the [4 8 12 x] TDD quad */
    dp  = dpc + Q*p7C_NSCELLS - 2;          /* init to point at D(i,q) vector    */
    for (q = Q-1; q >= 0; q--)
      {
        dcv = _mm256_mul_ps(dcv, *tp); tp--;
        *dp = _mm256_add_ps(*dp, dcv); dp -= p7C_NSCELLS;
      }
  }
    }
  else
    { /* With check for early convergence */
      __m256 sv;
      __m256 cv;  /* keeps track of whether any DD addition changes DQ(q) value */
      for (j = 1; j < 8; j++)
  {
    dcv = esl_avx_rightshift_ps(dcv);
    tp  = om->tfv_AVX + 8*Q - 1;  
    dp  = dpc + Q*p7C_NSCELLS - 2;
    cv  = zerov;
    for (q = Q-1; q >= 0; q--)
      { /* using cmpgt below tests if DD changed any DMO(q) without conditional branch (i.e. no if) */
        dcv  = _mm256_mul_ps(dcv, *tp); tp--;
        sv   = _mm256_add_ps(*dp, dcv);
        cv   = _mm256_or_ps(cv, _mm256_cmp_ps(sv, *dp, 14)); /* if DD path changed DQ(dpc,q), cv bits know it now */
        *dp  = sv; 
        dp  -= p7C_NSCELLS;
      }
    if (! _mm256_movemask_ps(cv)) break; /* if no DD path changed DQ(q) in this segment, then done, no more segments needed */
  }
    }

  /* Finally, M->D path contribution
   * these couldn't be added to M until we'd finished calculating D values on row.
   */
  dcv = esl_avx_rightshift_ps(P7C_DQ(dpc, 0));
  tp  = om->tfv_AVX + 7*Q - 3;   
  dp  = dpc + (Q-1)*p7C_NSCELLS; 
  for (q = Q-1; q >= 0; q--)
    {
      *dp  = _mm256_add_ps(*dp, _mm256_mul_ps(dcv, *tp)); tp -= 7; 
      dcv  = *(dp+1);                               dp -= p7C_NSCELLS;
    }
}
#endif
/* backward_row_rescale()
 * 
 * Sparse rescaling, using the scalefactor that Forward set and 
 * Backward shares.
 * 
 * xc          - ptr to specials on row (floats)
 * dpc         - ptr to vectors for row      
 * Q           - number of vectors; ox->Qf
 * scalefactor - scalefactor that Forward set for this row
 *
 * Upon return, values in current row <dpc> have been rescaled.
 */
#ifdef p7_build_SSE
static inline void
backward_row_rescale(float *xc, __m128 *dpc, int Q, float scalefactor)
{
  if (scalefactor > 1.0f)
    {
      __m128  sv = _mm_set1_ps(1.0 / scalefactor);
      __m128 *dp = dpc;
      int     q;

      xc[p7C_E]  /= scalefactor;
      xc[p7C_N]  /= scalefactor;
      xc[p7C_JJ] /= scalefactor;
      xc[p7C_J]  /= scalefactor;
      xc[p7C_B]  /= scalefactor;
      xc[p7C_CC] /= scalefactor;
      xc[p7C_C]  /= scalefactor;

      for (q = 0; q < Q; q++) 
	{
	  *dp = _mm_mul_ps(*dp, sv); dp++; /* M */
	  *dp = _mm_mul_ps(*dp, sv); dp++; /* D */
	  *dp = _mm_mul_ps(*dp, sv); dp++; /* I */
	}
    }
  xc[p7C_SCALE] = scalefactor;
}
#endif
#ifdef p7_build_AVX2
static inline void
backward_row_rescale_AVX(float *xc, __m256 *dpc, int Q, float scalefactor)
{
  if (scalefactor > 1.0f)
    {
      __m256  sv = _mm256_set1_ps(1.0 / scalefactor);
      __m256 *dp = dpc;
      int     q;

      xc[p7C_E]  /= scalefactor;
      xc[p7C_N]  /= scalefactor;
      xc[p7C_JJ] /= scalefactor;
      xc[p7C_J]  /= scalefactor;
      xc[p7C_B]  /= scalefactor;
      xc[p7C_CC] /= scalefactor;
      xc[p7C_C]  /= scalefactor;

      for (q = 0; q < Q; q++) 
  {
    *dp = _mm256_mul_ps(*dp, sv); dp++; /* M */
    *dp = _mm256_mul_ps(*dp, sv); dp++; /* D */
    *dp = _mm256_mul_ps(*dp, sv); dp++; /* I */
  }
    }
  xc[p7C_SCALE] = scalefactor;
}
#endif

/* Only needed during development, for memory profiling; J10/29
 */
#if 0
static inline int
sse_countge(__m128 v, float thresh)
{
  union { __m128 v; float x[4]; } u;
  int z;
  int c = 0;
  u.v   = v;
  for (z = 0; z < 4; z++) if (u.x[z] >= thresh) c++;
  return c;
}
#endif


/* posterior_decode_row()
 *
 * In production code, we don't have to save results of posterior
 * decoding; we immediately use them to calculate sparse mask on the row.
 *
 * In debugging code, we save pp's, overwriting the forward row as we
 * go; when all pp's on the row have been calculated, if there's a
 * ox->pp matrix, we copy the pp values there.
 * 
 * Understanding how sparse rescaling interacts with posterior decoding:
 * Let F_i, B_i be the correct (unscaled) forward/backward terms for a
 * given cell on row i that we would've calculated, had we had
 * infinite numerical range. Let T be the overall Forward or Backward
 * score of the sequence comparison; we call it T, because we'll use
 * the Forward score in the terminal T state.
 * 
 * Let f_i, b_i, t be the forward/backward terms we actually calculate
 * using sparse rescaling. Let s_i be the scale terms on every row
 * i=1..L. (s_0 = 1). Sparse rescaling means that:
 *     f_i = F_i / \sum_j=1^i s_j 
 *     b_i = B_j / \sum_j=i^L s_j 
 *     t   = T   / \sum_j=1^L s_j 
 * Observe that the s_i term is present in both f_i and f_b.
 * 
 * The posterior decoding we want to calculate is:
 *     (F_i * B_i) / T
 * Substituting and cancelling scale terms gives:
 *     (f_i * b_i) * (s_i / t)
 * So in the code below: caller passes us <overall_sc> (i.e. t), and
 * we set a multiplicative "scaleterm" of (s_i / t).
 * 
 * Note that because the debugging version of posterior decoding
 * code overwrites the forward row, if the caller wants to 
 * dump, copy, or test anything on the forward row, it must do it
 * BEFORE calling posterior_decode_row().
 * 
 * The threshold for determining 'significant' posterior alignment
 * probability is passed as <sm_thresh> (typically
 * <p7_SPARSIFY_THRESH>, in <p7_config.h.in>). We take
 * advantage of the fact that if N/C/J emission postprobs exceed
 * <1.0 - sm_thresh>, we don't even need to look at the vectorized
 * row; no cell can exceed threshold.
 *
 * Can throw <eslEINVAL> on bad code (something awry in data structure initialization)
 *           <eslEMEM> on allocation failure in sparsemask
 */
#ifdef p7_build_SSE
static inline int
posterior_decode_row(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc)
{
  int             Q        = ox->Qf;
  __m128        *fwd       = (__m128 *) ox->dpf[ox->R0 + ox->R]; /* a calculated fwd row R has been popped off */
  const  __m128 *bck       = (__m128 *) ox->dpf[rowi%2];
  float         *xf        = (float *) (fwd + Q*p7C_NSCELLS);
  const  float  *xb        = (float *) (bck + Q*p7C_NSCELLS);
  const __m128   threshv   = _mm_set1_ps(sm_thresh);
  float          scaleterm = xf[p7C_SCALE] / overall_sc; /* see comments above, on how rescaling affects posterior decoding equations */
  const __m128   cv        = _mm_set1_ps(scaleterm);
  float  pnonhomology;
  __m128 mask;
  int    maskbits;		/* xxxx 4-bit mask for which cells 0..3 have passed threshold (if any) */
  __m128 pv;
  int    q,r;
  int    status;

  /* test to see if *any* cells can meet threshold, before wasting
   * time looking at them all.
   * 
   * Useful side effect: row 0 automatically fails this test (all pp
   * in S->N->B), so posterior_decode_row() can be called on row 0 (in
   * debugging code, we need to decode and store row zero specials),
   * without triggering contract check failures in p7_sparsemask_* API
   * functions that are checking for i=1..L.
   * 
   * This code block MAY NOT CHANGE the contents of fwd, bck vectors;
   * in debugging, we will use them again to recalculate and store the
   * decoding.
   */
  pnonhomology = (xf[p7C_N] * xb[p7C_N] + xf[p7C_JJ] * xb[p7C_JJ] + xf[p7C_CC] * xb[p7C_CC]) * scaleterm;
  if (pnonhomology <= 1.0f - sm_thresh)
    {
      if ((status = p7_sparsemask_StartRow(sm, rowi)) != eslOK) return status;
      for (q = Q-1; q >= 0; q--)	           // reverse, because SPARSEMASK is entirely in reversed order 
	{
	  pv       =                _mm_mul_ps(P7C_MQ(fwd, q), P7C_MQ(bck, q));
	  pv       = _mm_add_ps(pv, _mm_mul_ps(P7C_IQ(fwd, q), P7C_IQ(bck, q)));
	  pv       = _mm_add_ps(pv, _mm_mul_ps(P7C_DQ(fwd, q), P7C_DQ(bck, q)));
	  pv       = _mm_mul_ps(pv, cv);           // pv is now the posterior probability of elements q,r=0..3 
	  mask     = _mm_cmpge_ps(pv, threshv);    // mask now has all 0's in elems r that failed thresh; all 1's for r that passed 
	  maskbits = _mm_movemask_ps(mask);	   // maskbits is now something like 0100: 1's indicate which cell passed. 
	  
	  for (r = 0; r < p7_VNF; r++) 
	    if ( maskbits & (1<<r)) 
	      if ((status = p7_sparsemask_Add(sm, q, r)) != eslOK) return status;
	}
      if ((status = p7_sparsemask_FinishRow(sm)) != eslOK) return status;
    }

#ifdef p7_DEBUGGING
  xf[p7C_E]  = xf[p7C_E]  * xb[p7C_E]  * scaleterm;
  xf[p7C_N]  = (rowi == 0 ? 1.0f : xf[p7C_N]  * xb[p7C_N]  * scaleterm);
  xf[p7C_JJ] = xf[p7C_JJ] * xb[p7C_JJ] * scaleterm;
  xf[p7C_J]  = xf[p7C_J]  * xb[p7C_J]  * scaleterm;
  xf[p7C_B]  = xf[p7C_B]  * xb[p7C_B]  * scaleterm;
  xf[p7C_CC] = xf[p7C_CC] * xb[p7C_CC] * scaleterm;
  xf[p7C_C]  = xf[p7C_C]  * xb[p7C_C]  * scaleterm;

  for (q = 0; q < Q; q++)	
    {
      P7C_MQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(P7C_MQ(fwd, q), P7C_MQ(bck, q)));
      P7C_DQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(P7C_DQ(fwd, q), P7C_DQ(bck, q)));
      P7C_IQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(P7C_IQ(fwd, q), P7C_IQ(bck, q)));
    }

  if (ox->pp)  save_debug_row_pp(ox, fwd, rowi);
#endif
  return eslOK;
}
#endif
#ifdef p7_build_AVX2
static inline int
posterior_decode_row_AVX(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc)
{
  int             Q        = ox->Qf_AVX;
  __m256        *fwd       = (__m256 *) ox->dpf_AVX[ox->R0 + ox->R_AVX]; /* a calculated fwd row R has been popped off */
  const  __m256 *bck       = (__m256 *) ox->dpf_AVX[rowi%2];
  float         *xf        = (float *) (fwd + Q*p7C_NSCELLS);
  const  float  *xb        = (float *) (bck + Q*p7C_NSCELLS);
  const __m256   threshv   = _mm256_set1_ps(sm_thresh);
  float          scaleterm = xf[p7C_SCALE] / overall_sc; /* see comments above, on how rescaling affects posterior decoding equations */
  const __m256   cv        = _mm256_set1_ps(scaleterm);
  float  pnonhomology;
  __m256 mask;
  int    maskbits;    /* xxxx 4-bit mask for which cells 0..3 have passed threshold (if any) */
  __m256 pv;
  int    q,r;
  int    status;

  /* test to see if *any* cells can meet threshold, before wasting
   * time looking at them all.
   * 
   * Useful side effect: row 0 automatically fails this test (all pp
   * in S->N->B), so posterior_decode_row() can be called on row 0 (in
   * debugging code, we need to decode and store row zero specials),
   * without triggering contract check failures in p7_sparsemask_* API
   * functions that are checking for i=1..L.
   * 
   * This code block MAY NOT CHANGE the contents of fwd, bck vectors;
   * in debugging, we will use them again to recalculate and store the
   * decoding.
   */
  // Need to figure out how to put the sparsemask stuff back when we want to just run AVX 
  pnonhomology = (xf[p7C_N] * xb[p7C_N] + xf[p7C_JJ] * xb[p7C_JJ] + xf[p7C_CC] * xb[p7C_CC]) * scaleterm;
  if (pnonhomology <= 1.0f - sm_thresh)
    {  
   //   if ((status = p7_sparsemask_StartRow(sm, rowi)) != eslOK) return status;
      for (q = Q-1; q >= 0; q--)             // reverse, because SPARSEMASK is entirely in reversed order 
  {
    pv       =                _mm256_mul_ps(P7C_MQ(fwd, q), P7C_MQ(bck, q));
    pv       = _mm256_add_ps(pv, _mm256_mul_ps(P7C_IQ(fwd, q), P7C_IQ(bck, q)));
    pv       = _mm256_add_ps(pv, _mm256_mul_ps(P7C_DQ(fwd, q), P7C_DQ(bck, q)));
    pv       = _mm256_mul_ps(pv, cv);           // pv is now the posterior probability of elements q,r=0..3 
    mask     = _mm256_cmp_ps(pv, threshv, 13);  // 13 is magic number for >= 
     // mask now has all 0's in elems r that failed thresh; all 1's for r that passed 
    maskbits = _mm256_movemask_ps(mask);    // maskbits is now something like 0100: 1's indicate which cell passed. 
  
  /*  for (r = 0; r < p7_VNF; r++) 
      if ( maskbits & (1<<r)) 
     //   if ((status = p7_sparsemask_Add(sm, q, r)) != eslOK) return status; */
  }
    //  if ((status = p7_sparsemask_FinishRow(sm)) != eslOK) return status;
    }
#ifdef p7_DEBUGGING
  xf[p7C_E]  = xf[p7C_E]  * xb[p7C_E]  * scaleterm;
  xf[p7C_N]  = (rowi == 0 ? 1.0f : xf[p7C_N]  * xb[p7C_N]  * scaleterm);
  xf[p7C_JJ] = xf[p7C_JJ] * xb[p7C_JJ] * scaleterm;
  xf[p7C_J]  = xf[p7C_J]  * xb[p7C_J]  * scaleterm;
  xf[p7C_B]  = xf[p7C_B]  * xb[p7C_B]  * scaleterm;
  xf[p7C_CC] = xf[p7C_CC] * xb[p7C_CC] * scaleterm;
  xf[p7C_C]  = xf[p7C_C]  * xb[p7C_C]  * scaleterm;

  for (q = 0; q < Q; q++) 
    {
      P7C_MQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(P7C_MQ(fwd, q), P7C_MQ(bck, q)));
      P7C_DQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(P7C_DQ(fwd, q), P7C_DQ(bck, q)));
      P7C_IQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(P7C_IQ(fwd, q), P7C_IQ(bck, q)));
    }

  if (ox->pp)  save_debug_row_pp(ox, fwd, rowi);
#endif
  return eslOK;
}
#endif

/*------------------ end, inlined recursions -------------------*/




/*****************************************************************
 * 3. Debugging tools
 *****************************************************************/
#ifdef p7_DEBUGGING

/* backward_row_zero()
 * 
 * Slightly peculiar but true: in production code we don't
 * need to calculate backward row 0, because we're only doing
 * posterior decoding on residues 1..L.  We only need row 0
 * if we need a complete Backward calculation -- which happens
 * when we're in debugging mode, and we're going to test the
 * Backwards score or compare the Backwards matrix to the
 * reference implementation.
 * 
 * x1  - residue dsq[1] on row 1
 * om  - query model
 * ox  - DP matrix; dpf[0] is the current Backward row
 * 
 * Upon return, ox->dpf[0] (backward row 0) has been calculated.
 * Returns the log of the value in the N cell at row 0. When this 
 * is added to the sum of the logs of all the scalefactors (which
 * the caller has been accumulating in <ox->bcksc>), then
 * the <ox->bcksc> value is finished and equal to the Backwards
 * raw score in nats.
 */
static inline float
backward_row_zero(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox)
{
  int          Q     = ox->Qf;
  __m128       *dpc  = (__m128 *) ox->dpf[0];
  __m128       *dpp  = (__m128 *) ox->dpf[1];
  const __m128 *rp   = om->rfv[x1];
  const __m128 zerov = _mm_setzero_ps();
  float        *xc   = (float *) (dpc + Q * p7C_NSCELLS); /* special states on current row i  */
  float        *xp   = (float *) (dpp + Q * p7C_NSCELLS); /* special states on "previous" row i+1 */
  __m128       *dp;
  __m128       *tp;
  __m128        xBv  = zerov;
  int           q;

  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. */
  dp  = dpp;
  tp  = om->tfv;
  for (q = 0; q < Q; q++)
    {
      *dp = _mm_mul_ps(*dp, *rp); rp++;
      xBv = _mm_add_ps(xBv, _mm_mul_ps(*dp, *tp)); dp+= p7C_NSCELLS; tp += 7;
    }

  /* Only B,N,E will decode to nonzero posterior probability; C,J,E
   * can't be reached.  However, Backwards recursion gives C,J,E
   * nonzero values, because it's the Forward term that makes them
   * impossible. Though it's tempting to set these values to 0.0 (since
   * we know they're impossible), just do the calculation anyway;
   * by convention, all our DP algorithms (sparse, reference, and here)
   * do this, and we have unit tests comparing the values.
   */
  xc[p7C_C] = xc[p7C_CC] = xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  esl_sse_hsum_ps(xBv, &(xc[p7C_B]));
  xc[p7C_J] = xc[p7C_JJ] = xc[p7C_B] * om->xf[p7O_J][p7O_MOVE] + xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_N]              = xc[p7C_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_E]              = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_E][p7O_LOOP];
  
  /* not needed in production code: */
  for (q = 0; q < Q; q++)
    P7C_MQ(dpc, q) = P7C_IQ(dpc, q) = P7C_DQ(dpc, q) = zerov;

  return logf(xc[p7C_N]);
}


/* save_debug_row_pp()
 * 
 * Debugging only. Transfer posterior decoding values from a
 * vectorized row to appropriate row of <ox->pp>, as probabilities.
 *
 * Zero wherever there's no emission:
 *   all values on row 0 
 *   all D states
 *   all E,B states
 * Zero in row positions that can't be reached in some state:
 *   all insert states on row 1
 *   J(1),C(1)
 *   N(L), J(L) on row L
 * Zero wherever there's no state:
 *   D1, IM states
 */
static void
save_debug_row_pp(P7_CHECKPTMX *ox, __m128 *dpc, int i)
{
  union { __m128 v; float x[p7_VNF]; } u;
  int      Q  = ox->Qf;
  float  *xc  = (float *) (dpc + Q*p7C_NSCELLS);
  int     q,k,z,s;

  if (! ox->pp) return;
  
  P7R_XMX(ox->pp,i,p7R_E)  = xc[p7C_E];
  P7R_XMX(ox->pp,i,p7R_N)  = xc[p7C_N];
  P7R_XMX(ox->pp,i,p7R_J)  = xc[p7C_J];
  P7R_XMX(ox->pp,i,p7R_B)  = xc[p7C_B];
  P7R_XMX(ox->pp,i,p7R_L)  = xc[p7C_B];	/* all mass in local path */
  P7R_XMX(ox->pp,i,p7R_G)  = 0.0;	/* ... none in glocal     */
  P7R_XMX(ox->pp,i,p7R_C)  = xc[p7C_C];
  P7R_XMX(ox->pp,i,p7R_JJ) = xc[p7C_JJ];
  P7R_XMX(ox->pp,i,p7R_CC) = xc[p7C_CC];
  
  /* in a posterior decoding matrix (prob-space), all k=0 cells are 0.0 */
  for (s = 0; s < p7R_NSCELLS; s++) P7R_MX(ox->pp,i,0,s) = 0.0f;

  /* ... all mass is on local path for the filter, so all glocal cells are 0.0 */
  for (k =1; k <= ox->M; k++)
    {
      P7R_MX(ox->pp,i,k,p7R_MG) = 0.0f;
      P7R_MX(ox->pp,i,k,p7R_IG) = 0.0f;
      P7R_MX(ox->pp,i,k,p7R_DG) = 0.0f;
    }

  /* now the transfer from checkptmx decoding to local cells of refmx */
  for (q = 0; q < Q; q++)
    {
      u.v = P7C_MQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_ML) = u.x[z]; }
      u.v = P7C_DQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_DL) = u.x[z]; }
      u.v = P7C_IQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_IL) = u.x[z]; }
    }
}

/* save_debug_row_fb()
 * 
 * Debugging only. Transfer posterior decoding values (sparse scaled,
 * prob space) from a vectorized row, to appropriate row of <ox->fwd>
 * or <ox->bck> (log space, inclusive of partial sum of scalefactors);
 * <ox->fwd> and <ox->bck> should be identical (within numerical error
 * tolerance) to a reference implementation Forward/Backward in log
 * space.
 */
static void
save_debug_row_fb(P7_CHECKPTMX *ox, P7_REFMX *gx, __m128 *dpc, int i, float totscale)
{
  union { __m128 v; float x[p7_VNF]; } u;
  int      Q  = ox->Qf;
  float  *xc  = (float *) (dpc + Q*p7C_NSCELLS);
  int     q,k,z;

  if (! gx) return;
  
  P7R_XMX(gx,i,p7R_E)  = logf(xc[p7C_E]) + totscale;
  P7R_XMX(gx,i,p7R_N)  = logf(xc[p7C_N]) + totscale;
  P7R_XMX(gx,i,p7R_J)  = logf(xc[p7C_J]) + totscale;
  P7R_XMX(gx,i,p7R_B)  = logf(xc[p7C_B]) + totscale;
  P7R_XMX(gx,i,p7R_L)  = P7R_XMX(gx,i,p7R_B);         /* filter is local-mode. all mass assigned to local path */
  P7R_XMX(gx,i,p7R_G)  = -eslINFINITY;		      /* ... and no mass assigned to glocal path               */
  P7R_XMX(gx,i,p7R_C)  = logf(xc[p7C_C]) + totscale;
  P7R_XMX(gx,i,p7R_JJ) = -eslINFINITY;                /* JJ only saved in decoding, not fwd/bck */
  P7R_XMX(gx,i,p7R_CC) = -eslINFINITY;                /* ... CC, ditto                          */
  
  /* in P7_REFMX, all k=0 cells are initialized to -eslINFINITY;
   * set all glocal cells to -eslINFINITY too: */
  for (k =1; k <= ox->M; k++)
    {
      P7R_MX(gx,i,k,p7R_MG) = -eslINFINITY;
      P7R_MX(gx,i,k,p7R_IG) = -eslINFINITY;
      P7R_MX(gx,i,k,p7R_DG) = -eslINFINITY;
    }

  /* now the transfer from checkptmx (scaled prob-space) to local cells of refmx (log-space): */
  for (q = 0; q < Q; q++)
    {
      u.v = P7C_MQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_ML) = logf(u.x[z]) + totscale; }
      u.v = P7C_DQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_DL) = logf(u.x[z]) + totscale; }
      u.v = P7C_IQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_IL) = logf(u.x[z]) + totscale; }
    }
}
#endif
/*---------------- end, debugging tools -------------------------*/


/*****************************************************************
 * 4. Stats driver: memory requirement profiling 
 *****************************************************************/

/* In production pipeline, there are up to four sorts of DP matrices
 * in play:
 *    MSV   MSV filter
 *    VF    Viterbi filter
 *    CHK   checkpointed Forward/Backward
 *    SP    sparse DP
 * and the reference implementation has one, for our baseline:
 *    REF   reference DP implementation of everything
 * and, the sparse DP implementation also requires a sparse mask,
 * which is of on the same memory usage order as the DP matrix:
 *    SM    sparse DP mask
 * We count it separately here because we will typically have one
 * mask, for up to four (V/F/B/D) matrices in play.
 *    
 * MSV and VF are O(M).
 * CHK is O(M \sqrt L)
 * SP and SM are O(L) when a nonzero sparsity posterior threshold is used in the f/b filter.
 * REF is O(ML). 
 * 
 * The goal of this driver is to characterize the minimal memory
 * requirements for each, and in particular, for sparse DP.  MSV, VF,
 * CHK, and REF minimal requirements are completely determined by
 * dimensions M and L. The requirement of SP, though, is an empirical
 * question of how many sparse cells get included by the f/b filter.
 * This driver runs the f/b filter, using one query profile, against a
 * target seq db, and for each profile/seq comparison, it prints a
 * line summarizing all the minimal memory requirements.
 * 
 * Note that even though we're empirically characterizing the typical
 * minimum requirement of SP, it nonetheless does have a proven O(L)
 * memory complexity for a nonzero posterior prob threshold in f/b
 * filter.
 * 
 * Note also that this is characterizing the "minimal" requirement,
 * not the actual allocated size of the various matrices, which tend
 * to grow depending on previous requirements. We're not testing our
 * reallocation strategy here. We might use these statistics to help
 * refine that strategy.
 */
#ifdef p7FWDFILTER_STATS
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "stats driver, ForwardFilter()";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_FILTERMX    *fx      = NULL;
  P7_CHECKPTMX   *ox      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fraw, nullsc, fsc, msvsc;
  float           msvmem, vfmem, fbmem, smmem, spmem, refmem;
  double          P;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  /* Open sequence file for reading */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);               
  gm = p7_profile_Create(hmm->M, abc); 
  p7_profile_Config(gm, hmm, bg);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  
  /* Initially allocate matrices for a default sequence size, 500; 
   * we resize as needed for each individual target seq 
   */
  fx  = p7_filtermx_Create  (gm->M);
  ox  = p7_checkptmx_Create (gm->M, 500, ESL_MBYTES(32));  
  sm  = p7_sparsemask_Create(gm->M, 500);

  printf("# %-28s %-30s %9s %9s %9s %9s %9s %9s %9s %9s\n",
	 "target seq", "query profile", "score", "P-value", "MSV (KB)", "VF (KB)", "CHK (MB)", "SM (MB)", "SP (MB)", "REF (MB)");
  printf("# %-28s %-30s %9s %9s %9s %9s %9s %9s %9s %9s\n",
	 "----------------------------", "------------------------------",  "---------", "---------", "---------", "---------", "---------", "---------", "---------", "---------");
  
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_bg_SetLength(bg,            sq->n);

      p7_checkptmx_GrowTo (ox, om->M, sq->n); 
      p7_sparsemask_Reinit(sm, gm->M, sq->n);

      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);

      /* Filter insig hits, partially simulating the real pipeline */
      p7_MSVFilter(sq->dsq, sq->n, om, fx, &msvsc);
      msvsc = (msvsc - nullsc) / eslCONST_LOG2;
      P     =  esl_gumbel_surv(msvsc,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
      if (P > 0.02) goto NEXT_SEQ;

      p7_ForwardFilter (sq->dsq, sq->n, om, ox, &fraw);
      p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSIFY_THRESH);

      /* Calculate minimum memory requirements for each step */
      msvmem = (double) ( P7_NVB(om->M) * sizeof(__m128i))    / 1024.;  
      vfmem  = (double) p7_filtermx_MinSizeof(om->M)          / 1024.;
      fbmem  = (double) p7_checkptmx_MinSizeof(om->M, sq->n)  / 1024. / 1024.;
      smmem  = (double) p7_sparsemask_MinSizeof(sm)           / 1024. / 1024.;
      spmem  = (double) p7_sparsemx_MinSizeof(sm)             / 1024. / 1024.;
      refmem = (double) p7_refmx_MinSizeof(om->M, sq->n)      / 1024. / 1024.;

      fsc  =  (fraw-nullsc) / eslCONST_LOG2;
      P    = esl_exp_surv(fsc,   om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);

      printf("%-30s %-30s %9.2f %9.2g %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	     sq->name,
	     hmm->name,
	     fsc, P,
	     msvmem,
	     vfmem,
	     fbmem, 
	     smmem,
	     spmem,
	     refmem);
	     
    NEXT_SEQ:
      esl_sq_Reuse(sq);
      p7_sparsemask_Reuse(sm);
      p7_checkptmx_Reuse(ox);
      p7_filtermx_Reuse(fx);
    }

  printf("# SPARSEMASK: kmem reallocs: %d\n", sm->n_krealloc);
  printf("#             seg reallocs:  %d\n", sm->n_srealloc);
  printf("#             row reallocs:  %d\n", sm->n_rrealloc);

  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(ox);
  p7_filtermx_Destroy(fx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}


#endif
/*--------------- end, stats driver -----------------------------*/


/*****************************************************************
 * 5. Benchmark
 *****************************************************************/
/* Difference between gcc -g vs. icc -O3 is large! */

#ifdef p7FWDFILTER_BENCHMARK
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Forward",                         0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,   "2000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for checkpointed ForwardFilter()";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_CHECKPTMX   *ox      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);

  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_profile_SetLength(gm, L);

  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  ox  = p7_checkptmx_Create (om->M, L, ESL_MBYTES(32));
  sm  = p7_sparsemask_Create(om->M, L);

  /* Baseline time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_sparsemask_Reinit(sm, om->M, L);

      p7_ForwardFilter(dsq, L, om, ox, &sc);
      if (! esl_opt_GetBoolean(go, "-F")) 
	{
	  p7_BackwardFilter(dsq, L, om, ox, sm, p7_SPARSIFY_THRESH);
	  esl_vec_IReverse(sm->kmem, sm->kmem, sm->ncells);
	}

      p7_checkptmx_Reuse(ox);
      p7_sparsemask_Reuse(sm);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_checkptmx_Destroy(ox);
  p7_sparsemask_Destroy(sm);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7FWDFILTER_BENCHMARK*/
/*-------------------- end, benchmark ---------------------------*/


/*****************************************************************
 * 6. Unit tests
 *****************************************************************/
#ifdef p7FWDFILTER_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sqio.h"

#include "search/modelconfig.h"
#include "misc/logsum.h"
#include "misc/emit.h"

#include "dp_reference/reference_fwdback.h"
#include "dp_reference/reference_decoding.h"

/* Compare scores of Forward, Backward to those from the reference
 * implementation.
 * 
 */
static void
utest_scores(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, int64_t ramlimit)
{
  char         msg[]  = "fbfilter scores unit test failed";
  P7_HMM      *hmm    = NULL;
  P7_PROFILE  *gm     = NULL;
  P7_OPROFILE *om     = NULL;
  ESL_DSQ     *dsqmem = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_DSQ     *dsq    = NULL;
  int          tL     = 0;
  ESL_SQ      *sq     = esl_sq_CreateDigital(abc);
  P7_CHECKPTMX *ox    = p7_checkptmx_Create(M, L, ramlimit);
  P7_REFMX    *fwd    = p7_refmx_Create   (M, L);
  P7_REFMX    *bck    = p7_refmx_Create   (M, L);
  P7_REFMX    *pp     = p7_refmx_Create   (M, L);
  P7_SPARSEMASK *sm   = p7_sparsemask_Create(M, L);
  float        tol2   = ( p7_logsum_IsSlowExact() ? 0.001  : 0.1);   /* absolute agreement of reference (log-space) and vector (prob-space) depends on whether we're using LUT-based logsum() */
  float fsc1, fsc2;
  float bsc2;
#ifdef p7_DEBUGGING
  float        bsc1;
  float        tol1   = 0.0001;	                                     /* forward and backward scores from same implementation type should agree with high tolerance */
  float        ptol   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01);  /* posterior decoding values differ by no more than this */
#endif

#ifdef p7_DEBUGGING
  /* We set the debugging tools to record full pp, fwd, bck matrices
   * for comparison to reference implementation: 
   */
  ox->pp  = p7_refmx_Create(M,L);	
  ox->fwd = p7_refmx_Create(M,L);	
  ox->bck = p7_refmx_Create(M,L);	
#endif

  if ( p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om) != eslOK) esl_fatal(msg);
  /* note the <gm> is config'ed local-only; both <om>,<gm> have length model set to <L> */

  while (N--)
    {
      p7_profile_SetLength(gm, L);       /* because it may have been reset to tL by last emitted seq */
      p7_oprofile_ReconfigLength(om, L); 

      /* A mix of generated (homologous) and random (nonhomologous) sequences */
      if (esl_rnd_Roll(r, 2)) 
	{
	  esl_rsq_xfIID(r, bg->f, abc->K, L, dsqmem);  
	  dsq = dsqmem;  
	  tL = L;     
	  /* fixed-length random emission: length config on <gm>,<om> don't need to change  */
	}
      else  
	{
	  do {
	    esl_sq_Reuse(sq);
	    p7_ProfileEmit(r, hmm, gm, bg, sq, NULL);
	  } while (sq->n > L * 3); /* keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
	  dsq = sq->dsq; 
	  tL = sq->n; 
	  /* variable-length seq emission: length config on <gm>,<om> will change */
	}
	
      if ( p7_profile_SetLength(gm, tL)       != eslOK) esl_fatal(msg);
      if ( p7_oprofile_ReconfigLength(om, tL) != eslOK) esl_fatal(msg); 
      if ( p7_checkptmx_GrowTo(ox,  M, tL)    != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_Reinit(sm,M, tL)     != eslOK) esl_fatal(msg);

      p7_ForwardFilter (dsq, tL, om, ox, &fsc1);
      p7_BackwardFilter(dsq, tL, om, ox,  sm, p7_SPARSIFY_THRESH);

      p7_ReferenceForward (dsq, tL, gm, fwd,  &fsc2);
      p7_ReferenceBackward(dsq, tL, gm, bck,  &bsc2);
      p7_ReferenceDecoding(dsq, tL, gm, fwd, bck, pp);

#ifdef p7_DEBUGGING
      /* vector Forward and Backward scores should agree with high tolerance.
       * Backward score is only available in debugging mode 
       */
      bsc1 = ox->bcksc;
      if (fabs(fsc1-bsc1) > tol1)    esl_fatal(msg);
#endif

      /* reference scores should agree with tolerance depending on whether logsum compiled to use LUT or not */
      if (fabs(fsc2-bsc2) > tol2)    esl_fatal(msg);

      /* Reference and vector implementations should agree depending on logsum */
      if (fabs(fsc1-fsc2) > tol2) esl_fatal(msg);

#ifdef p7_DEBUGGING
      /* Compare all DP cell values to reference implementation,
       * in fwd, bck, and pp matrices. Note the need for CompareLocal()
       * for the Backward matrix comparison, because the zero B->G
       * transition isn't evaluated until the *end* of glocal paths;
       * thus Backward values along the glocal paths are finite for 
       * the reference implementation, whereas in the ForwardFilter
       * they're -inf by construction (ForwardFilter is local-only).
       */
      if (p7_refmx_Compare     (pp,  ox->pp,  ptol) != eslOK) esl_fatal(msg);
      if (p7_refmx_Compare     (fwd, ox->fwd, tol2) != eslOK) esl_fatal(msg);
      if (p7_refmx_CompareLocal(bck, ox->bck, tol2) != eslOK) esl_fatal(msg);
#endif
      esl_sq_Reuse(sq);
      p7_refmx_Reuse(fwd);
      p7_refmx_Reuse(bck);
      p7_refmx_Reuse(pp);
      p7_checkptmx_Reuse(ox);
      p7_sparsemask_Reuse(sm);
    }

  free(dsqmem);
  esl_sq_Destroy(sq);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(ox);
  p7_refmx_Destroy(fwd);
  p7_refmx_Destroy(bck);
  p7_refmx_Destroy(pp);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7FWDFILTER_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/


/*****************************************************************
 * 7. Test driver
 *****************************************************************/
#ifdef p7FWDFILTER_TESTDRIVE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for checkpointed vector SSE Forward, Backward implementations";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  P7_BG          *bg   = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  utest_scores(r, abc, bg, M, L, N,  ESL_MBYTES(32));   /* normal sized models              */
  utest_scores(r, abc, bg, M, L, N,  ESL_MBYTES(0));    /* zero memory: force checkpointing */
  utest_scores(r, abc, bg, 1, L, 10, ESL_MBYTES(32));   /* size 1 models                    */
  utest_scores(r, abc, bg, M, 1, 10, ESL_MBYTES(32));   /* size 1 sequences                 */

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  utest_scores(r, abc, bg, M, L, N,  ESL_MBYTES(32));  
  utest_scores(r, abc, bg, M, L, N,  ESL_MBYTES(0));   
  utest_scores(r, abc, bg, 1, L, 10, ESL_MBYTES(32));  
  utest_scores(r, abc, bg, M, 1, 10, ESL_MBYTES(32));  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);

  fprintf(stderr, "#  status = ok\n");
  return eslOK;
}
#endif /*p7FWDFILTER_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/


/*****************************************************************
 * 8. Example
 *****************************************************************/
#ifdef p7FWDFILTER_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "output in one line awkable format",                0 },
#ifdef p7_DEBUGGING
  { "-D",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump vector DP matrices for examination (verbose)",0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump recorded forward matrix",                     0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump recorded backward matrix",                    0 },
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump recorded posterior prob matrix",              0 },
  { "--diplot",  eslARG_OUTFILE, NULL, NULL, NULL,   NULL,  NULL, NULL, "save domain inference plot to <f>",                0 },
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example driver, ForwardFilter()";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_REFMX       *gx      = NULL;
  P7_CHECKPTMX   *ox      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fraw, nullsc, fsc, bsc;
  float           gfraw, gbraw, gfsc;
  float           gmem, cmem, bmem;
  double          P, gP;
  int             status;
#ifdef p7_DEBUGGING
  int             store_pp   = FALSE;
  char           *diplotfile = esl_opt_GetString(go, "--diplot");
  FILE           *difp       = NULL;
#endif

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  /* Open sequence file for reading */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

#ifdef p7_DEBUGGING
  if (diplotfile && (difp = fopen(diplotfile, "w")) == NULL) p7_Fail("couldn't open %s for writing", diplotfile);
  if (difp || esl_opt_GetBoolean(go, "-P")) store_pp = TRUE;
#endif

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);               
  gm = p7_profile_Create(hmm->M, abc); 
  p7_profile_Config(gm, hmm, bg);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  /* p7_oprofile_Dump(stdout, om);  */

  /* Initially allocate matrices for a default sequence size, 500; 
   * we resize as needed for each individual target seq 
   */
  ox  = p7_checkptmx_Create (gm->M, 500, ESL_MBYTES(32));  
  gx  = p7_refmx_Create     (gm->M, 500);
  sm  = p7_sparsemask_Create(gm->M, 500);
#ifdef p7_DEBUGGING
  /* When the p7_DEBUGGING flag is up, <ox> matrix has the ability to
   * record generic, complete <fwd>, <bck>, and <pp> matrices, for
   * comparison to reference implementation, even when checkpointing.
   */
  if (esl_opt_GetBoolean(go, "-D")) p7_checkptmx_SetDumpMode(ox, stdout, TRUE);
  if (esl_opt_GetBoolean(go, "-F")) ox->fwd = p7_refmx_Create(gm->M, 100);
  if (esl_opt_GetBoolean(go, "-B")) ox->bck = p7_refmx_Create(gm->M, 100);
  if (store_pp)                     ox->pp  = p7_refmx_Create(gm->M, 100);
#endif

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_bg_SetLength(bg,            sq->n);

      p7_checkptmx_GrowTo  (ox, om->M, sq->n); 
      p7_sparsemask_Reinit(sm, gm->M, sq->n);

      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);
    
      p7_ForwardFilter (sq->dsq, sq->n, om, ox, &fraw);
      p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSIFY_THRESH);

      p7_ReferenceForward (sq->dsq, sq->n, gm, gx, &gfraw);
      p7_ReferenceBackward(sq->dsq, sq->n, gm, gx, &gbraw);

      bsc = 0.0;		/* Backward score only available in debugging mode */
#ifdef p7_DEBUGGING
      if (esl_opt_GetBoolean(go, "-F")) p7_refmx_Dump(stdout, ox->fwd);
      if (esl_opt_GetBoolean(go, "-B")) p7_refmx_Dump(stdout, ox->bck);
      if (esl_opt_GetBoolean(go, "-P")) p7_refmx_Dump(stdout, ox->pp);
      if (difp)                         p7_refmx_PlotDomainInference(difp, ox->pp, 1, sq->n, NULL);
      bsc  =  (ox->bcksc-nullsc) / eslCONST_LOG2;
#endif

      fsc  =  (fraw-nullsc) / eslCONST_LOG2;
      gfsc = (gfraw-nullsc) / eslCONST_LOG2;
      P  = esl_exp_surv(fsc,   om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
      gP = esl_exp_surv(gfsc,  gm->evparam[p7_FTAU],  gm->evparam[p7_FLAMBDA]);

      gmem = (float) p7_refmx_Sizeof(gx)     / 1024 / 1024;
      cmem = (float) p7_checkptmx_Sizeof(ox) / 1024 / 1024;
      bmem = (float) sm->ncells * 6. * sizeof(float) / 1024 / 1024;

      if (esl_opt_GetBoolean(go, "-1")) 
	printf("%-30s\t%-20s\t%9.2g\t%7.4f\t%7.4f\t%9.2g\t%6.1f\t%6.2fM\t%6.2fM\t%6.2fM\n", sq->name, hmm->name, P, fsc, bsc, gP, gfsc, gmem, cmem, bmem);
      else
	{
	  printf("query model:               %s\n",        hmm->name);
	  printf("target sequence:           %s\n",        sq->name);
	  printf("fwd filter raw score:      %.4f nats\n", fraw);
#ifdef p7_DEBUGGING
	  printf("bck filter raw score:      %.4f nats\n", ox->bcksc);
#endif
	  printf("null score:                %.2f nats\n", nullsc);
	  printf("per-seq score:             %.2f bits\n", fsc);
	  printf("P-value:                   %g\n",        P);
	  printf("Reference fwd raw score:   %.2f nats\n", gfraw);
	  printf("Reference bck raw score:   %.2f nats\n", gbraw);
	  printf("Reference Fwd bit score:   %.2f bits\n", gfsc);
	  printf("Reference Forward P-val:   %g\n",        gP);
	  printf("RAM, f/b filter:           %.2fM\n",    cmem);
	  printf("RAM, generic:              %.2fM\n",    gmem);
	  printf("RAM, sparse:               %.2fM\n",    bmem);
	}

      esl_sq_Reuse(sq);
      p7_refmx_Reuse(gx);
      p7_sparsemask_Reuse(sm);
      p7_checkptmx_Reuse(ox);
    }

#ifdef p7_DEBUGGING
  if (difp) fclose(difp);
#endif
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(ox);
  p7_refmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7FWDFILTER_EXAMPLE*/
/*---------------------- end, example ---------------------------*/

/*****************************************************************
 * 9. Notes
 *****************************************************************/

/* [a.]  On debugging and testing methods.
 * 
 *    When compiled with <p7_DEBUGGING> (specifically, when
 *    fwdfilter.c and p7_checkptmx.[ch] are thus compiled), the
 *    <P7_CHECKPTMX> structure is augmented with additional fields for
 *    debugging dumps and unit test comparisons.
 *    
 *   %% Dumping vector matrices for examination
 *      The values from the vectorized <P7_CHECKPTMX> can be dumped
 *      during DP calculations by calling <p7_checkptmx_SetDumpMode()>
 *      on the object. Dumping has to happen during DP, not after,
 *      because of the way checkpointing discards rows as it goes (and
 *      for posterior decoding, the implementation never stores a row
 *      at all). Dumped rows are prefixed by a tag "f1 O", "f1 X", "f2
 *      O", "f2 X", or "bck", indicating backward (bck), first pass
 *      Forward (f1), second pass Forward (f2), checkpointed rows (O)
 *      that get saved and recalled, and discarded rows (X) that get
 *      recalculated in the second Forward pass.
 *      
 *      This capability is most useful for examining small DP matrices
 *      by hand; see fwdfilter_example, -D option.
 *      
 *   %% Saving matrices for comparison to reference
 *      With the <p7_DEBUGGING> compile flag, a caller may
 *      additionally provide an allocated <P7_REFMX> to the
 *      <P7_CHECKPTMX>, to enable storage of all DP matrix values in a
 *      form suitable for a <p7_refmx_Compare()> call against <P7_REFMX>
 *      DP matrices calculated by the reference implementation.
 *      Caller does something like <ox->fwd = p7_refmx_Create(M,L)> to
 *      save a Forward matrix, and/or analogously for <ox->bck> and/or
 *      <ox->pp> for Backward and Decoding.
 *      
 *      This capability is most useful for unit tests and automated
 *      comparison to the reference implementation. See utest_scores().
 *      
 *   %% High-precision comparison to reference
 *      Normally the reference implementation uses a table-driven
 *      log-sum-exp approximation (see misc/logsum.c), in order to do
 *      stable numerical calculations in log space. This introduces
 *      nonnegligible numerical error into DP calculations, so
 *      comparisons between a probability space vector implementation
 *      and the reference implementation must allow a large amount of
 *      numeric slop. At a cost of about 20x in speed, if the
 *      p7_LOGSUM_SLOWEXACT flag is compiled in, the p7_FLogsum()
 *      function uses the (more) exact calculation, allowing DP cells
 *      values to be compared more stringently.
 *      
 *      This capability is reflected in unit tests that set tolerances
 *      for floating-point comparison, after checking the flag with a
 *      <p7_logsum_IsSlowExact()> call. See utest_scores() for
 *      example.
 */

/* [b.] Running time, in theory and in practice.
 *
 *    Checkpointing requires more time in theory, but in practice you
 *    probably won't notice. The checkpointing method requires
 *    recalculation of Forward rows that weren't checkpointed, meaning
 *    up to two Forwards passes (plus one Backwards and one posterior
 *    decoding pass) over the DP matrix, rather than one each. First,
 *    using the Forward score as a filter minimizes this penalty:
 *    relatively few sequences pass on to the BackwardFilter
 *    call. Second, memory management in the checkpointed P7_CHECKPTMX
 *    structure uses "partial checkpointing" to minimize the use of
 *    checkpointing; for most comparisons, of all but the longest
 *    query/target combinations, DP calculations will fit in the
 *    available memory, and checkpointing is not invoked. And finally,
 *    the implementation here has been further optimized, such that
 *    it's actually slightly faster (with checkpointing) than the
 *    original HMMER3.0 implementation of Forward and Backward without
 *    checkpointing.
 * 
 */

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/




                                          
