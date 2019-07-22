/* H4_CHECKPTMX: checkpointed DP matrix for F/B filter
 * 
 * Independent of vector ISA. Do not add ISA-specific code.
 *
 * See h4_checkptmx.md for notes.
 * 
 * Contents:
 *    1. The <H4_CHECKPTMX> object
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_alloc.h"

#include "simdvec.h"
#include "h4_checkptmx.h"



/*****************************************************************
 * 1. The <H4_CHECKPTMX> object
 *****************************************************************/

/* minimum_rows()
 * 
 * Calculate the minimum number of rows <R> needed to checkpoint a
 * forward matrix for a sequence of length <L>, exclusive of other
 * constant row overhead (R0, consisting of two backwards rows and the
 * fwd[0] initial row).
 * 
 * This value is a double. If it has a fraction, a partial checkpoint
 * block ("b" region) is needed, as in this example code:
 *    Rbc  = minimum_rows(L);
 *    Rc   = floor(Rbc);
 *    Rb   = (Rbc > Rc ? 1 : 0);
 *    minR = (int) ceil(Rbc);    // or, Rc+Rb
 */
double 
minimum_rows(int L)
{
  return (sqrt(9. + 8. * (double) L) - 3.) / 2.;
}


/* checkpointed_rows()
 * 
 * Given <L> and <maxR>, where <maxR> is the number of rows we
 * have available for Ra+Rb+Rc (exclusive of R0), and given that we
 * already know that <maxR> is sufficient to hold the calculation
 * (maxR >= ceil(minimum_rows(L))), solve for the number of
 * checkpointed rows (Rb+Rc) we need. The value is a double, and if it
 * has a fractional part, then a partial checkpoint block is needed,
 * Rb==1.
 * 
 * This equation is obtained by solving 
 *     L = Ra + (Rbc +2)(Rbc+1)/2 - 1
 * for Rbc (i.e. Rb+Rc) using the quadratic equation, after
 * substitution Ra = (maxR - Rbc) to get Rbc in terms of L,maxR.
 */
double
checkpointed_rows(int L, int maxR)
{
  return (sqrt(1. + 8. * (double) (L - maxR)) - 1.) / 2.;
}


/* A "full" matrix is easy: Ra = La = L, using Ra+R0 <= maxR rows total. */
static void
set_full(H4_CHECKPTMX *cpx, int L)
{
  cpx->Ra     = L;
  cpx->Rb     = 0;
  cpx->Rc     = 0;
  cpx->La     = L;
  cpx->Lb     = 0;
  cpx->Lc     = 0;
}


/* If we can fit a checkpointed matrix into the <R_avail> rows we have,
 * then the trick is to use all of thems, making the "a" region (all
 * rows kept) as large as possible, to minimize computation. This
 * means solving a little quadratic equation for Rb+Rc given L and
 * R_avail - R0: see <checkpointed_rows()> for that solution.
 */
void
set_checkpointed(H4_CHECKPTMX *cpx, int L, int R_avail)
{
  double Rbc = checkpointed_rows(L, R_avail - cpx->R0);   // checkpointed_rows() calculation only deals with R=Ra+Rb+Rc part, not R0 fixed overhead.
  double Rc  = floor(Rbc);

  cpx->Rc     = (int) Rc;
  cpx->Rb     = (Rbc > Rc ? 1 : 0);
  cpx->Ra     = R_avail - cpx->Rb - cpx->Rc - cpx->R0;
  cpx->Lc     = ((cpx->Rc + 2) * (cpx->Rc + 1)) / 2 - 1;
  cpx->La     = cpx->Ra;
  cpx->Lb     = L - cpx->La - cpx->Lc;
}   



/* If we can't fit in R_avail rows, then we checkpoint
 * the entire matrix; R0+Rb+Rc > R_avail, so we'll
 * need to reallocate (redline) temporarily.
 */
void
set_redlined(H4_CHECKPTMX *cpx, int L, double minR)
{
  double Rc = floor(minR);

  cpx->Rc     = (int) Rc;
  cpx->Rb     = (minR > Rc ? 1 : 0);
  cpx->Ra     = 0;
  cpx->Lc     = ((cpx->Rc + 2) * (cpx->Rc + 1)) / 2 - 1;
  cpx->La     = 0;
  cpx->Lb     = L - cpx->La - cpx->Lc;
}


/* set_row_layout()
 *
 * Determines the size of the a,b,c regions ("all", "between",
 * "checkpointed") of rows in the DP matrix. 
 *
 * Caller must have already set <cpx->allocW> and <cpx->R0>; they are
 * needed here.
 * 
 * Upon return, we've set the R{0abc} and L{abc} fields in the <cpx>
 * structure.
 * 
 * <R_avail> is the maximum number of rows the caller *wants* to use,
 * either because a <ramlimit>'ed allocation fits that number of rows,
 * or because an existing matrix has that number of valid rows.  We
 * will exceed this for one comparison if absolutely necessary, but
 * the next <_Reuse()> call will bring the allocation back down.
 * 
 * So there's three possibilities:
 *  1. A full matrix fits into our recommended max memory use.
 *     In this case, don't checkpoint - just assign a row r to 
 *     every position i (La = Ra).
 *  2. A checkpointed matrix fits into our recommended memory.
 *     Make as much of the matrix uncheckpointed as we can,
 *     using every row in maxR.
 *  3. We can't satisfy the recommended max memory, even fully
 *     checkpointed. Make a fully checkpointed matrix, in which
 *     Ra=0 and R0+Rb+Rc will exceed R_avail. Caller will have to 
 *     allocate for a "redlined" larger space.
 */
static void
set_row_layout(H4_CHECKPTMX *cpx, int allocL, int R_avail)
{
  double Rbc      = minimum_rows(allocL);               
  int    minR_chk = cpx->R0 + (int) ceil(Rbc);            // min # rows we need for checkpointing
  int    minR_all = cpx->R0 + allocL;                     // min # rows we need for full matrix 

  if      (minR_all <= R_avail) set_full        (cpx, allocL);
  else if (minR_chk <= R_avail) set_checkpointed(cpx, allocL, R_avail);
  else                          set_redlined    (cpx, allocL, Rbc);
}


/* Function:  h4_checkptmx_Create()
 * Synopsis:  Allocate a new <H4_CHECKPTMX> object.
 *
 * Purpose:   Allocate a new <H4_CHECKPTMX> checkpointed, striped vector
 *            DP matrix sufficient for the Forward/Backward local
 *            decoding calculation for a query model of length <M> and
 *            a target sequence of length <L>.
 *            
 *            Try to keep the allocation under <redline> bytes in
 *            memory. For example, <redline=ESL_MBYTES(128)>, sets a
 *            recommended memory cap of 128 MiB. Allocation can exceed
 *            this, if even a fully checkpointed <MxL> comparison
 *            requires it -- but in this "redlined" case, any
 *            subsequent <h4_checkptmx_Reinit()> call will try to
 *            reallocate it back downwards to the <redline>.
 *            
 *            Choice of <redline> should take into account how many
 *            parallel threads there are, because each one will likely
 *            have its own <H4_CHECKPTMX> allocation.
 *            
 *            By design spec, <M> and <L> are $\leq$ 100K.
 *
 * Args:      M        - query profile size, consensus positions (<=100000)
 *            L        - target sequence length, residues (<=100000)
 *            ramlimit - recommended memory limit, bytes
 *
 * Returns:   ptr to new <H4_CHECKPTMX> object on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
H4_CHECKPTMX *
h4_checkptmx_Create(int M, int L, int64_t redline)
{
  H4_CHECKPTMX *cpx = NULL;
  int           maxR;
  int           r;
  int           status;
  
  /* Validity of integer variable ranges may depend on design spec: */
  ESL_DASSERT1( (M <= 100000) );       // design spec says  model length M <= 100000
  ESL_DASSERT1( (L <= 100000) );       //             ... and seq length L <= 100000
  ESL_DASSERT1( (L >  0) );
  ESL_DASSERT1( (M >  0) );

  /* Level 1 allocation: the structure itself */
  ESL_ALLOC(cpx, sizeof(H4_CHECKPTMX));
  cpx->dp_mem  = NULL;
  cpx->dpf     = NULL;

  /* Set checkpointed row layout: allocR, R{abc}, L{abc} fields */
  cpx->R0      = 3;       // fwd[0]; bck[prv,cur] 

  /* Calculate allocW, the allocated width (in cells, floats) of each vector row */
  /*                    maxQf        *  3 MID cells * floats/vector              */
  cpx->allocW  = H4_Q(M, h4_VMAX_FB) * h4C_NSCELLS  * h4_VMAX_FB;        
  cpx->allocW += ESL_UPROUND(h4C_NXCELLS, h4_VMAX_FB);    // specials are scalar. upround to maintain vector mem alignment
  cpx->redline = redline / sizeof(float);                 // redline, allocW, allocN are all in <float> units, internally.                             

  /* Caller said to keep within <redline> bytes; let's go ahead and use all of what we're allowed */
  maxR            = cpx->redline / cpx->allocW;
  set_row_layout(cpx, L, maxR);
  cpx->allocR      = cpx->R0 + cpx->Ra + cpx->Rb + cpx->Rc;
  cpx->validR      = cpx->allocR;

  ESL_DASSERT1(( cpx->allocW % h4_VMAX_FB == 0 )); /* verify that concat'ed rows will stay aligned */

  /* Level 2 allocations: row pointers and dp cell memory */
  cpx->allocN = cpx->allocR * cpx->allocW; 
  cpx->dp_mem = esl_alloc_aligned(sizeof(float) * cpx->allocN, h4_VALIGN);  // dp_mem is aligned vector memory
  if (cpx->dp_mem == NULL) { status = eslEMEM; goto ERROR; }

  ESL_DASSERT1(( cpx->allocN <= cpx->redline ));  // remember, internally, allocN and redline are both in units of # of floats

  ESL_ALLOC( cpx->dpf, sizeof(float **) * cpx->allocR);    
  for (r = 0; r < cpx->validR; r++)
    cpx->dpf[r] = cpx->dp_mem + (r * cpx->allocW);  
  
#if eslDEBUGLEVEL > 0
  cpx->do_dumping     = FALSE;
  cpx->dfp            = NULL;
  cpx->dump_maxpfx    = 5;       
  cpx->dump_width     = 9;
  cpx->dump_precision = 4;
  cpx->fwd            = NULL;
  cpx->bck            = NULL;
  cpx->pp             = NULL;
  cpx->bcksc          = 0.0f;
#endif

  cpx->M  = 0;
  cpx->L  = 0;
  cpx->R  = 0; 
  cpx->Q  = 0;
  cpx->Vf = 0;  // until striped data are valid. Gets set by vector DP code.
  return cpx;

 ERROR:
  h4_checkptmx_Destroy(cpx);
  return NULL;
}


/* Function:  h4_checkptmx_Reinit()
 * Synopsis:  Reinitialize checkpointed DP matrix for a new comparison.
 *
 * Purpose:   Given an existing checkpointed matrix structure <cpx>,
 *            and the dimensions <M> and <L> of a new comparison,
 *            reinitialize <cpx> for that comparison, reallocating
 *            if needed.
 *
 *            Essentially equivalent to _Create(), but while reusing
 *            previous space and minimizing expensive reallocations.
 *            Should only be called in a context similar to a
 *            _Create() call. For efficiency, old data are not copied
 *            to new internal allocations, so any existing data may be
 *            lost.
 *
 *            Usually <cpx> only grows. However, if <cpx> is redlined
 *            (over its recommended allocation) and the new problem
 *            size <M,L> can fit in the preset recommended allocation,
 *            then <cpx> is reallocated down to the smaller recommended
 *            size.
 *            
 * Args:      cpx   - existing checkpointed matrix
 *            M     - new query profile length
 *            L     - new target sequence length         
 * 
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> if an allocation fails. The state of <cpx> is
 *            now undefined, and the caller should not use it. 
 */
int
h4_checkptmx_Reinit(H4_CHECKPTMX *cpx, int M, int L)
{
  int     minR_chk      = (int) ceil(minimum_rows(L)) + cpx->R0; // minimum number of DP rows needed  
  int     reset_dp_ptrs = FALSE;
  int     maxR;
  int64_t W;                    // minimum row width needed, bytes 
  int     r;
  int     status;

  ESL_DASSERT1( (M <= 100000) );       // design spec says, model length M <= 100000 
  ESL_DASSERT1( (L <= 100000) );       //           ... and,  seq length L <= 100000 
  ESL_DASSERT1( (L >  0) );
  ESL_DASSERT1( (M >  0) );

  /* Reset data-specific stuff for a new comparison */
  cpx->M  = 0;
  cpx->L  = 0;
  cpx->R  = 0;
  cpx->Q  = 0;
  cpx->Vf = 0;

  /* If we're debugging and we have stored copies of any matrices,
   * reset and realloc them too.  Must do this first, because we have an early 
   * successful exit condition coming below.
   */
#if eslDEBUGLEVEL > 0
  if (cpx->fwd && (status = h4_refmx_GrowTo(cpx->fwd, M, L)) != eslOK) goto ERROR;
  if (cpx->bck && (status = h4_refmx_GrowTo(cpx->bck, M, L)) != eslOK) goto ERROR;
  if (cpx->pp  && (status = h4_refmx_GrowTo(cpx->pp,  M, L)) != eslOK) goto ERROR;
  cpx->bcksc = 0.0f;
#endif
 
  /* Calculate W, minimum row width needed by NEW allocation, in floats (cells) */
  /*         maxQf        *  3 MID cells * cells/vector                         */
  W  = H4_Q(M,h4_VMAX_FB) * h4C_NSCELLS  * h4_VMAX_FB;       // vector part, aligned for any vector ISA
  W += ESL_UPROUND(h4C_NXCELLS, h4_VMAX_FB);                 // specials are scalar. maintain vector mem alignment

  /* Are the current allocations already satisfactory? */
  if (W <= cpx->allocW && cpx->allocN <= cpx->redline)
    {
      if      (L + cpx->R0 <= cpx->validR) { set_full        (cpx, L);              return eslOK; }
      else if (minR_chk    <= cpx->validR) { set_checkpointed(cpx, L, cpx->validR); return eslOK; }
    }

  /* Do individual matrix rows need to be expanded, by resetting cpx->dpf ptrs into dp_mem? */
  if ( W > cpx->allocW) 
    {
      cpx->allocW    = W;
      cpx->validR    = cpx->allocN / cpx->allocW; // validR is still guaranteed to be <= allocR after this
      reset_dp_ptrs = TRUE;
    }

  /* Does matrix dp_mem need reallocation, either up or down? */
  maxR  = cpx->allocN / cpx->allocW;                              // max rows if we use up to the current allocation size.      
  if ( (cpx->allocN > cpx->redline && minR_chk <= maxR) ||        // we were redlined, and current alloc will work: so downsize 
       minR_chk > cpx->validR)                                    // not enough memory for needed rows: so upsize                   
    {
      set_row_layout(cpx, L, maxR); 
      cpx->validR = cpx->R0 + cpx->Ra + cpx->Rb + cpx->Rc;      // this may be > allocR now; we'll reallocate dp[] next, if so 
      cpx->allocN = cpx->validR * cpx->allocW;

      esl_alloc_free(cpx->dp_mem);
      cpx->dp_mem = esl_alloc_aligned(sizeof(float) * cpx->allocN, h4_VALIGN);
      if (cpx->dp_mem == NULL) { status = eslEMEM; goto ERROR; }

      reset_dp_ptrs = TRUE;
    }
  else  // current validR will suffice, either full or checkpointed; but we still need to calculate a layout 
    {
      if   (L+cpx->R0 <= cpx->validR) set_full(cpx, L); 
      else                            set_checkpointed(cpx, L, cpx->validR);
    }
  
  /* Does the array of row ptrs need reallocation? */
  if (cpx->validR > cpx->allocR)
    {
      ESL_REALLOC(cpx->dpf, sizeof(float *) * cpx->validR);
      cpx->allocR    = cpx->validR;
      reset_dp_ptrs = TRUE;
    }

  /* Do the row ptrs need to be reset? */
  if (reset_dp_ptrs)
    {
      for (r = 0; r < cpx->validR; r++)
        cpx->dpf[r] = cpx->dp_mem + (r * cpx->allocW);
    }

  return eslOK;

 ERROR:
  return status;
}



/* Function:  h4_checkptmx_Sizeof()
 * Synopsis:  Returns size of checkpointed vector DP matrix, in bytes.
 * 
 * Purpose:   Returns the size of the checkpointed vector DP matrix
 *            in bytes. 
 *            
 *            If code has been compiled in debugging mode, the
 *            returned size includes a negligible amount of extra
 *            space for debugging fields in the structure (about 5
 *            ints, 4 pointers, and a float - around 56 bytes). The
 *            returned size does not include the use of any full
 *            Forward, Backward, or decoding matrices in the debugging
 *            part of the structure. This is because when we're in
 *            debugging mode asking about memory usage, we're usually
 *            interested in the estimated usage of the production
 *            code, because we're optimizing some parameter choices
 *            for example.
 */
size_t
h4_checkptmx_Sizeof(const H4_CHECKPTMX *cpx)
{
  size_t n = sizeof(H4_CHECKPTMX);
  n += cpx->allocN  * sizeof(float);     // neglects alignment overhead of up to h4_VALIGN bytes
  n += cpx->allocR  * sizeof(float *);     
  return n;
}

/* Function:  h4_checkptmx_MinSizeof()
 * Synopsis:  Returns minimum required size of a <H4_CHECKPTMX>, in bytes.
 *
 * Purpose:   Calculate and return the minimal required size, in bytes,
 *            of a checkpointed f/b matrix, for a comparison of a profile
 *            of length <M> to a sequence of length <L>.
 *            
 *            Does not require having an actual DP matrix allocated.
 *            We use this function when planning/profiling memory
 *            allocation strategies.
 */
size_t
h4_checkptmx_MinSizeof(int M, int L)
{
  size_t n    = sizeof(H4_CHECKPTMX);
  int    Q    = H4_Q(M,h4_VMAX_FB);               // number of vectors needed
  int    minR = 3 + (int) ceil(minimum_rows(L));  // 3 = Ra, 2 rows for backwards, 1 for fwd[0]

  n += minR * (Q * h4C_NSCELLS * h4_VMAX_FB);         // each row: Q supervectors. Each supervector: h4C_NSCELLS=3 vectors MDI. Each vector: h4_VMAX_FB cells
  n += minR * (ESL_UPROUND(h4C_NXCELLS, h4_VMAX_FB)); // dp_mem, specials: maintaining vector memory alignment 
  n += minR * sizeof(float *);                        // dpf[] row ptrs
  return n;
}



/* Function:  h4_checkptmx_Destroy()
 * Synopsis:  Frees a <H4_CHECKPTMX>.
 *
 * Purpose:   Free the <H4_CHECKPTMX> <cpx>. <cpx> may be <NULL>,
 *            or incompletely allocated.
 */
void
h4_checkptmx_Destroy(H4_CHECKPTMX *cpx)
{
  if (cpx)
    {
      if (cpx->dp_mem) esl_alloc_free(cpx->dp_mem);
      if (cpx->dpf)    free(cpx->dpf);

#if eslDEBUGLEVEL > 0
      if (cpx->fwd)    h4_refmx_Destroy(cpx->fwd);
      if (cpx->bck)    h4_refmx_Destroy(cpx->bck);
      if (cpx->pp)     h4_refmx_Destroy(cpx->pp);
#endif
    }
}
/*--------------- end, H4_CHECKPTMX object -----------------------*/

