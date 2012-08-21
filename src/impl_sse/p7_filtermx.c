/* Implementation of P7_FILTERMX: checkpointed, striped vector DP matrix.
 * Shared by all the filter implementations (SSV, MSV, Vit, Fwd, Bck).
 * 
 * Contents:
 *    1. API for the P7_FILTERMX object
 *    2. Debugging, development routines.
 *    3. Internal routines.
 *    4. Copyright and license information.
 */
#include "p7_config.h"

#include <stdlib.h>

#include "easel.h"

#include "hmmer.h"
#include "p7_filtermx.h"

static void set_row_layout  (P7_FILTERMX *ox, int allocL, int maxR); 
static void set_full        (P7_FILTERMX *ox, int L);
static void set_checkpointed(P7_FILTERMX *ox, int L, int R);
static void set_redlined    (P7_FILTERMX *ox, int L, double minR);

static double minimum_rows     (int L);
static double checkpointed_rows(int L, int R);


/*****************************************************************
 * 1. API for the <P7_FILTERMX> object
 *****************************************************************/

/* Function:  p7_filtermx_Create()
 * Synopsis:  Allocate a new <P7_FILTERMX> object.
 *
 * Purpose:   Allocate a new <P7_FILTERMX> checkpointed, striped vector
 *            DP matrix sufficient for all necessary filter
 *            calculations (from SSV to Backwards) for a query model
 *            of up to length <M> and a target sequence of up to
 *            length <L>.
 *            
 *            Try to keep the allocation within <ramlimit> bytes in
 *            memory.  For example, <ramlimit=ESL_MBYTES(128)>, sets a
 *            recommended memory limit of 128 MiB. Allocation can
 *            exceed this, if even a fully checkpointed <MxL>
 *            comparison requires it -- but in this case, any
 *            subsequent <p7_filtermx_GrowTo()> call that attempts to
 *            reuse the matrix will try to reallocated it back
 *            downwards to the <ramlimit>.
 *            
 *            Choice of <ramlimit> should take into account how many
 *            parallel threads there are, because each one will likely
 *            have its own <P7_FILTERMX> allocation.
 *            
 *            By design spec, <M> and <L> are $\leq$ 100K.
 *            
 *            When the <P7_FILTERMX> is only being used for
 *            <MSVFilter()>, <VitFilter()>, or some other single-row
 *            DP implementation (as opposed to a checkpointed matrix),
 *            it suffices to call <p7_filtermx_Create(M, 0, 0)>.  Even
 *            with <L=0>, a minimum of 3 fwd/bck rows are allocated,
 *            which suffices for any single-row DP routine.
 *            
 *
 * Args:      M        - query profile size, consensus positions (<=100000)
 *            L        - target sequence length, residues (<=100000)
 *            ramlimit - recommended memory limit, bytes
 *
 * Returns:   ptr to new <P7_FILTERMX> object on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_FILTERMX *
p7_filtermx_Create(int M, int L, int64_t ramlimit)
{
  P7_FILTERMX *ox = NULL;
  int          maxR;
  int          r;
  int          status;
  
  /* Validity of integer variable ranges may depend on design spec:                  */
  ESL_DASSERT1( (M <= 100000) );       /* design spec says, model length M <= 100000 */
  ESL_DASSERT1( (L <= 100000) );       /*           ... and,  seq length L <= 100000 */
  ESL_DASSERT1( (L >  0) );
  ESL_DASSERT1( (M >  0) );

  /* Level 1 allocation: the structure itself */
  ESL_ALLOC(ox, sizeof(P7_FILTERMX));
  ox->dp_mem  = NULL;
  ox->dpf     = NULL;

  /* Set checkpointed row layout: allocR, R{abc}, L{abc} fields */
  ox->R0          = 3;	                                                   /* fwd[0]; bck[prv,cur] */
  ox->allocW      = sizeof(float) * P7F_NQF(M) * p7F_NSCELLS * p7_VNF;	   /* accounts for main vector part of the row       */
  ox->allocW     += ESL_UPROUND(sizeof(float) * p7F_NXCELLS, p7_VALIGN);  /* plus specials (must maintain memory alignment) */
  ox->ramlimit    = ramlimit;
  maxR            = (int) (ox->ramlimit / ox->allocW); 
  set_row_layout(ox, L, maxR);
  ox->allocR      = ox->R0 + ox->Ra + ox->Rb + ox->Rc;
  ox->validR      = ox->allocR;

  ESL_DASSERT1( (ox->allocW % p7_VALIGN == 0) ); /* verify alignment */

  /* Level 2 allocations: row pointers and dp cell memory */
  ox->nalloc = ox->allocR * ox->allocW;
  ESL_ALLOC( ox->dp_mem, ox->nalloc + (p7_VALIGN-1));    /* (p7_VALIGN-1) because we'll hand-align memory */
  ESL_ALLOC( ox->dpf,    sizeof(float *) * ox->allocR);

  ox->dpf[0] = (char *) ( ((uintptr_t) ox->dp_mem + p7_VALIGN - 1) & p7_VALIMASK); /* hand memory alignment */
  for (r = 1; r < ox->validR; r++)
    ox->dpf[r] = ox->dpf[0] + r * ox->allocW;

#ifdef p7_DEBUGGING
  ox->do_dumping     = FALSE;
  ox->dfp            = NULL;
  ox->dump_maxpfx    = 5;	
  ox->dump_width     = 9;
  ox->dump_precision = 4;
  ox->dump_flags     = p7_DEFAULT;
  ox->fwd            = NULL;
  ox->bck            = NULL;
  ox->pp             = NULL;
  ox->bcksc          = 0.0f;
#endif

  ox->M  = 0;
  ox->L  = 0;
  ox->R  = 0;
  ox->Qf = 0;
  return ox;

 ERROR:
  p7_filtermx_Destroy(ox);
  return NULL;
}

/* Function:  p7_filtermx_GrowTo()
 * Synopsis:  Resize checkpointed DP matrix for new seq/model comparison.
 *
 * Purpose:   Given an existing checkpointed matrix structure <ox>,
 *            and the dimensions <M> and <L> of a new comparison,
 *            reallocate and reinitialize <ox>.
 *
 *            Essentially the same as free'ing the previous matrix and
 *            creating a new one -- but minimizes expensive memory
 *            allocation/reallocation calls.
 *            
 *            Usually <ox> only grows. The exception is if <ox> is
 *            redlined (over its recommended allocation) and the new
 *            problem size <M,L> can fit in the preset recommended
 *            allocation, then <ox> is reallocated down to the smaller
 *            recommended size.
 *            
 * Args:      ox    - existing checkpointed matrix
 *            M     - new query profile length
 *            L     - new target sequence length         
 * 
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> if an allocation fails. The state of <ox> is
 *            now undefined, and the caller should not use it. 
 */
int
p7_filtermx_GrowTo(P7_FILTERMX *ox, int M, int L)
{
  int     minR_chk      = (int) ceil(minimum_rows(L)) + ox->R0; /* minimum number of DP rows needed  */
  int     reset_dp_ptrs = FALSE;
  int     maxR;
  int64_t W;			/* minimum row width needed, bytes */
  int     r;
  int     status;

  /* Validity of integer variable ranges may depend on design spec:                  */
  ESL_DASSERT1( (M <= 100000) );       /* design spec says, model length M <= 100000 */
  ESL_DASSERT1( (L <= 100000) );       /*           ... and,  seq length L <= 100000 */
  ESL_DASSERT1( (L >  0) );
  ESL_DASSERT1( (M >  0) );

  /* If we're debugging and we have stored copies of any matrices,
   * grow them too.  Must do this first, because we have an early exit
   * condition coming below.
   */
#ifdef p7_DEBUGGING
  if (ox->fwd && (status = p7_refmx_GrowTo(ox->fwd, M, L)) != eslOK) goto ERROR;
  if (ox->bck && (status = p7_refmx_GrowTo(ox->bck, M, L)) != eslOK) goto ERROR;
  if (ox->pp  && (status = p7_refmx_GrowTo(ox->pp,  M, L)) != eslOK) goto ERROR;
#endif

  /* Calculate W, the minimum row width needed, in bytes */
  W  = sizeof(float) * P7F_NQF(M) * p7F_NSCELLS * p7_VNF;    /* vector part of row (MDI)     */
  W += ESL_UPROUND(sizeof(float) * p7F_NXCELLS, p7_VALIGN);  /* float part of row (specials); must maintain p7_VALIGN-byte alignment */

  /* Are current allocations satisfactory ? */
  if (W <= ox->allocW && ox->nalloc <= ox->ramlimit)
    {
      if      (L + ox->R0 <= ox->validR) { set_full        (ox, L);             return eslOK; }
      else if (minR_chk   <= ox->validR) { set_checkpointed(ox, L, ox->validR); return eslOK; }
    }

  /* Do individual matrix rows need to expand? */
  if ( W > ox->allocW) 
    {
      ox->allocW    = W;
      ox->validR    = (int) (ox->nalloc / ox->allocW); /* validR must be <= allocR */
      reset_dp_ptrs = TRUE;
    }

  /* Does matrix dp_mem need reallocation, either up or down? */
  maxR  = (int) (ox->nalloc / ox->allocW);                      /* max rows if we use up to the recommended allocation size.      */
  if ( (ox->nalloc > ox->ramlimit && minR_chk <= maxR) ||       /* we were redlined, and recommended alloc will work: so downsize */
       minR_chk > ox->validR)				        /* not enough memory for needed rows: so upsize                   */
    {
      set_row_layout(ox, L, maxR); 
      ox->validR = ox->R0 + ox->Ra + ox->Rb + ox->Rc;   /* this may be > allocR now; we'll reallocate dp[] next, if so     */
      ox->nalloc = ox->validR * ox->allocW;
      ESL_REALLOC(ox->dp_mem, ox->nalloc + (p7_VALIGN-1)); /* (p7_VALIGN-1) because we will manually align dpf ptrs into dp_mem */
      reset_dp_ptrs = TRUE;
    }
  else  /* current validR will suffice, either full or checkpointed; we still need to calculate a layout */
    {
      if   (L+ox->R0 <= ox->validR) set_full(ox, L); 
      else                          set_checkpointed(ox, L, ox->validR);
    }
  
  /* Does the array of row ptrs need reallocation? */
  if (ox->validR > ox->allocR)
    {
      ESL_REALLOC(ox->dpf, sizeof(float *) * ox->validR);
      ox->allocR    = ox->validR;
      reset_dp_ptrs = TRUE;
    }

  /* Do the row ptrs need to be reset? */
  if (reset_dp_ptrs)
    {
      ox->dpf[0] = (char *) ( ( (uintptr_t) ox->dp_mem + p7_VALIGN - 1) & p7_VALIMASK); /* vectors must be aligned on p7_VALIGN-byte boundary */
      for (r = 1; r < ox->validR; r++)
	ox->dpf[r] = ox->dpf[0] + (r * ox->allocW);
    }

  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_filtermx_Sizeof()
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
p7_filtermx_Sizeof(const P7_FILTERMX *ox)
{
  size_t n = sizeof(P7_FILTERMX);
  n += ox->nalloc + (p7_VALIGN-1);	          /* +15 because of manual alignment */
  n += ox->allocR  * sizeof(float *);	  
  return n;
}

/* Function:  p7_filtermx_MinSizeof()
 * Synopsis:  Returns minimum required size of a <P7_FILTERMX>, in bytes.
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
p7_filtermx_MinSizeof(int M, int L)
{
  size_t n    = sizeof(P7_FILTERMX);
  int    Q    = P7F_NQF(M);                       // number of vectors needed
  int    minR = 3 + (int) ceil(minimum_rows(L));  // 3 = Ra, 2 rows for backwards, 1 for fwd[0]
  
  n += p7_VALIGN-1;                                                  // dp_mem has to be hand-aligned for vectors
  n += minR * (sizeof(float) * p7_VNF * Q * p7F_NSCELLS);            // dp_mem, main: QR supercells; each has p7F_NSCELLS=3 cells, MID; each cell is __m128 vector of four floats (p7_VNF=4 * float)
  n += minR * (ESL_UPROUND(sizeof(float) * p7F_NXCELLS, p7_VALIGN)); // dp_mem, specials: maintaining vector memory alignment 
  n += minR * sizeof(float *);                                       // dpf[] row ptrs
  return n;
}


/* Function:  p7_filtermx_Reuse()
 * Synopsis:  Recycle a checkpointed vector DP matrix.
 *
 * Purpose:   Resets the checkpointed vector DP matrix <ox> for reuse,
 *            minimizing free/malloc wastefulness. All information
 *            specific to the DP problem we just computed is
 *            reinitialized. All allocations (and information about
 *            those allocations) are preserved.
 *            
 *            Caller will still need to call <p7_filtermx_GrowTo()>
 *            before each new DP, to be sure that the allocations are
 *            sufficient, and checkpointed rows are laid out.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_filtermx_Reuse(P7_FILTERMX *ox)
{
  int status;

  ox->M  = 0;
  ox->L  = 0;
  ox->R  = 0;
  ox->Qf = 0;

#ifdef p7_DEBUGGING
  if (ox->fwd && (status = p7_refmx_Reuse(ox->fwd)) != eslOK) return status;
  if (ox->bck && (status = p7_refmx_Reuse(ox->bck)) != eslOK) return status;
  if (ox->pp  && (status = p7_refmx_Reuse(ox->pp))  != eslOK) return status;
  ox->bcksc = 0.0f;
#endif

  return eslOK;
}


/* Function:  p7_filtermx_Destroy()
 * Synopsis:  Frees a <P7_FILTERMX>.
 *
 * Purpose:   Free the <P7_FILTERMX> <ox>. <ox> may be <NULL>,
 *            or incompletely allocated.
 */
void
p7_filtermx_Destroy(P7_FILTERMX *ox)
{
 if (ox) {
   if (ox->dp_mem) free(ox->dp_mem);
   if (ox->dpf)    free(ox->dpf);
#ifdef p7_DEBUGGING
   if (ox->fwd)    p7_refmx_Destroy(ox->fwd);
   if (ox->bck)    p7_refmx_Destroy(ox->bck);
   if (ox->pp)     p7_refmx_Destroy(ox->pp);
#endif
   free(ox);
 }
}
/*--------------- end, P7_FILTERMX object -----------------------*/



/*****************************************************************
 * 2. Debugging, development routines
 *****************************************************************/

/* Function:  p7_filtermx_SetDumpMode()
 * Synopsis:  Toggle dump mode flag in a P7_FILTERMX.
 *
 * Purpose:   Toggles whether DP matrix rows will be dumped for examination
 *            during <p7_ForwardFilter()>, <p7_BackwardFilter()>.
 *            Dumping has to be done row by row, on the fly during the
 *            DP calculations, not afterwards, because these routines
 *            are memory efficient (checkpointed) and they do not save
 *            all their rows.
 * 
 *            To turn on dumping, <truefalse> is <TRUE>, and <dfp> is an
 *            open <FILE> pointer for dump output (such as <stderr>).
 *            
 *            To turn off dumping, <truefalse> is <FALSE>, and <dfp>
 *            should be <NULL>.
 *            
 *            For this to have effect, <p7_DEBUGGING> compile-time
 *            flag must be on. If it is not, no dumping will occur,
 *            and this call is a no-op.
 *
 * Args:      ox        - DP matrix object to set
 *            dfp       - open FILE * for debugging output, or NULL
 *            truefalse - TRUE to set dump mode, or FALSE to unset
 *
 * Returns:   <eslOK> on success.
 */
int
p7_filtermx_SetDumpMode(P7_FILTERMX *ox, FILE *dfp, int truefalse)
{
#ifdef p7_DEBUGGING
  ox->do_dumping    = truefalse;
  ox->dfp           = dfp;
#endif
  return eslOK;
}


#ifdef p7_DEBUGGING

/* Function:  p7_filtermx_DecodeX()
 * Synopsis:  Convert a special X statecode to a string.
 */
char *
p7_filtermx_DecodeX(enum p7f_xcells_e xcode)
{
  switch (xcode) {
  case p7F_E:     return "E"; 
  case p7F_N:     return "N"; 
  case p7F_JJ:    return "JJ"; 
  case p7F_J:     return "J"; 
  case p7F_B:     return "B"; 
  case p7F_CC:    return "CC"; 
  case p7F_C:     return "C"; 
  case p7F_SCALE: return "SCALE"; 
  default:        return "?";
  }
}

/* Function:  p7_filtermx_DumpFBHeader()
 * Synopsis:  Prints the header of the fwd/bck dumps.
 *
 * Purpose:   Print the header that accompanies 
 *            <p7_filtermx_DumpFBRow()>.
 */
int
p7_filtermx_DumpFBHeader(P7_FILTERMX *ox)
{
  int maxpfx = ox->dump_maxpfx;
  int width  = ox->dump_width;
  int M      = ox->M;
  int k;

  fprintf(ox->dfp, "%*s", maxpfx, "");
  fprintf(ox->dfp, "      ");
  for (k = 0; k <= M;          k++) fprintf(ox->dfp, " %*d", width, k);
  for (k = 0; k < p7F_NXCELLS; k++) fprintf(ox->dfp, " %*s", width, p7_filtermx_DecodeX(k));
  fputc('\n', ox->dfp);

  fprintf(ox->dfp, "%*s", maxpfx, "");
  fprintf(ox->dfp, "      ");
  for (k = 0; k <= M+p7F_NXCELLS;  k++) fprintf(ox->dfp, " %*s", width, "--------");
  fputc('\n', ox->dfp);

  return eslOK;
}

/* Function:  p7_filtermx_DumpFBRow()
 * Synopsis:  Dump one row from fwd or bck version of the matrix.
 *
 * Purpose:   Dump current row <dpc> of forward or backward calculations from
 *            DP matrix <ox> for diagnostics. The index <rowi> is used
 *            as a row label, along with an additional free-text label
 *            <pfx>.  (The checkpointed backward implementation
 *            interleaves backward row calculations with recalculated
 *            fwd rows, both of which it is dumping; they need to be
 *            labeled something like "fwd" and "bck" to distinguish
 *            them in the debugging dump.)
 */
int
p7_filtermx_DumpFBRow(P7_FILTERMX *ox, int rowi, __m128 *dpc, char *pfx)
{
  union { __m128 v; float x[p7_VNF]; } u;
  float *v         = NULL;		/*  */
  int    Q         = ox->Qf;
  int    M         = ox->M;
  float *xc        = (float *) (dpc + Q*p7F_NSCELLS);
  int    logify    = (ox->dump_flags & p7_SHOW_LOG) ? TRUE : FALSE;
  int    maxpfx    = ox->dump_maxpfx;
  int    width     = ox->dump_width;
  int    precision = ox->dump_precision;
  int    k,q,z;
  int    status;

  ESL_ALLOC(v, sizeof(float) * ( (Q*p7_VNF) + 1));
  v[0] = 0.;

  /* Line 1. M cells: unpack, unstripe, print */
  for (q = 0; q < Q; q++) {
    u.v = P7F_MQ(dpc, q);
    for (z = 0; z < p7_VNF; z++) v[q+Q*z+1] = u.x[z];
  }
  fprintf(ox->dfp, "%*s %3d M", maxpfx, pfx, rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, " %*.*f", width, precision, (logify ? esl_logf(v[k]) : v[k]));

  /* Line 1 end: Specials */
  for (z = 0; z < p7F_NXCELLS; z++)
    fprintf(ox->dfp, " %*.*f", width, precision, (logify ? esl_logf(xc[z]) : xc[z]));
  fputc('\n', ox->dfp);

  /* Line 2: I cells: unpack, unstripe, print */
  for (q = 0; q < Q; q++) {
    u.v = P7F_IQ(dpc, q);
    for (z = 0; z < p7_VNF; z++) v[q+Q*z+1] = u.x[z];
  }
  fprintf(ox->dfp, "%*s %3d I", maxpfx, pfx, rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, " %*.*f", width, precision, (logify ? esl_logf(v[k]) : v[k]));
  fputc('\n', ox->dfp);

  /* Line 3. D cells: unpack, unstripe, print */
  for (q = 0; q < Q; q++) {
    u.v = P7F_DQ(dpc, q);
    for (z = 0; z < p7_VNF; z++) v[q+Q*z+1] = u.x[z];
  }
  fprintf(ox->dfp, "%*s %3d D", maxpfx, pfx, rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, " %*.*f", width, precision, (logify ? esl_logf(v[k]) : v[k]));
  fputc('\n', ox->dfp);
  fputc('\n', ox->dfp);

  free(v);
  return eslOK;

 ERROR:
  if (v) free(v);
  return status;
}


/* Function:  p7_filtermx_DumpMFRow()
 * Synopsis:  Dump one row from MSV version of a DP matrix.
 *
 * Purpose:   Dump current row of MSV calculations from DP matrix <ox>
 *            for diagnostics, and include the values of specials
 *            <xE>, etc. The index <rowi> for the current row is used
 *            as a row label. This routine has to be specialized for
 *            the layout of the MSVFilter() row, because it's all
 *            match scores dp[0..q..Q-1], rather than triplets of
 *            M,D,I.
 * 
 *            If <rowi> is 0, print a header first too.
 * 
 *            The output format is coordinated with <p7_refmx_Dump()> to
 *            facilitate comparison to a known answer.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. 
 */
int
p7_filtermx_DumpMFRow(P7_FILTERMX *ox, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC)
{
  __m128i *dp = (__m128i *) ox->dpf[0];	/* we may use up to allocW*validR bytes on our one MSV row */
  int      Q  = P7F_NQB(ox->M);	        /* and we actually use Q * sizeof(__m128i)                 */
  int      M  = ox->M;
  uint8_t *v  = NULL;		/* array of unstriped scores  */
  int      q,z,k;
  union { __m128i v; uint8_t i[16]; } tmp;
  int      status;

  ESL_ALLOC(v, sizeof(unsigned char) * ((Q*16)+1));
  v[0] = 0;

  /* Header (if we're on the 0th row)  */
  if (rowi == 0)
    {
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M;  k++) fprintf(ox->dfp, "%3d ", k);
      fprintf(ox->dfp, "%3s %3s %3s %3s %3s\n", "E", "N", "J", "B", "C");
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M+5;  k++) fprintf(ox->dfp, "%3s ", "---");
      fprintf(ox->dfp, "\n");
    }

  /* Unpack and unstripe, then print M's. */
  for (q = 0; q < Q; q++) {
    tmp.v = dp[q];
    for (z = 0; z < 16; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d M ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", v[k]);

  /* The specials */
  fprintf(ox->dfp, "%3d %3d %3d %3d %3d\n", xE, xN, xJ, xB, xC);

  /* I's are all 0's; print just to facilitate comparison. */
  fprintf(ox->dfp, "%4d I ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", 0);
  fprintf(ox->dfp, "\n");

  /* D's are all 0's too */
  fprintf(ox->dfp, "%4d D ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", 0);
  fprintf(ox->dfp, "\n\n");

  free(v);
  return eslOK;

ERROR:
  free(v);
  return status;
}



/* Function:  p7_filtermx_DumpVFRow()
 * Synopsis:  Dump current row of ViterbiFilter (int16) part of <ox> matrix.
 *
 * Purpose:   Dump current row of ViterbiFilter (int16) part of DP
 *            matrix <ox> for diagnostics, and include the values of
 *            specials <xE>, etc. The index <rowi> for the current row
 *            is used as a row label.
 *
 *            If <rowi> is 0, print a header first too.
 * 
 *            The output format is coordinated with <p7_refmx_Dump()> to
 *            facilitate comparison to a known answer.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_filtermx_DumpVFRow(P7_FILTERMX *ox, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC)
{
  __m128i *dp = (__m128i *) ox->dpf[0]; /* we may use up to allocW*validR bytes on our one MSV row */
  int      Q  = P7F_NQW(ox->M);		/* and we actually use Q * sizeof(__m128i)                 */
  int      M  = ox->M;
  int16_t *v  = NULL;		/* array of unstriped, uninterleaved scores  */
  int      q,z,k;
  union { __m128i v; int16_t i[8]; } tmp;
  int      status;

  ESL_ALLOC(v, sizeof(int16_t) * ((Q*8)+1));
  v[0] = 0;

  /* Header (if we're on the 0th row)
   */
  if (rowi == 0)
    {
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M;  k++) fprintf(ox->dfp, "%6d ", k);
      fprintf(ox->dfp, "%6s %6s %6s %6s %6s\n", "E", "N", "J", "B", "C");
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M+5;  k++) fprintf(ox->dfp, "%6s ", "------");
      fprintf(ox->dfp, "\n");
    }

  /* Unpack and unstripe, then print M's. */
  for (q = 0; q < Q; q++) {
    tmp.v = P7F_MQ(dp, q);
    for (z = 0; z < 8; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d M ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%6d ", v[k]);

  /* The specials */
  fprintf(ox->dfp, "%6d %6d %6d %6d %6d\n", xE, xN, xJ, xB, xC);

  /* Unpack and unstripe, then print I's. */
  for (q = 0; q < Q; q++) {
    tmp.v = P7F_IQ(dp, q);
    for (z = 0; z < 8; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d I ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%6d ", v[k]);
  fprintf(ox->dfp, "\n");

  /* Unpack, unstripe, then print D's. */
  for (q = 0; q < Q; q++) {
    tmp.v = P7F_DQ(dp, q);
    for (z = 0; z < 8; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d D ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%6d ", v[k]);
  fprintf(ox->dfp, "\n\n");

  free(v);
  return eslOK;

ERROR:
  free(v);
  return status;

}


#endif /*p7_DEBUGGING*/
/*---------------- end, debugging -------------------------------*/




/*****************************************************************
 * 3. Internal routines
 *****************************************************************/

/* set_row_layout()
 *
 * Determines the size of the a,b,c regions ("all", "between",
 * "checkpointed") of rows in the DP matrix. 
 *
 * Caller must have already set <ox->allocW> and <ox->R0>; they are
 * needed here.
 * 
 * Upon return, we've set the R{0abc} and L{abc} fields in the <ox>
 * structure.
 * 
 * <maxR> is the maximum number of rows the caller *wants* to use,
 * either because a <ramlimit>'ed allocation fits that number of rows,
 * or because an existing matrix has that number of valid rows.  We
 * will exceed this for one comparison if absolutely necessary, but
 * the next <_Reuse()> call will bring the allocation back down.
 * 
 * So there's three possibilities:
 *  1. A full matrix fits into our recommended max memory use.
 *  2. A checkpointed matrix fits into our recommended memory.
 *     Make as much of the matrix uncheckpointed as we can,
 *     using every row in maxR.
 *  3. We can't satisfy the recommended max memory, even fully
 *     checkpointed. Make a fully checkpointed matrix, in which
 *     R0+Ra+Rb+Rc will exceed maxR, and caller will have to 
 *     allocate ("redlined").
 */
static void
set_row_layout(P7_FILTERMX *ox, int allocL, int maxR)
{
  double Rbc      = minimum_rows(allocL);               
  int    minR_chk = ox->R0 + (int) ceil(Rbc);	          /* min # rows we need for checkpointing          */
  int    minR_all = ox->R0 + allocL;		          /* min # rows we need for full matrix            */

  if      (minR_all <= maxR) set_full        (ox, allocL);
  else if (minR_chk <= maxR) set_checkpointed(ox, allocL, maxR);
  else                       set_redlined    (ox, allocL, Rbc);
}

/* A "full" matrix is easy: Ra = La = L, using Ra+R0 <= maxR rows total. */
static void
set_full(P7_FILTERMX *ox, int L)
{
  ox->Ra     = L;
  ox->Rb     = 0;
  ox->Rc     = 0;
  ox->La     = L;
  ox->Lb     = 0;
  ox->Lc     = 0;
}

/* If we can fit a checkpointed matrix into maxR rows, 
 * then the trick is to use all maxR rows, making the
 * "a" region (all rows kept) as large as possible, to
 * minimize computation. This means solving a little
 * quadratic equation for Rb+Rc given L and maxR: see
 * <checkpointed_rows()> for that solution.
 */
static void
set_checkpointed(P7_FILTERMX *ox, int L, int R)
{
  double Rbc = checkpointed_rows(L, R - ox->R0);
  double Rc  = floor(Rbc);

  ox->Rc     = (int) Rc;
  ox->Rb     = (Rbc > Rc ? 1 : 0);
  ox->Ra     = R - ox->Rb - ox->Rc - ox->R0;
  ox->Lc     = ((ox->Rc + 2) * (ox->Rc + 1)) / 2 - 1;
  ox->La     = ox->Ra;
  ox->Lb     = L - ox->La - ox->Lc;
}   

/* If we can't fit in maxR rows, then we checkpoint
 * the entire matrix; R0+Ra+Rb+Rc > maxR.
 */
static void
set_redlined(P7_FILTERMX *ox, int L, double minR)
{
  double Rc = floor(minR);

  ox->Rc     = (int) Rc;
  ox->Rb     = (minR > Rc ? 1 : 0);
  ox->Ra     = 0;
  ox->Lc     = ((ox->Rc + 2) * (ox->Rc + 1)) / 2 - 1;
  ox->La     = 0;
  ox->Lb     = L - ox->La - ox->Lc;
}

/* minimum_rows()
 * 
 * Calculate the minimum number of rows needed to checkpoint a
 * forward matrix for a sequence of length L, exclusive of 
 * other constant row overhead (R0: two backwards rows, fwd[0]
 * initial row).
 * 
 * This value is a double; if it has a fraction, a partial checkpoint
 * block ("b" region) is needed, as in this typical code:
 *    Rbc  = minimum_rows(L);
 *    Rc   = floor(Rbc);
 *    Rb   = (Rbc > Rc ? 1 : 0);
 *    minR = (int) ceil(Rbc);    // or, Rc+Rb
 */
static double 
minimum_rows(int L)
{
  return (sqrt(9. + 8. * (double) L) - 3.) / 2.;
}

/* checkpointed_rows()
 * 
 * Given L and maxR, solve for the number of checkpointed
 * rows (Rb+Rc) we need. The value is a double; if it has
 * a fractional part, then a partial checkpoint block is
 * needed, Rb==1.
 * 
 * This equation is obtained by solving 
 *     L = Ra + (Rbc +2)(Rbc+1)/2 - 1
 * for Rbc (i.e. Rb+Rc) using the quadratic equation,
 * after substitution Ra = (maxR - Rbc - R0) to get
 * Rbc in terms of L,maxR.
 */
static double
checkpointed_rows(int L, int R)
{
  return (sqrt(1. + 8. * (double) (L - R)) - 1.) / 2.;
}
/*----------------- end, internals ------------------------------*/


/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
