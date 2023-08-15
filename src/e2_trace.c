/* E2_TRACE, the traceback structure.
 *
 * Contents:
 *   1. The E2_TRACE structure
 *   2. Access routines
 *   3. Debugging tools
 *   4. Creating traces by DP traceback
 *   5. Creating faux traces from existing MSAs
 *   6. Counting traces into new HMMs
 *   7. Unit tests
 *   8. Test driver
 *   9. Copyright and license information
 * 
 * Stylistic note: elements in a trace path are usually indexed by z.
 */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "e2.h"
#include "e1_model.h"
#include "e2_trace.h"

/*****************************************************************
 * 1. The E2_TRACE structure.
 *****************************************************************/

static E2_TRACE *trace_create_engine(int initial_nalloc, int initial_ndomalloc, int with_posteriors);

/* Function:  e2_trace_Create()
 * Synopsis:  Allocates a (growable, reusable) traceback.
 *
 * Purpose:   Allocates a traceback. 
 *  
 *            Tracebacks are growable. A reasonable initial internal
 *            allocation is made here, and routines that generate
 *            tracebacks will dynamically grow the trace as needed.
 *            
 *            Tracebacks are reusable. Usually a routine only
 *            allocates one, and reuses its memory over and over as
 *            new target sequences are aligned.
 *
 * Returns:   a pointer to the new <E2_TRACE> structure on success.
 *
 * Throws:    <NULL> on allocation error.
 */
E2_TRACE *
e2_trace_Create(void)
{
  int       initial_nalloc    = 256;
  int       initial_ndomalloc = 16;
  int       with_posteriors   = FALSE;
  return trace_create_engine(initial_nalloc, initial_ndomalloc, with_posteriors);
}

/* Function:  e2_trace_CreateWithPP()
 * Synopsis:  Allocates a traceback that includes posterior probs.
 * Incept:    ER, Thu Dec 15 10:45:19 EST 2011 [Janelia]
 *
 * Purpose:   Allocates a traceback that includes <tr->pp[z]> fields
 *            for posterior probabilities of emitted residues; 
 *            otherwise identical to <e2_trace_Create()>.
 */
E2_TRACE *
e2_trace_CreateWithPP(void)
{
  int       initial_nalloc    = 256;
  int       initial_ndomalloc = 16;
  int       with_posteriors   = TRUE;
  return trace_create_engine(initial_nalloc, initial_ndomalloc, with_posteriors);
}

static E2_TRACE *
trace_create_engine(int initial_nalloc, int initial_ndomalloc, int with_posteriors)
{
  E2_TRACE *tr      = NULL;
  int       status;

  ESL_ALLOC(tr, sizeof(E2_TRACE));
  tr->st = NULL;
  tr->k  = NULL;
  tr->a  = NULL;
  tr->i  = NULL;
  tr->j  = NULL;
  tr->pp = NULL;
  tr->M    = 0;
  tr->Lrow = 0;
  tr->Lcol = 0;
  tr->tfrom  = tr->tto  = NULL;
  tr->sqfrom = tr->sqto = NULL;
  tr->anfrom = tr->anto = NULL;
  /* The trace data itself */
  ESL_ALLOC(tr->st, sizeof(int ) * initial_nalloc);
  ESL_ALLOC(tr->k,  sizeof(int)  * initial_nalloc);
  ESL_ALLOC(tr->a,  sizeof(int)  * initial_nalloc);
  ESL_ALLOC(tr->i,  sizeof(int)  * initial_nalloc);
  ESL_ALLOC(tr->j,  sizeof(int)  * initial_nalloc);
  if (with_posteriors)
    ESL_ALLOC(tr->pp, sizeof(float) * initial_nalloc);
  tr->N      = 0;
  tr->nalloc = initial_nalloc;

  /* The trace's index: table of domain start/stop coords */
  ESL_ALLOC(tr->tfrom,   sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->tto,     sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->sqfrom,  sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->sqto,    sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->anfrom, sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->anto,   sizeof(int) * initial_ndomalloc);
  tr->ndom      = 0;
  tr->ndomalloc = initial_ndomalloc;
  return tr;

 ERROR:
  if (tr != NULL) e2_trace_Destroy(tr);
  return NULL;
}

/* Function:  e2_trace_Reuse()
 * Synopsis:  Prepare a trace for reuse.
 * Incept:    SRE, Tue Jan  9 13:02:34 2007 [Janelia]
 *
 * Purpose:   Reinitializes an existing trace object, reusing its
 *            memory.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      STL11/124
 */
int
e2_trace_Reuse(E2_TRACE *tr)
{
  tr->N    = 0;
  tr->M    = 0;
  tr->Lrow = 0;
  tr->Lcol = 0;
  tr->ndom = 0;
  return eslOK;
}

/* Function:  e2_trace_Grow()
 * Synopsis:  Grow the allocation for trace data.
 *
 * Purpose:   If <tr> can't fit another state, double its allocation for
 *            traceback data.
 *            
 *            This doesn't reallocate the domain index; see
 *            <e2_trace_GrowIndex()> or <e2_trace_GrowIndexTo()> for
 *            that.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; in this case, the data in
 *            <tr> are unaffected.
 */
int
e2_trace_Grow(E2_TRACE *tr)
{
  void *tmp;
  int   status;
  
  if (tr->N < tr->nalloc) return eslOK;

  ESL_RALLOC(tr->st, tmp, sizeof(int)  *2*tr->nalloc);
  ESL_RALLOC(tr->k,  tmp, sizeof(int)  *2*tr->nalloc);
  ESL_RALLOC(tr->a,  tmp, sizeof(int)  *2*tr->nalloc);
  ESL_RALLOC(tr->i,  tmp, sizeof(int)  *2*tr->nalloc);
  ESL_RALLOC(tr->j,  tmp, sizeof(int)  *2*tr->nalloc);
  if (tr->pp != NULL) ESL_RALLOC(tr->pp,  tmp, sizeof(float) *2*tr->nalloc);
  tr->nalloc *= 2;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  e2_trace_GrowIndex()
 * Synopsis:  Grows the allocation of the trace's domain index.
 * Incept:    SRE, Fri Jan  4 10:40:02 2008 [Janelia]
 *
 * Purpose:   If <tr> can't fit another domain in its index,
 *            double the allocation of the index in <tr>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; in this case, the 
 *            data in <tr> are unaffected.
 */
int
e2_trace_GrowIndex(E2_TRACE *tr)
{
  void *p;
  int   status;

  if (tr->ndom < tr->ndomalloc) return eslOK;

  ESL_RALLOC(tr->tfrom,   p, sizeof(int)*2*tr->ndomalloc);
  ESL_RALLOC(tr->tto,     p, sizeof(int)*2*tr->ndomalloc);
  ESL_RALLOC(tr->sqfrom,  p, sizeof(int)*2*tr->ndomalloc);
  ESL_RALLOC(tr->sqto,    p, sizeof(int)*2*tr->ndomalloc);
  ESL_RALLOC(tr->anfrom,  p, sizeof(int)*2*tr->ndomalloc);
  ESL_RALLOC(tr->anto,    p, sizeof(int)*2*tr->ndomalloc);
  tr->ndomalloc *= 2;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  e2_trace_GrowTo()
 * Synopsis:  Reallocates trace to a given minimum size.
 *
 * Purpose:   Reallocates a trace structure <tr> to hold a trace
 *            of at least length <N> states.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; the data in <tr>
 *            are unaffected by failure.
 */
int
e2_trace_GrowTo(E2_TRACE *tr, int N)
{
  int status;
  void *tmp;

  if (N < tr->nalloc) return eslOK; /* no-op */
  
  ESL_RALLOC(tr->st, tmp, sizeof(int)  * N);
  ESL_RALLOC(tr->k,  tmp, sizeof(int)  * N);
  ESL_RALLOC(tr->a,  tmp, sizeof(int)  * N);
  ESL_RALLOC(tr->i,  tmp, sizeof(int)  * N);
  ESL_RALLOC(tr->j,  tmp, sizeof(int)  * N);
  if (tr->pp != NULL) ESL_RALLOC(tr->pp,  tmp, sizeof(float) * N);
  tr->nalloc = N;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  e2_trace_GrowIndexTo()
 * Synopsis:  Reallocates domain index for a given minimum number.
 * Incept:    SRE, Fri Jan  4 10:47:43 2008 [Janelia]
 *
 * Purpose:   Reallocates the domain index in <tr> to index
 *            at least <ndom> domains.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure, in which case
 *            the data in <tr> are unaffected.
 */
int
e2_trace_GrowIndexTo(E2_TRACE *tr, int ndom)
{
  void *p;
  int   status;

  if (ndom < tr->ndomalloc) return eslOK;

  ESL_RALLOC(tr->tfrom,   p, sizeof(int)*ndom);
  ESL_RALLOC(tr->tto,     p, sizeof(int)*ndom);
  ESL_RALLOC(tr->sqfrom,  p, sizeof(int)*ndom);
  ESL_RALLOC(tr->sqto,    p, sizeof(int)*ndom);
  ESL_RALLOC(tr->anfrom,  p, sizeof(int)*ndom);
  ESL_RALLOC(tr->anto,    p, sizeof(int)*ndom);
  tr->ndomalloc = ndom;
  return eslOK;
  
 ERROR:
  return status;
}


/* Function:  e2_trace_Destroy()
 * Synopsis:  Frees a trace.
 *
 * Purpose:   Frees a trace structure <tr>.
 *
 * Returns:   (void)
 */
void 
e2_trace_Destroy(E2_TRACE *tr)
{
  if (tr == NULL) return;
  if (tr->st      != NULL) free(tr->st);
  if (tr->k       != NULL) free(tr->k);
  if (tr->a       != NULL) free(tr->a);
  if (tr->i       != NULL) free(tr->i);
  if (tr->j       != NULL) free(tr->j);
  if (tr->pp      != NULL) free(tr->pp);
  if (tr->tfrom   != NULL) free(tr->tfrom);
  if (tr->tto     != NULL) free(tr->tto);
  if (tr->sqfrom  != NULL) free(tr->sqfrom);
  if (tr->sqto    != NULL) free(tr->sqto);
  if (tr->anfrom  != NULL) free(tr->anfrom);
  if (tr->anto    != NULL) free(tr->anto);
  free(tr);
  return;
}

/* Function:  e2_trace_DestroyArray()
 *
 * Purpose:   Frees an array of <N> trace structures, <tr[0..N-1]>.
 *
 * Returns:   (void)
 */
void 
e2_trace_DestroyArray(E2_TRACE **tr, int N)
{
  int idx;

  if (tr == NULL) return;
  for (idx = 0; idx < N; idx++)
    {
      if (tr[idx] == NULL) continue;
      e2_trace_Destroy(tr[idx]);
    }
  free(tr);
  return;
}

/*---------------------- end, E2_TRACE --------------------------*/

/*****************************************************************
 * 3. Debugging tools.
 *****************************************************************/
/* Function:  e2_trace_Dump()
 * Incept:    SRE, Fri Jan  5 09:26:04 2007 [Janelia]
 *
 * Purpose:   Dumps internals of a traceback structure <tr> to <fp>.
 *            If <gm> is non-NULL, also prints transition/emission scores.
 *            If <dsq> is non-NULL, also prints residues (using alphabet
 *            in the <gm>).
 *            
 * Args:      fp   - stream to dump to (often stdout)
 *            tr   - trace to dump
 *            gm   - NULL, or score profile corresponding to trace
 *            dsq  - NULL, or digitized seq corresponding to trace        
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEINVAL> if trace contains something corrupt or invalid;
 *            in this case, dump will be aborted, possibly after partial
 *            output.
 */
int
e2_trace_Dump(FILE *fp, const E2_TRACE *tr, const E2_PROFILE *gm, const E2HMMER_PROFILE *gm7, const PSQ *sql, const PSQ *sqr, const PSQ *sqa) 
{
  float  accuracy = 0.0f;
  float  sc       = 0.0f;
  float  tsc;
  PSQ   *sqrow;
  PSQ   *sqcol;
  float *prow = NULL;
  float *pcol = NULL;
  float *pa   = NULL;
  int    K;
  int    xi, xj, xa;
  int    z;		/* counter for trace position */
  int    status;
 
  if (tr == NULL) { fprintf(fp, " [ trace is NULL ]\n"); return eslOK; }

  if (gm == NULL) 
    {		/* Yes, this does get used: during model construction. */ 
      fprintf(fp, "st     a       i     j   - traceback len %d\n", tr->N);
      fprintf(fp, "--   ----    ----- -----\n");
      for (z = 0; z < tr->N; z++) {
	fprintf(fp, "%1s  %d/%d %4d %6d %6d \n", e2_model_DecodeStatetype(tr->st[z]), z, tr->N, tr->a[z], tr->i[z], tr->j[z]);
      } 
    } 
  else if (gm != NULL)
    {
      if (sql != NULL && sqr != NULL && sqa != NULL) {
	K = sql->abc->K;
	if (K != sqr->abc->K || K != sqa->abc->K) { status = eslFAIL; goto ERROR; }
	ESL_ALLOC(prow, sizeof(float) * (K+1));
	ESL_ALLOC(pcol, sizeof(float) * (K+1));
 	ESL_ALLOC(pa,   sizeof(float) * (K+1));
     }
      else { status = eslFAIL; goto ERROR; }

      sqrow = (tr->rowsq == e2P_SL)? (PSQ *)sql : (PSQ *)sqr;
      sqcol = (tr->rowsq == e2P_SL)? (PSQ *)sqr : (PSQ *)sql;

      fprintf(fp, "st     a       i    j       tsc    esc        pp       a     i   j   - traceback len %d\n", tr->N);
      fprintf(fp, "--   ----    ----- -----  ------  ------    ------    ---   --- ---\n");
      for (z = 0; z < tr->N; z++) 
	{
	  if (z < tr->N-1) 
	    {
	      status = e2_profile_GetT(gm, tr->st[z], tr->st[z+1], &tsc);
	      if (status != eslOK) goto ERROR;
	    }
	  else tsc = 0.0f;
	  
	  fprintf(fp, "%1s  %4d %6d %6d  %8.4f", e2_model_DecodeStatetype(tr->st[z]), tr->a[z], tr->i[z], tr->j[z], tsc);
	  sc += tsc;
	  if (sql != NULL && sqr != NULL && sqa != NULL) {	    
	    
	    switch(tr->st[z]) {
	    case e2T_S:
	      break;
	    case e2T_BB:
		if (tr->pp != NULL) { fprintf(fp, "          %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      break;

	    case e2T_N1:
	      if (z > 0 && tr->M && tr->st[z-1] == tr->st[z]) {
		psq_ProfProbs(tr->i[z], sqrow, prow);
		xi = esl_vec_FArgMax(prow, K);
		fprintf(fp, " %8.4f", e2P_FLSC(gm, xi));
		sc += e2P_FLSC(gm, xi);
		if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
		fprintf(fp, "     %c     %c  %c", '-', tolower(gm->abc->sym[xi]), '-');
	      }
	      break;

	    case e2T_J1:
	      if (z > 0 && tr->st[z-1] == tr->st[z]) {
		psq_ProfProbs(tr->i[z], sqrow, prow);
		xi = esl_vec_FArgMax(prow, K);
		fprintf(fp, " %8.4f", e2P_FLSC(gm, xi));
		sc += e2P_FLSC(gm, xi);
		if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
		fprintf(fp, "     %c     %c  %c", '-', tolower(gm->abc->sym[xi]), '-');
	      }
	      break;

	    case e2T_C1:
	      if (z > 0 && tr->st[z-1] == tr->st[z]) {
		psq_ProfProbs(tr->i[z], sqrow, prow);
		xi = esl_vec_FArgMax(prow, K);
		fprintf(fp, " %8.4f", e2P_FLSC(gm, xi));
		sc += e2P_FLSC(gm, xi);
		if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
		fprintf(fp, "     %c     %c  %c", '-', tolower(gm->abc->sym[xi]), '-');
	      }
	      break;

	    case e2T_N2:
	      if (z > 0 && tr->st[z-1] == tr->st[z]) {
		psq_ProfProbs(tr->j[z], sqcol, pcol);
		xj = esl_vec_FArgMax(pcol, K);
		fprintf(fp, " %8.4f", e2P_FRSC(gm, xj));
		sc += e2P_FRSC(gm, xj);
		if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
		fprintf(fp, "     %c     %c  %c", '-', '-', tolower(gm->abc->sym[xj]));
	      }
	      break;
	      
	    case e2T_J2:
	      if (z > 0 && tr->st[z-1] == tr->st[z]) {
		psq_ProfProbs(tr->j[z], sqcol, pcol);
		xj = esl_vec_FArgMax(pcol, K);
		fprintf(fp, " %8.4f", e2P_FRSC(gm, xj));
		sc += e2P_FRSC(gm, xj);
		if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
		fprintf(fp, "     %c     %c  %c", '-', '-', tolower(gm->abc->sym[xj]));
	      }
	      break;
	      
	    case e2T_C2:
	      if (z > 0 && tr->st[z-1] == tr->st[z]) {
		psq_ProfProbs(tr->j[z], sqcol, pcol);
		xj = esl_vec_FArgMax(pcol, K);
		fprintf(fp, " %8.4f", e2P_FRSC(gm, xj));
		sc += e2P_FRSC(gm, xj);
		if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
		fprintf(fp, "     %c     %c  %c", '-', '-', tolower(gm->abc->sym[xj]));
	      }
	      break;

	    case e2T_IB:
	      psq_ProfProbs(tr->i[z], sqrow, prow);
	      xi = esl_vec_FArgMax(prow, K);
	      fprintf(fp, " %8.4f", e2P_ILSC(gm, xi));
	      sc += e2P_ILSC(gm, xi);
	      if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      fprintf(fp, "     %c     %c  %c", '-', gm->abc->sym[xi], '-');
	      break;
	      
	    case e2T_SS:
	      psq_ProfProbs(tr->i[z], sqrow, prow);
	      psq_ProfProbs(tr->j[z], sqcol, pcol);
	      psq_ProfProbs(tr->a[z], sqa, pa);
 	      xi = esl_vec_FArgMax(prow, K);
	      xj = esl_vec_FArgMax(pcol, K);
	      xa = esl_vec_FArgMax(pa, K);
	      fprintf(fp, " %8.4f", e2P_SSSC(gm, xi, xj));
	      sc += e2P_SSSC(gm, xi, xj);
	      if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      fprintf(fp, "     %c     %c  %c", gm->abc->sym[xa], gm->abc->sym[xi], gm->abc->sym[xj]);
	      break;

	    case e2T_DS:
	      psq_ProfProbs(tr->j[z], sqcol, pcol);
	      psq_ProfProbs(tr->a[z], sqa, pa);
	      xj = esl_vec_FArgMax(pcol, K);
	      xa = esl_vec_FArgMax(pa, K);
	      fprintf(fp, " %8.4f", e2P_SRSC(gm, xj));
	      sc += e2P_SRSC(gm, xj);
	      if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      fprintf(fp, "     %c     %c  %c", gm->abc->sym[xa], '-', gm->abc->sym[xj]);
	      break;

	    case e2T_IS:
	      psq_ProfProbs(tr->i[z], sqrow, prow);
	      xi = esl_vec_FArgMax(prow, K);
	      fprintf(fp, " %8.4f", e2P_ILSC(gm, xi));
	      sc += e2P_ILSC(gm, xi);
	      if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      fprintf(fp, "     %c     %c  %c", '-', gm->abc->sym[xi], '-');
	      break;
	    case e2T_SD:
	      psq_ProfProbs(tr->i[z], sqrow, prow);
	      psq_ProfProbs(tr->a[z], sqa, pa);
	      xi = esl_vec_FArgMax(prow, K);
	      xa = esl_vec_FArgMax(pa, K);
	      fprintf(fp, " %8.4f", e2P_SLSC(gm, xi));
	      sc += e2P_SLSC(gm, xi);
	      if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      fprintf(fp, "     %c     %c  %c", gm->abc->sym[xa], gm->abc->sym[xi], '-');
	      break;

	    case e2T_DD:
	      psq_ProfProbs(tr->a[z], sqa, pa);
	      xa = esl_vec_FArgMax(pa, K);
	      fprintf(fp, " %8s", "-");
	      if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      fprintf(fp, "     %c     %c  %c", gm->abc->sym[xa], '-', '-');
	      break;
	    case e2T_ID:
	      psq_ProfProbs(tr->i[z], sqrow, prow);
	      xi = esl_vec_FArgMax(prow, K);
	      fprintf(fp, " %8.4f", e2P_ILSC(gm, xi));
	      sc += e2P_ILSC(gm, xi);
	      if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      fprintf(fp, "     %c     %c  %c", '-', gm->abc->sym[xi], '-');
	      break;
	    case e2T_ii:
	      psq_ProfProbs(tr->i[z], sqrow, prow);
	      xi = esl_vec_FArgMax(prow, K);
	      fprintf(fp, " %8.4f", e2P_ILSC(gm, xi));
	      sc += e2P_ILSC(gm, xi);
	      if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      fprintf(fp, "     %c     %c  %c", '-', '-', gm->abc->sym[xi]);
	      break;
	    case e2T_BI:
	      psq_ProfProbs(tr->j[z], sqcol, pcol);
	      xj = esl_vec_FArgMax(pcol, K);
	      fprintf(fp, " %8.4f", e2P_IRSC(gm, xj));
	      sc += e2P_IRSC(gm, xj);
	      if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      fprintf(fp, "     %c     %c  %c", '-', '-', gm->abc->sym[xj]);
	      break;
	    case e2T_SI:
	      psq_ProfProbs(tr->j[z], sqcol, pcol);
	      xj = esl_vec_FArgMax(pcol, K);
	      fprintf(fp, " %8.4f", e2P_IRSC(gm, xj));
	      sc += e2P_IRSC(gm, xj);
	      if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      fprintf(fp, "     %c     %c  %c", '-', '-', gm->abc->sym[xj]);
	      break;
	    case e2T_DI:
	      psq_ProfProbs(tr->j[z], sqcol, pcol);
	      xj = esl_vec_FArgMax(pcol, K);
	      fprintf(fp, " %8.4f", e2P_IRSC(gm, xj));
	      sc += e2P_IRSC(gm, xj);
	      if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      fprintf(fp, "     %c     %c  %c", '-', '-', gm->abc->sym[xj]);
	      break;
	    case e2T_II:
	      psq_ProfProbs(tr->j[z], sqcol, pcol);
	      xj = esl_vec_FArgMax(pcol, K);
	      fprintf(fp, " %8.4f", e2P_IRSC(gm, xj));
	      sc += e2P_IRSC(gm, xj);
	      if (tr->pp != NULL) { fprintf(fp, " %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      fprintf(fp, "     %c     %c  %c", '-', '-', gm->abc->sym[xj]);
	      break;
	    case e2T_EE:
		if (tr->pp != NULL) { fprintf(fp, "          %8.4f", tr->pp[z]); accuracy += tr->pp[z]; }
	      break;
	    case e2T_T:
	      break;
	    default: status = eslFAIL; goto ERROR;
	    }
	    
	  }
	  else  fprintf(fp, " %8s %8s %c %c %c", "-", "-", '-', '-', '-');
	  fputs("\n", fp);
	}
      fprintf(fp, "                -------- ----------------- --------\n");
      fprintf(fp, "                  total:     %8.4f     %8.4f\n\n", sc, accuracy);

      if (prow) free(prow);
      if (pcol) free(pcol);
      if (pa) free(pa);         
   }
  else if (gm7 != NULL)
    {
    }

  
  return eslOK;

 ERROR:
  if (prow) free(prow);
  if (pcol) free(pcol);
  if (pa) free(pa);
  return status;
}


/*------------------ end, debugging tools -----------------------*/




/*****************************************************************
 * 4. Creating traces by DP traceback
 *****************************************************************/

/* Function:  e2_trace_Append()
 * Synopsis:  Add an element (state/residue) to a growing trace.
 *
 * Purpose:   Adds an element to a trace <tr> that is growing
 *            left-to-right. The element is defined by a state type
 *            <st> (such as <E2T_M>); a node index <k> (1..M for
 *            M,D,I main states; else 0); and a dsq position <i> (1..L
 *            for emitters, else 0).
 *            
 *            For CNJ states, which emit on transition, by convention
 *            we associate the emission with the downstream state; therefore
 *            the first state in any run of CNJ states has i=0. 
 *            
 *            Reallocates the trace (by doubling) if necessary.
 *            
 *            Caller can grow a trace right-to-left too, if it
 *            plans to call <E2_trace_Reverse()>. 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure. The element is successfully
 *            added, but no more elements can be added before this trace is
 *            destroyed.
 *            
 *            <eslEINVAL> if you try to add an element to a trace whose
 *            reallocation has already failed.
 */
int
e2_trace_Append(E2_TRACE *tr, int st, int k, int a, int i, int j)
{
  int N = tr->N;
  int status;

  if ((status = e2_trace_Grow(tr)) != eslOK) return status;

  switch (st) {
    /* Nonemitting states */
  case e2T_S:  tr->st[N] = e2T_S;  tr->i[N] = 0; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = k; break;
  case e2T_BB: tr->st[N] = e2T_BB; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; break;
  case e2T_EE: tr->st[N] = e2T_EE; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; break;
  case e2T_T:  tr->st[N] = e2T_T;  tr->i[N] = 0; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = k; break;

  case e2T_XX: tr->st[N] = e2T_XX; tr->i[N] = 0; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = k; break;

   /* Substitution/Deletion/Ancestral states: */
  case e2T_SS: tr->st[N] = e2T_SS; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; break;
  case e2T_DS: tr->st[N] = e2T_DS; tr->i[N] = 0; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; break;
  case e2T_SD: tr->st[N] = e2T_SD; tr->i[N] = i; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = k; break;
  case e2T_DD: tr->st[N] = e2T_DD; tr->i[N] = 0; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = k; break;
  
    /* Insertion states, but k valid */
  case e2T_IB: tr->st[N] = e2T_IB; tr->i[N] = i; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = k; break;
  case e2T_IS: tr->st[N] = e2T_IS; tr->i[N] = i; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = k; break;
  case e2T_ID: tr->st[N] = e2T_ID; tr->i[N] = i; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = k; break;
  case e2T_ii: tr->st[N] = e2T_ii; tr->i[N] = i; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = k; break;

  case e2T_BI: tr->st[N] = e2T_BI; tr->i[N] = 0; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; break;
  case e2T_SI: tr->st[N] = e2T_SI; tr->i[N] = 0; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; break;
  case e2T_DI: tr->st[N] = e2T_DI; tr->i[N] = 0; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; break;
  case e2T_II: tr->st[N] = e2T_II; tr->i[N] = 0; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; break;

    /* extra states, not valide k */
  case e2T_N1: tr->st[N] = e2T_N1; tr->i[N] = i; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = 0; break;
  case e2T_J1: tr->st[N] = e2T_J1; tr->i[N] = i; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = 0; break;
  case e2T_C1: tr->st[N] = e2T_C1; tr->i[N] = i; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = 0; break;
  case e2T_N2: tr->st[N] = e2T_N2; tr->i[N] = 0; tr->j[N] = j; tr->a[N] = a; tr->k[N] = 0; break;
  case e2T_J2: tr->st[N] = e2T_J2; tr->i[N] = 0; tr->j[N] = j; tr->a[N] = a; tr->k[N] = 0; break;
  case e2T_C2: tr->st[N] = e2T_C2; tr->i[N] = 0; tr->j[N] = j; tr->a[N] = a; tr->k[N] = 0; break;

  default:    ESL_EXCEPTION(eslEINVAL, "no such state; can't append");
  }

  tr->N ++;

  return eslOK;
}

/* Function:  e2_trace_AppendWithPP()
 * Synopsis:  Add element to growing trace, with posterior probability.
 *
 * Purpose:   Same as <e2_trace_Append()>, but also records a posterior
 *            probability estimate for emitted residues. <pp> is assumed to be
 *            zero for nonemitting states even if a nonzero argument is
 *            mistakenly passed. 
 */
int
e2_trace_AppendWithPP(E2_TRACE *tr, int st, int k, int a, int i, int j, float pp)
{
  int N  = tr->N;
  int status;

  if ((status = e2_trace_Grow(tr)) != eslOK) return status;

  switch (st) {
    /* Nonemitting states, but a valid: */
  case e2T_S:  tr->st[N] = e2T_S;  tr->i[N] = 0; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_BB: tr->st[N] = e2T_BB; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_EE: tr->st[N] = e2T_EE; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_T:  tr->st[N] = e2T_T;  tr->i[N] = 0; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_XX: tr->st[N] = e2T_XX; tr->i[N] = 0; tr->j[N] = 0; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;

  /* Substitution/Deletion/Ancestral states: */
  case e2T_SS: tr->st[N] = e2T_SS; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_DS: tr->st[N] = e2T_DS; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_SD: tr->st[N] = e2T_SD; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_DD: tr->st[N] = e2T_DD; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  
    /* Insertion states, but a valid k */
  case e2T_IB: tr->st[N] = e2T_IB; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_IS: tr->st[N] = e2T_IS; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_ID: tr->st[N] = e2T_ID; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_ii: tr->st[N] = e2T_ii; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;

  case e2T_BI: tr->st[N] = e2T_BI; tr->i[N] = 0; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_SI: tr->st[N] = e2T_SI; tr->i[N] = 0; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_DI: tr->st[N] = e2T_DI; tr->i[N] = 0; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
  case e2T_II: tr->st[N] = e2T_II; tr->i[N] = 0; tr->j[N] = j; tr->a[N] = a; tr->k[N] = k; tr->pp[N] = pp; break;
    
    /* extra states, not valide k */
  case e2T_N1: tr->st[N] = e2T_N1; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = 0; tr->pp[N] = pp; break;
  case e2T_J1: tr->st[N] = e2T_J1; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = 0; tr->pp[N] = pp; break;
  case e2T_C1: tr->st[N] = e2T_C1; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = 0; tr->pp[N] = pp; break;
  case e2T_N2: tr->st[N] = e2T_N2; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = 0; tr->pp[N] = pp; break;
  case e2T_J2: tr->st[N] = e2T_J2; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = 0; tr->pp[N] = pp; break;
  case e2T_C2: tr->st[N] = e2T_C2; tr->i[N] = i; tr->j[N] = j; tr->a[N] = a; tr->k[N] = 0; tr->pp[N] = pp; break;

  default:    ESL_EXCEPTION(eslEINVAL, "no such state (%d); can't append", st);
  }

  tr->N ++;

   return eslOK;
}

/* Function: e2_trace_Reverse()
 * Synopsis: Reverse the arrays in a traceback structure.
 * 
 * Purpose:  Reverse the arrays in a traceback structure.  Tracebacks
 *           from DP algorithms are collected backwards, and they call this
 *           function when they're done.
 *           
 * Args:     tr - the traceback to reverse. tr->N must be set.
 *                
 * Return:   <eslOK> on success; <tr> is modified.
 */                
int
e2_trace_Reverse(E2_TRACE *tr)
{
  int    z;
  int    tmp;
  float  tmpf;
 
  /* Reverse the trace in place. */ 
  for (z = 0; z < tr->N/2; z++)
    {
      tmp = tr->st[tr->N-z-1];        tr->st[tr->N-z-1] = tr->st[z];        tr->st[z] = tmp;
      tmp = tr->i[tr->N-z-1];         tr->i[tr->N-z-1]  = tr->i[z];         tr->i[z]  = tmp;
      tmp = tr->j[tr->N-z-1];         tr->j[tr->N-z-1]  = tr->j[z];         tr->j[z]  = tmp;
      tmp = tr->M - tr->a[tr->N-z-1]; tr->a[tr->N-z-1]  = tr->M - tr->a[z]; tr->a[z]  = tmp; /* need to reverse the order of a here too */
      tmp = tr->k[tr->N-z-1];         tr->k[tr->N-z-1] = tr->k[z];          tr->k[z]  = tmp; /* need to reverse the order of k here too */
      if (tr->pp != NULL) {
	tmpf = tr->pp[tr->N-z-1];     tr->pp[tr->N-z-1] = tr->pp[z];        tr->pp[z] = tmpf;
      }  
    }
  
  /* We need to worry about the middle residue in odd-length N, because 
     while we're in-place, the value of k needs to be reversed. 
   */
  if (tr->N%2 != 0) {
    tmp = tr->M - tr->a[tr->N-z-1]; tr->a[tr->N-z-1]  = tr->M - tr->a[z]; tr->a[z]  = tmp;
    tmp = tr->k[tr->N-z-1];         tr->k[tr->N-z-1]  = tr->k[z];         tr->k[z]  = tmp;
   }

   return eslOK;
}

/* Function:  e2_trace_Index()
 * Synopsis:  Internally index the domains in a trace.
 * Incept:    SRE, Fri Jan  4 11:12:24 2008 [Janelia]
 *
 * Purpose:   Create an internal index of the domains in <tr>.
 *            This makes calls to <GetDomainCount()> and 
 *            <GetDomainCoords()> more efficient, and it is
 *            a necessary prerequisite for creating alignments
 *            of any individual domains in a multidomain trace with
 *            <e2_alidisplay_Create()>.
 *
 * Returns:   <eslOK> on success. 
 *
 * Throws:    <eslEMEM> on allocation failure, in which case the
 *            data in the trace is still fine, but the domain index
 *            table isn't constructed.
 */
int
e2_trace_Index(E2_TRACE *tr)
{
  int curst, prvst; /* current and previous states */
  int z;
  int status;

  tr->ndom = 0;
  prvst = e1T_B;

  for (z = 0; z < tr->N; z++)
    {
      curst = tr->st[z];
      switch (curst) {
      case e1T_B:
	if (z != 0) return eslFAIL;
	break;

      case e1T_S:
	if (prvst == e1T_I) { /* begin of a domain */
	  if ((status = e2_trace_GrowIndex(tr)) != eslOK) goto ERROR;
	  tr->tfrom[tr->ndom]  = z;
	  tr->sqfrom[tr->ndom] = tr->i[z];
	  tr->anfrom[tr->ndom] = tr->a[z];
	}
	else {
	  tr->sqto[tr->ndom]   = tr->i[z];
	  tr->anto[tr->ndom]   = tr->a[z];
	}
	break;

      case e1T_D:
	if (prvst == e1T_I) { /* begin of a domain */
	  if ((status = e2_trace_GrowIndex(tr)) != eslOK) goto ERROR;
	  tr->tfrom[tr->ndom]  = z;
	  tr->sqfrom[tr->ndom] = tr->i[z];
	  tr->anfrom[tr->ndom] = tr->a[z];
	}
	else {
	  tr->anto[tr->ndom]   = tr->a[z];
	}
	break;

      case e1T_I:
	if (prvst == e1T_S || prvst == e1T_D) { /* end of a domain */
	  tr->tto[tr->ndom] = z;
	  tr->ndom++;
	}
	else {
	  tr->sqto[tr->ndom] = tr->i[z];
	}
	break;

      case e1T_E:
	if (z != tr->N) return eslFAIL;
	break;
      }
      prvst = curst;
    }
  return eslOK;
  
 ERROR:
  return status;
}



/*----------- end, creating traces by DP traceback ---------------*/


