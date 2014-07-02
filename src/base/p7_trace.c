/* P7_TRACE: the traceback structure.
 *
 * Contents:
 *   1. The P7_TRACE structure
 *   2. Access routines
 *   3. Debugging tools
 *   4. Visualization tools
 *   5. Creating traces by DP traceback
 *   6. Creating faux traces from existing MSAs
 *   7. Counting traces into new HMMs
 *   8. Unit tests
 *   9. Test driver
 *  10. Example
 *  11. Copyright and license information
 * 
 * Stylistic note: elements in a trace path are usually indexed by z.
 */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"	/* used in _PlotHeatMap() */
#include "esl_msa.h"
#include "esl_random.h"

#include "base/p7_hmm.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

/*****************************************************************
 * 1. The P7_TRACE structure
 *****************************************************************/

static P7_TRACE *trace_create_engine(int initial_nalloc, int initial_ndomalloc, int with_posteriors);

/* Function:  p7_trace_Create()
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
 * Returns:   a pointer to the new <P7_TRACE> structure on success.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_TRACE *
p7_trace_Create(void)
{
  int       initial_nalloc    = 256;
  int       initial_ndomalloc = 16;
  int       with_posteriors   = FALSE;
  return trace_create_engine(initial_nalloc, initial_ndomalloc, with_posteriors);
}

/* Function:  p7_trace_CreateWithPP()
 * Synopsis:  Allocates a traceback that includes posterior probs.
 *
 * Purpose:   Allocates a traceback that includes <tr->pp[z]> fields
 *            for posterior probabilities of emitted residues; 
 *            otherwise identical to <p7_trace_Create()>.
 */
P7_TRACE *
p7_trace_CreateWithPP(void)
{
  int       initial_nalloc    = 256;
  int       initial_ndomalloc = 16;
  int       with_posteriors   = TRUE;
  return trace_create_engine(initial_nalloc, initial_ndomalloc, with_posteriors);
}

static P7_TRACE *
trace_create_engine(int initial_nalloc, int initial_ndomalloc, int with_posteriors)
{
  P7_TRACE *tr      = NULL;
  int       status;

  ESL_ALLOC(tr, sizeof(P7_TRACE));
  tr->st = NULL;
  tr->k  = NULL;
  tr->i  = NULL;
  tr->pp = NULL;
  tr->M  = 0;
  tr->L  = 0;
  tr->tfrom   = tr->tto   = NULL;
  tr->sqfrom  = tr->sqto  = NULL;
  tr->hmmfrom = tr->hmmto = NULL;
  tr->anch                = NULL;

  /* The trace data itself */
  ESL_ALLOC(tr->st,   sizeof(char) * initial_nalloc);
  ESL_ALLOC(tr->k,    sizeof(int)  * initial_nalloc);
  ESL_ALLOC(tr->i,    sizeof(int)  * initial_nalloc);
  if (with_posteriors)
    ESL_ALLOC(tr->pp, sizeof(float) * initial_nalloc);
  tr->N      = 0;
  tr->nalloc = initial_nalloc;

  /* The trace's index: table of domain start/stop coords */
  ESL_ALLOC(tr->tfrom,   sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->tto,     sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->sqfrom,  sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->sqto,    sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->hmmfrom, sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->hmmto,   sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->anch,    sizeof(int) * initial_ndomalloc);
  tr->ndom      = 0;
  tr->ndomalloc = initial_ndomalloc;
  return tr;

 ERROR:
  if (tr != NULL) p7_trace_Destroy(tr);
  return NULL;
}


/* Function:  p7_trace_Reuse()
 * Synopsis:  Prepare a trace for reuse.
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
p7_trace_Reuse(P7_TRACE *tr)
{
  tr->N    = 0;
  tr->M    = 0;
  tr->L    = 0;
  tr->ndom = 0;
  return eslOK;
}

/* Function:  p7_trace_Grow()
 * Synopsis:  Grow the allocation for trace data.
 *
 * Purpose:   If <tr> can't fit another state, double its allocation for
 *            traceback data.
 *            
 *            This doesn't reallocate the domain index; see
 *            <p7_trace_GrowIndex()> or <p7_trace_GrowIndexTo()> for
 *            that.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; in this case, the data in
 *            <tr> are unaffected.
 */
int
p7_trace_Grow(P7_TRACE *tr)
{
  void *tmp;
  int   status;
  
  if (tr->N < tr->nalloc) return eslOK;

  ESL_RALLOC(tr->st, tmp, sizeof(char) *2*tr->nalloc);
  ESL_RALLOC(tr->k,  tmp, sizeof(int)  *2*tr->nalloc);
  ESL_RALLOC(tr->i,  tmp, sizeof(int)  *2*tr->nalloc);
  if (tr->pp != NULL) ESL_RALLOC(tr->pp,  tmp, sizeof(float) *2*tr->nalloc);
  tr->nalloc *= 2;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_trace_GrowIndex()
 * Synopsis:  Grows the allocation of the trace's domain index.
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
p7_trace_GrowIndex(P7_TRACE *tr)
{
  int   status;

  if (tr->ndom < tr->ndomalloc) return eslOK;

  ESL_REALLOC(tr->tfrom,   sizeof(int)*2*tr->ndomalloc);
  ESL_REALLOC(tr->tto,     sizeof(int)*2*tr->ndomalloc);
  ESL_REALLOC(tr->sqfrom,  sizeof(int)*2*tr->ndomalloc);
  ESL_REALLOC(tr->sqto,    sizeof(int)*2*tr->ndomalloc);
  ESL_REALLOC(tr->hmmfrom, sizeof(int)*2*tr->ndomalloc);
  ESL_REALLOC(tr->hmmto,   sizeof(int)*2*tr->ndomalloc);
  ESL_REALLOC(tr->anch,    sizeof(int)*2*tr->ndomalloc);
  tr->ndomalloc *= 2;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_trace_GrowTo()
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
p7_trace_GrowTo(P7_TRACE *tr, int N)
{
  int status;
  void *tmp;

  if (N < tr->nalloc) return eslOK; /* no-op */
  
  ESL_RALLOC(tr->st, tmp, sizeof(char) *N);
  ESL_RALLOC(tr->k,  tmp, sizeof(int)  *N);
  ESL_RALLOC(tr->i,  tmp, sizeof(int)  *N);
  if (tr->pp != NULL) ESL_RALLOC(tr->pp,  tmp, sizeof(float) *N);
  tr->nalloc = N;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_trace_GrowIndexTo()
 * Synopsis:  Reallocates domain index for a given minimum number.
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
p7_trace_GrowIndexTo(P7_TRACE *tr, int ndom)
{
  int   status;

  if (ndom < tr->ndomalloc) return eslOK;

  ESL_REALLOC(tr->tfrom,   sizeof(int)*ndom);
  ESL_REALLOC(tr->tto,     sizeof(int)*ndom);
  ESL_REALLOC(tr->sqfrom,  sizeof(int)*ndom);
  ESL_REALLOC(tr->sqto,    sizeof(int)*ndom);
  ESL_REALLOC(tr->hmmfrom, sizeof(int)*ndom);
  ESL_REALLOC(tr->hmmto,   sizeof(int)*ndom);
  ESL_REALLOC(tr->hmmto,   sizeof(int)*ndom);
  ESL_REALLOC(tr->anch,    sizeof(int)*ndom);
  tr->ndomalloc = ndom;
  return eslOK;
  
 ERROR:
  return status;
}


/* Function:  p7_trace_Destroy()
 * Synopsis:  Frees a trace.
 *
 * Purpose:   Frees a trace structure <tr>.
 *
 * Returns:   (void)
 */
void 
p7_trace_Destroy(P7_TRACE *tr)
{
  if (tr == NULL) return;
  if (tr->st)      free(tr->st);
  if (tr->k)       free(tr->k);
  if (tr->i)       free(tr->i);
  if (tr->pp)      free(tr->pp);
  if (tr->tfrom)   free(tr->tfrom);
  if (tr->tto)     free(tr->tto);
  if (tr->sqfrom)  free(tr->sqfrom);
  if (tr->sqto)    free(tr->sqto);
  if (tr->hmmfrom) free(tr->hmmfrom);
  if (tr->hmmto)   free(tr->hmmto);
  if (tr->anch)    free(tr->anch);
  free(tr);
  return;
}

/* Function:  p7_trace_DestroyArray()
 *
 * Purpose:   Frees an array of <N> trace structures, <tr[0..N-1]>.
 *
 * Returns:   (void)
 */
void 
p7_trace_DestroyArray(P7_TRACE **tr, int N)
{
  int idx;

  if (tr == NULL) return;
  for (idx = 0; idx < N; idx++)
    {
      if (tr[idx] == NULL) continue;
      p7_trace_Destroy(tr[idx]);
    }
  free(tr);
  return;
}
/*---------------------- end, P7_TRACE --------------------------*/




/*****************************************************************
 * 2. Access routines
 *****************************************************************/

/* Function:  p7_trace_GetDomainCount()
 *
 * Purpose:   Determine the number of hits in the trace <tr> -- that is,
 *            the number of times the trace traverses the model from
 *            B...E.  Return that number in <ret_ndom>.
 *            
 *            Done simply by counting the number of B states used in
 *            the trace.
 *            
 *            Only sensible on profile traces. Core traces have 1
 *            domain by definition.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_trace_GetDomainCount(const P7_TRACE *tr, int *ret_ndom)
{
  int z;
  int ndom = 0;

  if (tr->ndom > 0) 
    ndom = tr->ndom; /* if we already indexed the domains, we know the answer */
  else {
    for (z = 0; z < tr->N; z++)
      if (tr->st[z] == p7T_B) ndom++;
  }
  *ret_ndom = ndom;
  return eslOK;
}

/* Function:  p7_trace_GetStateUseCounts()
 *
 * Purpose:   Accumulate counts of each different state type in trace <tr>. 
 *
 *            <counts[]> is allocated for at least <p7T_NSTATETYPES>
 *            integers, indexed by statetype. Upon return,
 *            <counts[p7T_MG]> contains the number of glocal match states
 *            in the trace, for example.
 */
int
p7_trace_GetStateUseCounts(const P7_TRACE *tr, int *counts)
{
  int x,z;

  for (x = 0; x < p7T_NSTATETYPES; x++) counts[x] = 0;

  for (z = 0; z < tr->N; z++) {
    x = tr->st[z];
    if (x < 0 || x >= p7T_NSTATETYPES) ESL_EXCEPTION(eslEINVAL, "bad state type");
    counts[x]++;
  }
  return eslOK;
}

/* Function:  p7_trace_GetDomainCoords()
 *
 * Purpose:   Retrieve the bounds of domain alignment number <which> in
 *            traceback <tr>. <which> starts from 0. The total number
 *            of domains in a trace can be obtained from
 *            <p7_trace_GetDomainCount()>, or caller can just loop
 *            an increasing <which> count until <eslEOD> is returned.
 *            
 *            Start/end in the sequence are returned in <ret_i1>,
 *            <ret_i2>. Start/end in the model are returned in <ret_k1>, 
 *            <ret_k2>.
 *
 *            It only makes sense to call this function on profile
 *            traces.
 *            
 *            By local alignment bounds convention, the domain
 *            alignment is defined as bounded by match states, so <k1>
 *            and <k2> are the coords of the first and last match
 *            state (in range 1..M), and <i1> and <i2> are the coords
 *            of the residues aligned to those match states. Profiles
 *            do allow a Mk->DDD->E trailer; nonetheless, if such a
 *            trailer occurs, the k2 coord still refers to the last
 *            match state's coordinate. Note that such trailers would
 *            only occur in generated or sampled paths, not Viterbi
 *            paths; in Viterbi alignments with exit probabilities of
 *            1.0, the direct Mk->E path will always have higher
 *            probability than a Mk->DDD->E path.
 *
 * Returns:   <eslOK> on success, and the coords are returned.
 *            <eslEOD> if the trace doesn't contain a <which>'th
 *            domain, and the coords are all returned as 0.
 *            
 * Throws:    <eslEINVAL> if you stupidly pass a <which> less than 0;
 *            <eslECORRUPT> if something is grievously wrong with <tr>.           
 */
int
p7_trace_GetDomainCoords(const P7_TRACE *tr, int which,
			 int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2)
{
  int z;
  int status;

  if (which < 0) ESL_XEXCEPTION(eslEINVAL, "bad which < 0");

  if (tr->ndom) 		/* do we have an index? then this'll be easy */
    {
      if (which >= tr->ndom) { status = eslEOD; goto ERROR; }
      *ret_i1 = tr->sqfrom[which];
      *ret_i2 = tr->sqto[which];
      *ret_k1 = tr->hmmfrom[which];
      *ret_k2 = tr->hmmto[which];
      return eslOK;
    }

  /* else, the hard way.
   * skip z to one state past the which'th B state. 
   */
  for (z = 0; which >= 0 && z < tr->N; z++)
    if (tr->st[z] == p7T_B) which--;
  if (z == tr->N) { status = eslEOD; goto ERROR; }

  /* Find start coord on model and seq.
   * B->G->M1/D1 or B->L->Mk; z is currently on G/L.
   */
  z++;
  *ret_k1 = tr->k[z];		   /* z is on {MD}1 for glocal, MLk local   */
  while (tr->st[z] == p7T_DG) z++; 
  *ret_i1 = tr->i[z];		   /* z is on MGk for glocal, MLk for local */

  /* Skip ahead to the state before E. */
  while (tr->st[z+1] != p7T_E) z++;

  /* Find end coord on model, seq 
   * glocal: we're on {MD}m; local: we're on {MD}k (Dk only happens on suboptimal path)
   */
  *ret_k2 = tr->k[z];
  while (tr->st[z] == p7T_DG || tr->st[z] == p7T_DL) z--;  
  *ret_i2 = tr->i[z];
  return eslOK;

 ERROR:
  *ret_i1 = 0;
  *ret_i2 = 0;
  *ret_k1 = 0;
  *ret_k2 = 0;
  return status;
}
/*---------------- end, access routines -------------------------*/




/*****************************************************************
 * 3. Debugging tools.
 *****************************************************************/

/* Function:  p7_trace_DecodeStatetype()
 * Synopsis:  Convert an internal state type code to a string.
 *
 * Purpose:   Returns the state type in text, as a string.
 *            For example, <p7_trace_DecodeStatetype(p7T_S)>
 *            returns "S".
 *            
 * Throws:    an internal <eslEINVAL> exception if the code doesn't 
 *            exist, and returns <NULL>.           
 */
char *
p7_trace_DecodeStatetype(char st)
{
  switch (st) {
  case p7T_ML: return "ML";
  case p7T_MG: return "MG";
  case p7T_IL: return "IL";
  case p7T_IG: return "IG";
  case p7T_DL: return "DL";
  case p7T_DG: return "DG";
  case p7T_S:  return "S";
  case p7T_N:  return "N";
  case p7T_B:  return "B";
  case p7T_L:  return "L";
  case p7T_G:  return "G";
  case p7T_E:  return "E";
  case p7T_C:  return "C";
  case p7T_J:  return "J";
  case p7T_T:  return "T";
  default:     break;
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such statetype code %d", st);
  return NULL;
}

/* Function:  p7_trace_Validate()
 *
 * Purpose:   Validate the internal data in a trace structure <tr>
 *            representing an alignment of an HMM to a 
 *            digital sequence <sq>. The digital sequence may be either
 *            unaligned (usually) or aligned (in the case of "fake"
 *            tracebacks generated from an MSA during a
 *            model construction process). 
 *            
 *            We don't pass the HMM that the trace is associated with,
 *            because we might have constructed the trace during
 *            HMM construction when we don't have an HMM yet; but 
 *            we always have a digital sequence.
 *
 *            Intended for debugging/development/testing only.
 *            
 * Args:      tr     - trace to validate
 *            abc    - alphabet corresponding to sequence <sq>
 *            dsq    - digital sequence that <tr> is explaining
 *            errbuf - NULL, or an error message buffer allocated
 *                     for at least eslERRBUFSIZE chars.           
 *
 * Returns:   <eslOK> if trace appears fine.
 *            Returns <eslFAIL> if a problem is detected; if <errbuf> is
 *            provided (non-<NULL>), an informative message is formatted
 *            there to indicate the reason for the failure.
 */
int
p7_trace_Validate(const P7_TRACE *tr, const ESL_ALPHABET *abc, const ESL_DSQ *dsq, char *errbuf)
{
                                           /* X  ML MG IL IG DL DG S  N  B  L  G  E  C  J  T */
  static int is_kindexed[p7T_NSTATETYPES] = { 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; /* main states that have a k index */
  static int is_kbumped [p7T_NSTATETYPES] = { 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; /* states that advance k index on transition */
  static int is_memitter[p7T_NSTATETYPES] = { 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; /* states that emit on state */
  static int is_temitter[p7T_NSTATETYPES] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0 }; /* states that emit on transition */
  static int is_valid   [p7T_NSTATETYPES] = { 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0 }; /* valid in nonterminal positions 1..N-2 in trace */

  static int valid_transition[p7T_NSTATETYPES][p7T_NSTATETYPES] = {
    /*         -  ML MG IL IG DL DG S  N  B  L  G  E  C  J  T */
    /* -  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, 
    /* ML */ { 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, 
    /* MG */ { 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /* MG->E only if k==M */
    /* IL */ { 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /* IL->E only on construction, special case where MSA implies it */
    /* IG */ { 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, 
    /* DL */ { 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, 
    /* DG */ { 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /* DG->E only if k==M */
    /* S  */ { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, 
    /* N  */ { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0 }, 
    /* B  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0 }, 
    /* L  */ { 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* L->{DI}L only on construction, special case where MSA implies it */
    /* G  */ { 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* G->{MD} only for k=1 */
    /* E  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0 }, 
    /* C  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1 }, 
    /* J  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, 
    /* T  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* T is terminal, has no transitions */
  };
  int  z, i, k;			/* position in trace, sequence, model */

  /* minimum trace length is 8, S->N->B->L->Mk->E->C->T. If we don't have at least that,
   * we're definitely in trouble. Exception: N=0 can happen in the case of an "impossible"
   * trace.
   */
  if (tr == NULL)         return eslOK;
  if (tr->N == 0)         return eslOK;
  if (tr->N < 8)          ESL_FAIL(eslFAIL, errbuf, "trace is too short");
  if (tr->N > tr->nalloc) ESL_FAIL(eslFAIL, errbuf, "N of %d isn't sensible", tr->N);

  /* Verify known terminal states */
  if (tr->st[0]       != p7T_S) ESL_FAIL(eslFAIL, errbuf, "first state should be S");
  if (tr->st[1]       != p7T_N) ESL_FAIL(eslFAIL, errbuf, "second state should be N");
  if (tr->st[tr->N-2] != p7T_C) ESL_FAIL(eslFAIL, errbuf, "penultimate state should be C");
  if (tr->st[tr->N-1] != p7T_T) ESL_FAIL(eslFAIL, errbuf, "last state should be T");

  /* Main validation loop. */
  k = 0; 
  i = 1;
  for (z = 1; z < tr->N-1; z++)
    {
      /* on special case of an aligned dsq[] that includes gaps: increment i to past them */
      while (esl_abc_XIsGap(abc, dsq[i]) || esl_abc_XIsMissing(abc, dsq[i])) i++;

      /* Validate state type st[] */
      if (tr->st[z] < 0 || tr->st[z] >= p7T_NSTATETYPES) ESL_FAIL(eslFAIL, errbuf, "state type out of range, at %d", z);
      if (! is_valid[(int) tr->st[z]])                   ESL_FAIL(eslFAIL, errbuf, "invalid state type at %d", z);

      /* Validate state transition z-1,z */
      if (! valid_transition[(int) tr->st[z-1]][(int) tr->st[z]])             ESL_FAIL(eslFAIL, errbuf, "invalid transition, at %d", z);
      if (tr->st[z-1] == p7T_MG && tr->st[z] == p7T_E && tr->k[z-1] != tr->M) ESL_FAIL(eslFAIL, errbuf, "MG->E glocal exit only if k==M");
      if (tr->st[z-1] == p7T_DG && tr->st[z] == p7T_E && tr->k[z-1] != tr->M) ESL_FAIL(eslFAIL, errbuf, "DG->E glocal exit only if k==M");
      if (tr->st[z-1] == p7T_G  && tr->k[z]  != 1)                            ESL_FAIL(eslFAIL, errbuf, "G->{MD}k glocal entry only k==1");

      /* Validate state index k */
      if (is_kindexed[(int) tr->st[z]]) 
	{
	  if      (tr->st[z-1] == p7T_L)        k = tr->k[z]; /* on a L->Mk entry, trust k; else verify based on prev k */
	  else if (is_kbumped[(int) tr->st[z]]) k++;	      /* M,D states increment previous k; I's don't */

	  if (tr->k[z] < 1 || tr->k[z] > tr->M) ESL_FAIL(eslFAIL, errbuf, "invalid k[] at %d", z);
	  if (tr->k[z] != k)                    ESL_FAIL(eslFAIL, errbuf, "expected k doesn't match trace's k");
	}
      else
	{
	  k = 0;		/* reset expected k */
	  if (tr->k[z] != 0)                  ESL_FAIL(eslFAIL, errbuf, "invalid k[z] at %d", z);
	}

      /* Validate emission data, i[] and pp[] */
      if (is_memitter[(int) tr->st[z]] || (is_temitter[(int) tr->st[z]] && tr->st[z-1] == tr->st[z]))
	{
	  if (tr->i[z] < 1 || tr->i[z] > tr->L)                 ESL_FAIL(eslFAIL, errbuf, "invalid i[] at %d", z);
	  if (tr->pp && (tr->pp[z] < 0.0 || tr->pp[z] > 1.001)) ESL_FAIL(eslFAIL, errbuf, "invalid pp[] at %d", z);
	  if (tr->i[z] != i)                                    ESL_FAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	  i++;
	}
      else
	{
	  if (tr->i[z] != 0)               ESL_FAIL(eslFAIL, errbuf, "invalid i[] at %d", z);
	  if (tr->pp && tr->pp[z] != 0.0f) ESL_FAIL(eslFAIL, errbuf, "invalid pp[] at %d", z);
	}
    }

  /* Trace should have accounted for all residues in the dsq */
  for (; dsq[i] != eslDSQ_SENTINEL; i++) 
    if (! (esl_abc_XIsGap(abc, dsq[i]) || esl_abc_XIsMissing(abc, dsq[i])))
      ESL_FAIL(eslFAIL, errbuf, "trace didn't account for all residues in dsq");

  /* i should be sitting on dsq[L+1] sentinel right now */
  if (i != tr->L+1) ESL_FAIL(eslFAIL, errbuf, "L=%d, but i went to %d\n", tr->L, i-1);
  return eslOK;
}


/* Function:  p7_trace_Dump()
 * Synopsis:  Dump a trace's internals to a stream.
 *
 * Purpose:   Dump the internals of trace <tr> to stream <fp>,
 *            for debugging/examination.
 *            
 *            If <tr> is <NULL>, just print "[null trace]". We
 *            tolerate dumping null traces because
 *            <p7_trace_FauxFromMSA()> may normally return them.
 *
 *            See <p7_trace_DumpAnnotated()> for a more verbose
 *            and useful version. This version is used when we
 *            don't necessarily have a profile and/or sequence
 *            to go with a trace.
 */
int 
p7_trace_Dump(FILE *fp, const P7_TRACE *tr)
{
  int z;

  if (tr == NULL) { fprintf(fp, "[null trace]\n");                     return eslOK; }
  if (tr->N == 0) { fprintf(fp, "[no trace: all paths impossible]\n"); return eslOK; }

  fprintf(fp, "z     st   k     i  \n");  
  fprintf(fp, "----- -- ----- -----\n"); 

  for (z = 0; z < tr->N; z++)
    fprintf(fp, "%5d %2s %5d %5d\n", z, p7_trace_DecodeStatetype(tr->st[z]), tr->k[z], tr->i[z]);

  fprintf(fp, "\n#          M = %d\n", tr->M);
  fprintf(fp,   "#          L = %d\n", tr->L);
  fprintf(fp,   "# allocation = %d\n", tr->nalloc);
  return eslOK;
}

/* Function:  p7_trace_DumpAnnotated()
 *
 * Purpose:   Dumps internals of a traceback structure <tr> to <fp>,
 *            annotating residues, scores, and posterior probabilities using
 *            the profile <gm> and the sequence <dsq> that correspond
 *            to this trace.
 *            
 *            If <tr> is <NULL>, just print "[null trace]". We
 *            tolerate dumping null traces because
 *            <p7_trace_FauxFromMSA()> may normally return them.
 *            
 * Args:      fp   - stream to dump to (often stdout)
 *            tr   - trace to dump (may be NULL)
 *            gm   - score profile corresponding to trace
 *            dsq  - digitized seq corresponding to trace        
 *
 * Returns:   <eslOK> on success.
 */
int
p7_trace_DumpAnnotated(FILE *fp, const P7_TRACE *tr, const P7_PROFILE *gm, const ESL_DSQ *dsq)
{
  float accuracy = 0.0;
  float sc;
  float esc, tsc;
  int   z;

  if (tr == NULL) { fprintf(fp, "[null trace]\n");                     return eslOK; }
  if (tr->N == 0) { fprintf(fp, "[no trace: all paths impossible]\n"); return eslOK; }

  fprintf(fp, "#  z   st   k     i   x_i  transit  emission postprob\n");
  fprintf(fp, "#----- -- ----- ----- ---  -------- -------- --------\n");

  for (z = 0; z < tr->N-1; z++)	/* not including T state at end: */
    {
      tsc = p7_profile_GetT(gm, tr->st[z], tr->k[z], tr->st[z+1], tr->k[z+1]);

      esc = 0.;
      if (tr->i[z]) {
	if      (tr->st[z] == p7T_ML || tr->st[z] == p7T_MG) esc = P7P_MSC(gm, tr->k[z], dsq[tr->i[z]]);
	else if (tr->st[z] == p7T_IL || tr->st[z] == p7T_IG) esc = P7P_ISC(gm, tr->k[z], dsq[tr->i[z]]);
      }

      fprintf(fp, "%5d %2s %5d %5d  %c  %8.4f %8.4f %8.4f\n",
	      z,
	      p7_trace_DecodeStatetype(tr->st[z]),
	      tr->k[z],
	      tr->i[z],
	      (tr->i[z] ? gm->abc->sym[dsq[tr->i[z]]] : '-'),
	      tsc, 
	      esc, 
	      (tr->pp ? tr->pp[z] : 0.0f));

      accuracy += (tr->pp ? tr->pp[z] : 0.0f);
    }
  p7_trace_Score(tr, dsq, gm, &sc);

  fprintf(fp, "%5d %2s\n", z, p7_trace_DecodeStatetype(tr->st[z])); /* T state */
  fprintf(fp, "                                   -------- --------\n");
  fprintf(fp, "                          total:   %8.4f %8.4f\n", sc, accuracy);

  fprintf(fp, "\n#          M = %d\n", tr->M);
  fprintf(fp,   "#          L = %d\n", tr->L);
  fprintf(fp,   "# allocation = %d\n", tr->nalloc);
  return eslOK;
}


/* Function:  p7_trace_Compare()
 * Synopsis:  Compare two traces for identity
 *
 * Purpose:   Compare two tracebacks; return <eslOK> if they
 *            are identical, throw <eslFAIL> if not.
 *            
 *            If posterior probability annotation is present in 
 *            both traces, they are compared using <esl_FCompare()>
 *            and a relative tolerance of <pptol>.
 *            
 *            If domain indices are present in both traces,
 *            the two indexes are compared.
 */
int
p7_trace_Compare(P7_TRACE *tr1, P7_TRACE *tr2, float pptol)
{
  int z,d;
  
  if (tr1->N != tr2->N) ESL_EXCEPTION(eslFAIL, "traces' N differ");
  if (tr1->M != tr2->M) ESL_EXCEPTION(eslFAIL, "traces' M differ");
  if (tr1->L != tr2->L) ESL_EXCEPTION(eslFAIL, "traces' L differ");
  
  /* Main data in the trace */
  for (z = 0; z < tr1->N; z++)
    {
      if (tr1->st[z] != tr2->st[z]) ESL_EXCEPTION(eslFAIL, "traces' state types differ at %d", z);
      if (tr1->k[z]  != tr2->k[z])  ESL_EXCEPTION(eslFAIL, "traces' k indices differ at %d", z);
      if (tr1->i[z]  != tr2->i[z])  ESL_EXCEPTION(eslFAIL, "traces' i indices differ at %d", z);
      if (tr1->pp != NULL && tr2->pp != NULL)
	if (esl_FCompare(tr1->pp[z], tr2->pp[z], pptol) != eslOK)  /* comparison of 0.0 for positions without pp will succeed */
	  ESL_EXCEPTION(eslFAIL, "traces' posterior probs differ at %d", z);
    }

  /* Optional domain index */
  if (tr1->ndom > 0 && tr2->ndom > 0)
    {
      if (tr1->ndom != tr2->ndom) ESL_EXCEPTION(eslFAIL, "traces' domain table # differ");

      for (d = 0; d < tr1->ndom; d++)
	{
	  if (tr1->tfrom[d]   != tr2->tfrom[d])    ESL_EXCEPTION(eslFAIL, "traces' tfrom differs, domain %d",     d);
	  if (tr1->tto[d]     != tr2->tto[d])      ESL_EXCEPTION(eslFAIL, "traces' tto differs, domain %d",       d);
	  if (tr1->sqfrom[d]  != tr2->sqfrom[d])   ESL_EXCEPTION(eslFAIL, "traces' sqfrom differs, domain %d",    d);
	  if (tr1->sqto[d]    != tr2->sqto[d])     ESL_EXCEPTION(eslFAIL, "traces' sqto differs, domain %d",      d);
	  if (tr1->hmmfrom[d] != tr2->hmmfrom[d])  ESL_EXCEPTION(eslFAIL, "traces' hmmfrom differs, domain %d",   d);
	  if (tr1->hmmto[d]   != tr2->hmmto[d])    ESL_EXCEPTION(eslFAIL, "traces' hmmto differs, domain %d",     d);
	  //if (tr1->anch[d]    != tr2->anch[d])     ESL_EXCEPTION(eslFAIL, "traces' anchors differ for domain %d", d);  // I think <anch> is only used by mass trace.
	}
    }
  return eslOK;
}


/* Function:  p7_trace_CompareLoosely()
 * Synopsis:  Weaker version of p7_trace_Compare().
 *
 * Purpose:   When we say "reference and sparse traces must be
 *            identical" in some of the sparse DP unit tests, we don't
 *            really mean it.  It is possible to have two (or more)
 *            possible traces with exactly the same score, such that
 *            any of them are valid Viterbi paths. It is possible for
 *            sparse trace to find one, and reference trace to find
 *            another.
 * 
 *            The most common example is on a single-residue
 *            alignment. Suppose the reference trace has state ML31
 *            aligned to residue Y64, and that's the only aligned
 *            residue, with all other residues explained by N/C. That
 *            is, SN...NB->ML31->EC...CT. Now any and all other Y
 *            residues in the target sequence can also be aligned to
 *            ML31, necessarily receiving the same emission score, and
 *            the trace necessarily receives the same overall score.
 * 
 *            Of course, this is a pathological case. I didn't expect
 *            alternative traces with identical scores, for
 *            position-specific floating-point scores, but an example
 *            like that above showed up in a unit test failure. If
 *            <p7_trace_Compare()> is used in sparse DP unit tests
 *            that are subject to this issue, they will sometimes fail.
 * 
 *            Thus the <p7_trace_CompareLoosely() variant, which
 *            allows emitting M/I states to emit exactly the same
 *            subsequence in the target: that is, it checks that the
 *            residue identities match (as opposed to more stringently
 *            requiring the i indices to match), for all M/I states.
 */
int
p7_trace_CompareLoosely(P7_TRACE *tr1, P7_TRACE *tr2, ESL_DSQ *dsq)
{
  int z1 = 0;			/* position in <tr1>: 0..tr1->N-1 */
  int z2 = 0;			/* position in <tr2>: 0..tr2->N-1 */

  for (z1 = 0; z1 < tr1->N; z1++)
    {
      if (p7_trace_IsM(tr1->st[z1]) || p7_trace_IsI(tr1->st[z1]))
	{
	  while (z2 < tr2->N-1 && tr1->st[z1] != tr2->st[z2]) z2++;
	  
	  if (tr1->st[z1]     != tr2->st[z2])      return eslFAIL;
	  if (tr1->k[z1]      != tr2->k[z2])       return eslFAIL;
	  if (dsq[tr1->i[z1]] != dsq[tr2->i[z2]])  return eslFAIL;
	  z2++;
	}
    }
  return eslOK;
}



/* Function:  p7_trace_Score()
 *
 * Purpose:   Score path <tr> for digital target sequence <dsq> 
 *            using profile <gm>. Return the raw lod score (nats) in
 *            <*ret_sc>.
 *            
 *            If <tr> is empty (an 'impossible' trace, with <tr->N=0>),
 *            then <*ret_sc = -eslINFINITY>. This can arise for example
 *            when there is no possible path at all that can generate
 *            a sequence, hence the sequence correctly scores <-eslINFINITY>;
 *            by convention, traces in this situation have <tr->N=0>.
 *            
 * Args:      tr     - traceback path to score
 *            dsq    - digitized sequence
 *            gm     - score profile
 *            ret_sc - RETURN: lod score of trace <tr>
 *
 * Returns:   <eslOK> on success, and <*ret_sc> contains the
 *            lod score for the trace.
 *
 * Throws:    <eslEINVAL> if something's wrong with the trace.
 *            Now <*ret_sc> is returned as $-\infty$.
 *            
 * Note:      On numerical roundoff error: Because we sum terms in the
 *            forward direction, exactly the same order as the Viterbi
 *            algorithm does, the trace score of an optimal Viterbi
 *            trace should be identical to the Viterbi score
 *            calculated by the Viterbi algorithm. To calculate the
 *            trace score with minimal roundoff error, see
 *            <p7_trace_ScoreKahan()>. To calculate the trace score in
 *            the backwards direction, see <p7_trace_ScoreBackwards()>.
 *
 * Xref:      [SRE:J13/121 22-Jun-14]: numerical accuracy considerations
 */
int 
p7_trace_Score(const P7_TRACE *tr, const ESL_DSQ *dsq, const P7_PROFILE *gm, float *ret_sc)
{
  float sc = (tr->N ? 0.0f : -eslINFINITY);
  int   z;
  int   xi;

  for (z = 0; z < tr->N-1; z++) {
    xi = dsq[tr->i[z]];
    
    if      (tr->st[z] == p7T_ML || tr->st[z] == p7T_MG) sc += P7P_MSC(gm, tr->k[z], xi);
    else if (tr->st[z] == p7T_IL || tr->st[z] == p7T_IG) sc += P7P_ISC(gm, tr->k[z], xi);

    sc += p7_profile_GetT(gm, tr->st[z], tr->k[z], tr->st[z+1], tr->k[z+1]);
  }
  *ret_sc = sc;
  return eslOK;
}


/* Function:  p7_trace_ScoreKahan()
 * Synopsis:  Calculate trace score with high numerical accuracy.
 *
 * Purpose:   Same as <p7_trace_Score()>, except using Kahan compensated
 *            summation to minimize roundoff error accumulation.
 *            
 *            Because Viterbi itself does not use compensated
 *            summation, trace score of an optimal Viterbi trace using
 *            <p7_trace_ScoreKahan()> is not necessarily identical to
 *            the Viterbi score; they will usually differ by a small
 *            amount of roundoff error accumulation.
 *            <p7_trace_ScoreKahan()> is used in development to
 *            characterize that numerical error.
 *
 * Xref:      [SRE:J13/121 22-Jun-14]: incept
 */
int
p7_trace_ScoreKahan(const P7_TRACE *tr, const ESL_DSQ *dsq, const P7_PROFILE *gm, float *ret_sc)
{
  float  sc;		/* total lod score   */
  int    z;             /* position in tr */
  int    xi;		/* digitized symbol in dsq */
  float  y,t,c;		/* Kahan compensated summation */

  sc = (tr->N ? 0.0f : -eslINFINITY);
  c  = 0.0;
  for (z = 0; z < tr->N-1; z++) {
    xi = dsq[tr->i[z]];

    if      (tr->st[z] == p7T_ML || tr->st[z] == p7T_MG) { y = P7P_MSC(gm, tr->k[z], xi) - c; t = sc + y; c = (t-sc)-y; sc = t; }
    else if (tr->st[z] == p7T_IL || tr->st[z] == p7T_IG) { y = P7P_ISC(gm, tr->k[z], xi) - c; t = sc + y; c = (t-sc)-y; sc = t; }

    y = p7_profile_GetT(gm, tr->st[z], tr->k[z], tr->st[z+1], tr->k[z+1]) - c; t = sc + y; c = (t-sc)-y; sc = t;
  }

  *ret_sc = sc;
  return eslOK;
}

/* Function:  p7_trace_ScoreBackwards()
 * Synopsis:  Calculate trace score, summing terms in reverse order.
 *
 * Purpose:   Same as <p7_trace_Score()> but we sum terms in reverse
 *            order.  Result will generally differ from that of
 *            <p7_trace_Score()> because of numerical roundoff. This
 *            routine is here for characterizing that numerical
 *            roundoff error during code development.
 *
 * Xref:      [SRE:J13/121 22-Jun-14]: incept
 */
int
p7_trace_ScoreBackwards(const P7_TRACE *tr, const ESL_DSQ *dsq, const P7_PROFILE *gm, float *ret_sc)
{
  float sc = (tr->N ? 0.0f : -eslINFINITY);
  int   z;
  int   xi;
 
  for (z = tr->N-2; z >= 0; z--) {  
    xi = dsq[tr->i[z]];

    if      (tr->st[z] == p7T_ML || tr->st[z] == p7T_MG) sc += P7P_MSC(gm, tr->k[z], xi);
    else if (tr->st[z] == p7T_IL || tr->st[z] == p7T_IG) sc += P7P_ISC(gm, tr->k[z], xi);

    sc += p7_profile_GetT(gm, tr->st[z], tr->k[z], tr->st[z+1], tr->k[z+1]);
  }

  *ret_sc = sc;
  return eslOK;
}
    

/* Function:  p7_trace_ScoreDomain()
 * Synopsis:  Return a per-domain trace score for a selected domain.
 *
 * Purpose:   For an indexed trace <tr> of sequence <dsq> aligned to
 *            model <gm>, find domain number <which> (0..tr->ndom-1)
 *            and calculate a "per-domain" score for this domain.
 *            Return this lod score (in nats) in <*ret_sc>. The caller
 *            still needs to subtract the null model score and convert
 *            to bits.
 *            
 *            The per-domain score is the score that you'd get if this
 *            were the only domain in the sequence. If there is only
 *            one domain, it is equal to the overall trace score.
 *            
 *            The score respects the given configuration of <gm>
 *            (whether it is in unihit or multihit mode); the
 *            per-domain score is not necessarily w.r.t. unihit mode
 *            configuration.
 *
 * Args:      tr     - trace
 *            dsq    - sequence (digital, 1..tr->L)
 *            gm     - profile
 *            which  - which domain to score, 0..tr->ndom-1
 *            ret_sc - RETURN: per-domain score in nats
 *
 * Returns:   <eslOK> on success, and <*ret_sc> is the per-domain score
 *            in nats.
 *
 * Throws:    <eslEINVAL> if the trace isn't indexed, or if it has
 *            no domain number <which>. Now <*ret_sc> is undefined.
 */
int
p7_trace_ScoreDomain(P7_TRACE *tr, ESL_DSQ *dsq, P7_PROFILE *gm, int which, float *ret_sc)
{
  float sc = 0.0;
  int   xi;
  int   z;

  if (! tr->ndom)        ESL_EXCEPTION(eslEINVAL, "p7_trace_ScoreDomain() requires an indexed trace structure");
  if (which >= tr->ndom) ESL_EXCEPTION(eslEINVAL, "<which> exceeds number of domains");

  /* S->N, N->N*(i-1), N->B ;  E->C, C->C*(L-j), C->T */
  sc += (float)(tr->sqfrom[which]-1)     * gm->xsc[p7P_N][p7P_LOOP] + gm->xsc[p7P_N][p7P_MOVE];
  sc += (float)(tr->L - tr->sqto[which]) * gm->xsc[p7P_C][p7P_LOOP] + gm->xsc[p7P_C][p7P_MOVE] + gm->xsc[p7P_E][p7P_MOVE];
  
  for (z = tr->tfrom[which]; z < tr->tto[which]; z++)
    {
      xi = dsq[tr->i[z]];

      if      (tr->st[z] == p7T_ML || tr->st[z] == p7T_MG) sc += P7P_MSC(gm, tr->k[z], xi);
      else if (tr->st[z] == p7T_IL || tr->st[z] == p7T_IG) sc += P7P_ISC(gm, tr->k[z], xi);

      sc += p7_profile_GetT(gm, tr->st[z], tr->k[z], tr->st[z+1], tr->k[z+1]);
    }
  *ret_sc = sc;
  return eslOK;
}


/* Function:  p7_trace_GetExpectedAccuracy()
 * Synopsis:  Returns the sum of the posterior residue decoding probs.
 */
float
p7_trace_GetExpectedAccuracy(const P7_TRACE *tr)
{
  float accuracy = 0.0;
  int   z;

  for (z = 0; z < tr->N; z++)
    accuracy += tr->pp[z];
  return accuracy;
}
/*------------------ end, debugging tools -----------------------*/



/*****************************************************************
 * 4. Visualization tools.
 *****************************************************************/

/* Function:  p7_trace_PlotDomainInference()
 * Synopsis:  Plot a discrete trace in the same style as we plot decoding matrix data.
 *
 * Purpose:   Write a "domain inference plot" in xmgrace XY format to
 *            the open stream <ofp>, for indexed trace <tr>, for sequence
 *            region <ia>..<ib>.
 *            
 *            If you want the whole sequence, pass <ia=1>,<ib=tr->L>.
 *
 *            Trace <tr> must be indexed with <p7_trace_Index()>.
 *            
 *            The output graph is in the same style as, for example,
 *            <p7_refmx_PlotDomainInference()>, the plots we use to
 *            study posterior decoding of domain position. Purpose is
 *            to enable comparing discrete traces (optimal or sampled)
 *            to posterior decoding results.
 *
 * Args:      ofp  : open output stream for writing xmgrace plot
 *            tr   : trace to plot
 *            ia   : start of plot, in sequence coords (1, perhaps)
 *            ib   : end of plot, in seq coords (tr->L, perhaps)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_trace_PlotDomainInference(FILE *ofp, const P7_TRACE *tr, int ia, int ib)
{
  int   i;
  int   d         = 0;
  float tr_height = 1.2;

  /* Special case: tr->N==0 means there's no trace at all.
   * This case gets handled correctly because ndom==0, so sets 0,1,2,3
   * will be flatline at 0, and there's no individual domain sets.
   */

  /* Set #0 is P(homology), which for a trace is discrete 1.0 or 0.0, in/out of domains */
  d = 0;
  for (i = 1; i <= tr->L; i++)
    {
      if (i >= ia && i <= ib)
	{
	  if (d >= tr->ndom || i < tr->sqfrom[d]) fprintf(ofp, "%-6d 0.0\n", i); /* nonhomologous position */
	  else			                  fprintf(ofp, "%-6d 1.0\n", i); /* inside homology region */
	}
      if (d < tr->ndom && i == tr->sqto[d]) d++;
    }
  fprintf(ofp, "&\n");

  /* Set #1 is P(B), the begin state. 1.0 at -1 position before a domain */
  d = 0;
  for (i = 0; i <= tr->L; i++)
    {
      if (i >= ia && i <= ib)
	{
	  if (d >= tr->ndom || i != tr->sqfrom[d]-1) fprintf(ofp, "%-6d 0.0\n", i); 
	  else                                       fprintf(ofp, "%-6d 1.0\n", i); /* begin = 1.0 at start-1 position */
	}
      if (d < tr->ndom && i == tr->sqto[d]) d++;
    }
  fprintf(ofp, "&\n");

  /* Set #2 is P(E), the end state. 1.0 at last position in domain */
  d = 0;
  for (i = 0; i <= tr->L; i++)
    {
      if (i >= ia && i <= ib)
	{
	  if (d >= tr->ndom || i != tr->sqto[d]) fprintf(ofp, "%-6d 0.0\n", i); 
	  else                                   fprintf(ofp, "%-6d 1.0\n", i); 
	}
      if (d < tr->ndom && i == tr->sqto[d]) d++;
    }
  fprintf(ofp, "&\n");

  /* Set #3 is max rather than sum of pp's for M,I states; here, for
   * Viterbi, same as set #0, and we include it only for completeness,
   * to exactly match the style of the main routine we're emulating,
   * <p7_refmx_PlotDomainInference()>.
   */
  d = 0;
  for (i = 1; i <= tr->L; i++) {
    if (i >= ia && i <= ib) {
      if (d >= tr->ndom || i < tr->sqfrom[d]) fprintf(ofp, "%-6d 0.0\n", i); /* nonhomologous position */
      else			              fprintf(ofp, "%-6d 1.0\n", i); /* inside homology region */
    }
    if (d < tr->ndom && i == tr->sqto[d]) d++;
  }
  fprintf(ofp, "&\n");

  /* Remaining sets are horiz lines, representing individual domains appearing in the optional <tr> */
  for (d = 0; d < tr->ndom; d++)
    {
      fprintf(ofp, "%-6d %.5f\n", tr->sqfrom[d], tr_height);
      fprintf(ofp, "%-6d %.5f\n", tr->sqto[d],   tr_height);
      fprintf(ofp, "&\n");
    }

  return eslOK;
}


/* Function:  p7_trace_PlotHeatMap()
 * Synopsis:  Plot heat map visualization of single trace in DP matrix, in PostScript
 *
 * Purpose:   Plot a heat map representation of a single trace <tr> in a
 *            window of its DP matrix, writing it to open stream <ofp>
 *            in PostScript format. Window is sequence position
 *            <ia>..<ib>, and model node <ka>..<kb>. To plot the whole
 *            matrix, pass <ia=1>, <ib=tr->L>, <ka=1>, <kb=tr->M>.
 *            
 *            Obviously this "heat map" is a trivialized one,
 *            discretized at 0.0/1.0 values, for the discrete
 *            trace. The purpose is to be able to compare discrete
 *            traces to the results of posterior decoding, plotted
 *            with <p7_refmx_PlotHeatMap()>.
 *            
 *            See <p7_refmx_PlotHeatMap()> for more information, and
 *            perhaps also <esl_dmatrix_PlotHeatMap()>.
 *            
 *            Note that the plot is constructed with the sequence
 *            coords running along the vertical axis from bottom to
 *            top (NOT top to bottom), and model coords running along
 *            the horizontal axis from left to right. This enables a
 *            rotate-right of the PostScript to make the plot have the
 *            sequence run left to right and model top to bottom,
 *            which has become sort of standard in how I show these
 *            matrices in slides and figures, but is different from
 *            how I think about them in the code.
 *
 * Args:      ofp : open stream for writing PostScript output
 *            tr  : trace to visualize in a DP matrix
 *            ia  : seq coord start (1..L; perhaps 1)
 *            ib  : seq coord end   (1..L; perhaps tr->L)
 *            ka  : model coord start (1..M; perhaps 1)
 *            kb  : model coord end   (1..M; perhaps tr->M)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_trace_PlotHeatMap(FILE *ofp, P7_TRACE *tr, int ia, int ib, int ka, int kb)
{
  ESL_DMATRIX *dmx  = NULL;
  int          nrow = ib - ia + 1;
  int          ncol = kb - ka + 1;
  int          i,k,z;
  int          status;
  
  /* Copy to an <ESL_DMATRIX>, which has heat mapping tools.
   * With a discrete trace, cell values are either 1.0 or 0.0.
   * Init the whole matrix to 0, then traverse the (linear)
   * trace and write 1.0 in every cell it visits.
   * Remember that i[] coord of D states is 0 in the trace.
   */
  if ((dmx = esl_dmatrix_Create(nrow, ncol)) == NULL) { status = eslEMEM; goto ERROR; }
  esl_dmatrix_SetZero(dmx);

  for (z = 0; z < tr->N; z++)
    {
      i = (tr->i[z] ? tr->i[z] : i); /* last i we emitted: coord we'll use for D states */
      k = tr->k[z];

      if (p7_trace_IsMain(tr->st[z]) ) /* any {MDI}{LG} */
	dmx->mx[ib-i][k-ka] = 1.0;
    }

  /* Now plot it using dmatrix heatmapping */
  if ((status = esl_dmatrix_PlotHeatMap(ofp, dmx, 0.0, 1.0)) != eslOK) goto ERROR;

  esl_dmatrix_Destroy(dmx);
  return eslOK;

 ERROR:
  if (dmx) esl_dmatrix_Destroy(dmx);
  return status;
}



/*****************************************************************
 * 5. Creating traces by DP traceback
 *****************************************************************/

/* Function:  p7_trace_Append()
 * Synopsis:  Add an element (state/residue) to a growing trace.
 *
 * Purpose:   Adds an element to a trace <tr> that is growing
 *            left-to-right. The element is defined by a state type
 *            <st> (such as <p7T_ML>); a node index <k> (1..M for
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
 *            plans to call <p7_trace_Reverse()>. 
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
p7_trace_Append(P7_TRACE *tr, char st, int k, int i)
{
  int status;

  if ((status = p7_trace_Grow(tr)) != eslOK) return status;

  switch (st) {
    /* Emit-on-transition states: */
  case p7T_N: case p7T_C: case p7T_J: 
    tr->i[tr->N] = ( (tr->st[tr->N-1] == st) ? i : 0);
    tr->k[tr->N] = 0;
    break;

    /* Nonemitting states, outside main model: */
  case p7T_S: case p7T_B: case p7T_G: case p7T_L: case p7T_E: case p7T_T:
    tr->i[tr->N] = 0; 
    tr->k[tr->N] = 0; 
    break;

    /* Nonemitting, but in main model (k valid) */
  case p7T_DL: case p7T_DG: 
    tr->i[tr->N] = 0; 
    tr->k[tr->N] = k; 
    break;

    /* Emitting states, with valid k position in model: */
  case p7T_ML: case p7T_MG: case p7T_IL: case p7T_IG:
    tr->i[tr->N] = i;
    tr->k[tr->N] = k; 
    break;

  default: 
    ESL_EXCEPTION(eslEINVAL, "no such state; can't append");
  }

  tr->st[tr->N] = st;
  tr->N++;
  return eslOK;
}

/* Function:  p7_trace_AppendWithPP()
 * Synopsis:  Add element to growing trace, with posterior probability.
 *
 * Purpose:   Same as <p7_trace_Append()>, but also records a posterior
 *            probability estimate for emitted residues. For nonemitting
 *            states, <pp> must be passed as zero.
 */
int
p7_trace_AppendWithPP(P7_TRACE *tr, char st, int k, int i, float pp)
{
  int status;

  if ((status = p7_trace_Append(tr, st, k, i)) != eslOK) return status;
  tr->pp[tr->N-1] = pp;
  return eslOK;
}

/* Function: p7_trace_Reverse()
 * Synopsis: Reverse the arrays in a traceback structure.
 * 
 * Purpose:  Reverse the arrays in a traceback structure.  Tracebacks
 *           from DP algorithms are collected backwards, and they call this
 *           function when they're done.
 *           
 *           At least for now, this invalidates any domain index
 *           table, if it exists. The expectd order of invocation is
 *           to create the traceback backwards, <Reverse()> it, then
 *           <Index()> it.
 *           
 * Args:     tr - the traceback to reverse. tr->N must be set.
 *                
 * Return:   <eslOK> on success; <tr> is modified.
 */                
int
p7_trace_Reverse(P7_TRACE *tr)
{
  int    z;
  int    tmp;
  float  tmpf;

  /* For emit-on-transition states N,C,J, traces always obey the
   * C-,Cx,Cx,Cx convention even when they were constructed backwards;
   * so we make them Cx,Cx,Cx,C- by pulling residues backwards by one,
   * just before reversing them. (Other ways of doing this would be
   * fine too.
   */
  for (z = 0; z < tr->N; z++)
    {
      if ( (tr->st[z] == p7T_N && tr->st[z+1] == p7T_N) ||
	   (tr->st[z] == p7T_C && tr->st[z+1] == p7T_C) ||
	   (tr->st[z] == p7T_J && tr->st[z+1] == p7T_J))
	{
	  if (tr->i[z] == 0 && tr->i[z+1] > 0) 
	    { 
	      tr->i[z]   = tr->i[z+1]; 
	      tr->i[z+1] = 0; 
	      if (tr->pp != NULL) {
		tr->pp[z]   = tr->pp[z+1];
		tr->pp[z+1] = 0.0;
	      }
	    }
	}
    }

  /* Reverse the trace in place. */
  for (z = 0; z < tr->N/2; z++)
    {
      tmp = tr->st[tr->N-z-1];  tr->st[tr->N-z-1] = tr->st[z];   tr->st[z] = tmp;
      tmp = tr->k[tr->N-z-1];   tr->k[tr->N-z-1]  = tr->k[z];    tr->k[z]  = tmp;
      tmp = tr->i[tr->N-z-1];   tr->i[tr->N-z-1]  = tr->i[z];    tr->i[z]  = tmp;
      if (tr->pp != NULL) {
	tmpf = tr->pp[tr->N-z-1];   tr->pp[tr->N-z-1]  = tr->pp[z];    tr->pp[z]  = tmpf;
      }
    }
  /* don't worry about the middle residue in odd-length N, since we're in-place  */
  return eslOK;
}


/* Function:  p7_trace_Index()
 * Synopsis:  Internally index the domains in a trace.
 *
 * Purpose:   Create an internal index of the domains in <tr>.  This
 *            indexing makes individual calls to
 *            <GetDomainCount()> and <GetDomainCoords()> more
 *            efficient, and it is a necessary prerequisite for
 *            creating alignments of any individual domains in a
 *            multidomain trace with <p7_alidisplay_Create()>.
 *            
 *            It's possible for alignment routines to result in empty
 *            tracebacks, where there is no possible path at all
 *            (Viterbi score -inf); empty tracebacks have <tr->N=0>. 
 *            In this case, <p7_trace_Index()> sets <tr->ndom> to 0.
 *
 * Returns:   <eslOK> on success. 
 *
 * Throws:    <eslEMEM> on allocation failure, in which case the
 *            data in the trace is still fine, but the domain index
 *            table isn't constructed.
 */
int
p7_trace_Index(P7_TRACE *tr)
{
  int   z,z2;
  float best_pp;
  int   anchor;
  int   status;

  tr->ndom = 0;
  z        = 0;
  while (z < tr->N && tr->st[z] != p7T_T)
    {
      /* Find start coord on model and seq.
       * B->G->M1/D1 or B->L->Mk; z is currently on B.
       */
      if (tr->st[z] == p7T_B)
	{
	  if ((status = p7_trace_GrowIndex(tr)) != eslOK) return status;
	  
	  tr->tfrom[tr->ndom]   = z;
	  z += 2;		/* skip B,{GL}; now we're either on M1/D1 for glocal, or MLk for local */

	  tr->hmmfrom[tr->ndom] = tr->k[z];
	  if (tr->pp) { best_pp = tr->pp[z]; anchor = z; }

	  while (tr->st[z] == p7T_DG) {
	    if (tr->pp && tr->pp[z] > best_pp) { best_pp = tr->pp[z]; anchor = z; }
	    z++;
	  }

	  tr->sqfrom[tr->ndom]  = tr->i[z];
	  if (tr->pp && tr->pp[z] > best_pp) { best_pp = tr->pp[z]; anchor = z; }

	  /* skip ahead to state before the E */
	  while (tr->st[z+1] != p7T_E) { 
	    if (tr->pp && tr->pp[z] > best_pp) { best_pp = tr->pp[z]; anchor = z; }
	    z++;
	  }

	  /* G... Mm/Dm -> E or L... Mk/Dk->E */
	  tr->hmmto[tr->ndom] = tr->k[z];
	  z2 = z;
	  while (tr->st[z2] == p7T_DG || tr->st[z2] == p7T_DL) z2--;
	  tr->sqto[tr->ndom] = tr->i[z2];
	  tr->tto[tr->ndom]   = ++z;
	  /* z is now on E state */
	  if (tr->pp) { tr->anch[tr->ndom] = anchor; }
	  tr->ndom++;
	}
      z++;
    }
  return eslOK;
}
/*----------- end, creating traces by DP traceback ---------------*/


/*****************************************************************
 * 6. Creating faux traces from MSAs
 *****************************************************************/

/* Function:  p7_trace_FauxFromMSA()
 * Synopsis:  Create array of faux tracebacks from an existing MSA.
 *
 * Purpose:   Given an existing <msa> and an array <matassign> that
 *            flags the alignment columns that are assigned to consensus
 *            match states (matassign[1..alen] = 1|0); create an array
 *            of faux traces <tr[0..msa->nseq-1]>. <optflags> controls 
 *            optional behavior; it can be <p7_DEFAULT> or <p7_MSA_COORDS>,
 *            as explained below. <matassign[]> must mark at least 
 *            one consensus column.
 *            
 *            "Faux" tracebacks may contain transitions that are
 *            illegal in H3, but are implied by an input MSA: they do
 *            not necessarily pass a Validate(). Read on for (much)
 *            more detail.
 *            
 *            This creates single-domain profile traces. Any flanking
 *            insertions (outside the first/last consensus column) are
 *            assigned to NN/CC transitions. No I0/Im states and no
 *            J states occur.
 *            
 *            Alignments are assumed to be glocal unless otherwise
 *            indicated to be local. Currently the convention for
 *            flagging a local alignment is that the leading/trailing
 *            gaps are converted to missing data symbols; see
 *            <esl_msa_MarkFragments()> for example. (But this is a
 *            problematic convention; it would be better to directly
 *            annotate the MSA somehow, instead of overloading the ~
 *            character this way.)
 *            
 *            By default (<optflags = p7_DEFAULT>), the <i> coordinate
 *            in the faux tracebacks is <1..L>, relative to the
 *            unaligned raw sequences in <msa>, the way most H3 traces
 *            are supposed to be. In some cases (such as model
 *            construction from an MSA) it is convenient to reference
 *            residues in the MSA cooordinate system directly; setting
 *            <optflags = p7_MSA_COORDS> makes the traces come out
 *            with <i=1..alen> coords for residues.
 *            
 *            The reason these are called "faux" tracebacks is the
 *            following, and it's very important. An input MSA did not
 *            necessarily come from H3. It may imply transitions that
 *            H3 disallows. Thus a returned trace may contain Ik->Dk+1
 *            and Dk->Ik transitions. The input MSA may also imply a
 *            local alignment with no consensus Mk entry; in this
 *            case, a returned <tr[idx]> will be <NULL>.
 *
 *            We use faux tracebacks in two ways. One is in model
 *            construction, where we will pass them to TraceCount(),
 *            for parameterization. Trace counting must either ignore
 *            I->D and D->I transitions and <NULL> traces, or "fix"
 *            them first; see <TraceDoctor()>. Here the goal is to
 *            only count model transitions we've observed; it seems
 *            fair enough to simply ignore events we observe that H3
 *            doesn't model. In effect, this amounts to H3 changing
 *            the input alignment. One extreme example of this is that
 *            in a local alignment path, any residues in "insert"
 *            columns before or after the first consensus residue will
 *            be assigned to N/C state emissions, regardless of their
 *            position in the alignment (that is, the MSA could imply
 *            Ik->Ik+1 transitions -- H3 in effect moves all such
 *            residues out of the homology region.)
 *            
 *            The second way is in merging a given MSA to an
 *            H3-constructed MSA, using <p7_tracealign_*()> functions;
 *            for example, <hmmalign --mapali>. Here the goal is to
 *            recover the exact input alignment, as much as possible.
 *            The <p7_tracealign_*> routines are explicitly designed
 *            to handle illegal D->I/I->D transitions, in order to
 *            enable reconstruction of almost exactly the same input
 *            alignment (the exception being unaligned
 *            insertions). For best results, the caller should not
 *            define any local fragments (because of the
 *            aforementioned issue of possibly pushing residues into
 *            N/C states); glocal alignment paths are more faithfully
 *            reconstructed as the original alignment. 
 *            
 *            If the caller wants to make traces that guarantee
 *            recovery of exactly the input alignment, it should
 *            define all sequences as glocal, and all columns as
 *            consensus (matassign[ai] = 1 for all ai=1..alen).
 *
 *            A related point: input protein alignments may include
 *            '*' characters, which H3 considers to be a "nonresidue"
 *            -- but if the user included them, the user is almost
 *            certainly counting them as residues at least for its
 *            unaligned residue coordinate system.  Therefore we
 *            include '*' characters as emitted residues here, though
 *            they are illegal residue characters for H3. Same issues
 *            apply as above: caller must either doctor these
 *            emissions, or assure that they will be tolerated or
 *            ignored.
 *
 * Args:      msa       - digital alignment
 *            matassign - flag for each alignment column, whether
 *                        it is consensus or not. matassign[1..alen] = 1|0; 
 *                        matassign[0] = 0
 *            optflags  - p7_DEFAULT | p7_MSA_COORDS 
 *            tr        - RETURN: caller provides 0..nseq-1 pointer 
 *                        array for holding returned traces.
 *
 * Returns:   <eslOK> on success, and tr[0..nseq-1] now point to newly
 *            created traces; caller is responsible for freeing these.
 *
 *            These traces are not necessarily "legal", and may not
 *            pass trace validation. Caller must be prepared to deal
 *            with some consequences of mapping an input MSA onto the
 *            restricted H3 profile architecture.  A <tr[idx]> may be
 *            <NULL> (local ali path had no consensus residues for seq
 *            <idx>). A <tr[idx]> may contain D->I and I->D
 *            transitions.  A <tr[idx]> may contain emissions of the
 *            nonresidue '*'.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Notes:     To facilitate future changes in the unsatisfying convention
 *            of flagging local/glocal alignments, we isolate its
 *            decipherment of local vs. glocal: the code that sets
 *            <is_local> can be replaced by whatever the local
 *            vs. glocal MSA annotation convention is.
 *
 *            Because the nonresidue character '*' is treated as a residue
 *            here, you'll see confusing  <if (...XIsResidue || ...XIsNonresidue)> 
 *            OR tests that might look like they would always succeed. Not so;
 *            the other possibilities are gap symbols and missing data symbols.
 *
 * Xref:      J5/17: build.c::fake_tracebacks() becomes p7_trace_FauxFromMSA();
 *                   ability to handle MSA or raw coords added.
 *            J9/42: upgraded to dual-mode local/glocal profile paths.
 */
int
p7_trace_FauxFromMSA(ESL_MSA *msa, int *matassign, int optflags, P7_TRACE **tr)
{		      
  int  idx;			/* counter over seqs in MSA */
  int  k;                       /* position in HMM                 */
  int  ai;                      /* position in alignment columns 1..alen */
  int  i;			/* position in unaligned sequence residues 1..L */
  int  showpos;			/* coord to actually record: ai or i */
  int  lpos_c, rpos_c;		/* leftmost, rightmost consensus column, 1..alen */
  int  lpos, rpos;		/* leftmost, rightmost {MD}k assignment in this seq: glocal: lpos_c,rpos_c; local:  */
  int  is_local;
  int  status = eslOK;

  for (idx = 0; idx < msa->nseq; idx++) tr[idx] = NULL;

  /* Find leftmost and rightmost consensus columns */
  for (lpos_c = 1;         lpos_c <= msa->alen; lpos_c++)  if (matassign[lpos_c]) break;
  for (rpos_c = msa->alen; rpos_c >= 1;         rpos_c--)  if (matassign[rpos_c]) break;
  /* if there were no consensus columns, lpos_c = alen+1 and rpos_c = 0 now; but we specified requirement of at least one consensus column in matassign[] */

  for (idx = 0; idx < msa->nseq; idx++)
    {
      /* Decipher whatever the convention is that tells us whether this is a local (fragment) or glocal sequence alignment */
      is_local = (esl_abc_XIsMissing(msa->abc, msa->ax[idx][1]) || esl_abc_XIsMissing(msa->abc, msa->ax[idx][msa->alen])) ? TRUE : FALSE;
      
      /* Find the first/last column for model entry/exit on this seq. 
       * For glocal, that's always simply the first/last consensus column.
       * For local, it's the first/last col with a residue in it (enter/exit on Mk; in construction, Ik is also possible as special case)
       */
      if (is_local)
	{
	  for (lpos = lpos_c; lpos <= rpos_c; lpos++) 
	    if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][lpos]) || esl_abc_XIsNonresidue(msa->abc, msa->ax[idx][lpos])) break; /* "nonresidue" is the '*' character */
	  for (rpos = rpos_c; rpos >= lpos_c; rpos--)
	    if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][rpos]) || esl_abc_XIsNonresidue(msa->abc, msa->ax[idx][rpos])) break;
	}
      else { lpos = lpos_c; rpos = rpos_c; }

      /* It's possible in a local path that there were no residues in
       * consensus columns, in which case now lpos > rpos (lpos = rpos_c+1,
       * rpos = lpos_c-1) and logic below will fail. Moreover,
       * without any Mk to enter on, we're headed for a B->L->E
       * empty path that isn't legal. So, leave such a trace NULL.
       */
      if (lpos > rpos) continue;

      if ((tr[idx] = p7_trace_Create())                      == NULL) goto ERROR; 
      if ((status  = p7_trace_Append(tr[idx], p7T_S, 0, 0)) != eslOK) goto ERROR;
      if ((status  = p7_trace_Append(tr[idx], p7T_N, 0, 0)) != eslOK) goto ERROR;

      k = 0;
      i = 1; 
      for (ai = 1; ai <= msa->alen; ai++)
	{
	  showpos = (optflags & p7_MSA_COORDS) ? ai : i;
	  if (matassign[ai]) k++;

	  if (ai == lpos) {	/* left edge of homology: model entry */
	    if ((status  = p7_trace_Append(tr[idx], p7T_B, 0, 0)) != eslOK) goto ERROR;
	    if (is_local) { if ((status  = p7_trace_Append(tr[idx], p7T_L, 0, 0)) != eslOK) goto ERROR; }
	    else          { if ((status  = p7_trace_Append(tr[idx], p7T_G, 0, 0)) != eslOK) goto ERROR; }
	  }

	  if  (esl_abc_XIsResidue(msa->abc, msa->ax[idx][ai]) || esl_abc_XIsNonresidue(msa->abc, msa->ax[idx][ai])) 
	    {
	      if      (ai < lpos)     { if ((status  = p7_trace_Append(tr[idx], p7T_N,  0, showpos)) != eslOK) goto ERROR; }
	      else if (ai > rpos)     { if ((status  = p7_trace_Append(tr[idx], p7T_C,  0, showpos)) != eslOK) goto ERROR; }
	      else if (is_local) {
		if (matassign[ai])    { if ((status  = p7_trace_Append(tr[idx], p7T_ML, k, showpos)) != eslOK) goto ERROR; }
		else                  { if ((status  = p7_trace_Append(tr[idx], p7T_IL, k, showpos)) != eslOK) goto ERROR; }
	      } else {
		if (matassign[ai])    { if ((status  = p7_trace_Append(tr[idx], p7T_MG, k, showpos)) != eslOK) goto ERROR; }
		else                  { if ((status  = p7_trace_Append(tr[idx], p7T_IG, k, showpos)) != eslOK) goto ERROR; }
	      }
	      i++;
	    }
	  else if (matassign[ai]) /* delete, or nothing. gaps or ~ in nonconsensus columns. */
	    {
	      if      (is_local && ai > lpos && ai < rpos) { if ((status  = p7_trace_Append(tr[idx], p7T_DL, k, 0)) != eslOK) goto ERROR; }
	      else if (! is_local)                         { if ((status  = p7_trace_Append(tr[idx], p7T_DG, k, 0)) != eslOK) goto ERROR; }
	    }

	  if (ai == rpos) { 	/* right edge of consensus homology region */
	    if ((status = p7_trace_Append(tr[idx], p7T_E, 0, 0)) != eslOK) goto ERROR; 
	    if ((status = p7_trace_Append(tr[idx], p7T_C, 0, 0)) != eslOK) goto ERROR;
	  }
	}
      if ((status = p7_trace_Append(tr[idx], p7T_T, 0, 0)) != eslOK) goto ERROR;

      tr[idx]->M = k;
      tr[idx]->L = (optflags & p7_MSA_COORDS) ? msa->alen : i-1;
    }
  return eslOK;

 ERROR:
  for (idx = 0; idx < msa->nseq; idx++) { p7_trace_Destroy(tr[idx]); tr[idx] = NULL; }
  return status; 
}



/* Function: p7_trace_Doctor()
 * Synopsis: Hack a trace to assure it is valid. 
 * 
 * Purpose:  Plan 7 disallows D->I and I->D "chatter" transitions.
 *           However, these transitions will be implied by many
 *           alignments. trace_doctor() arbitrarily collapses I->D or
 *           D->I into a single M position in the trace.
 *  
 *           Once <tr> has been doctored, it will pass validation
 *           by <p7_trace_Validate()>.
 *           
 *           trace_doctor does not examine any scores when it does
 *           this. In ambiguous situations (D->I->D) the symbol
 *           will be pulled arbitrarily to the left, regardless
 *           of whether that's the best column to put it in or not.
 *           
 * Args:     tr      - trace to doctor (may be NULL)
 *           opt_ndi - optRETURN: number of DI transitions doctored
 *           opt_nid - optRETURN: number of ID transitions doctored
 * 
 * Return:   <eslOK> on success, and the trace <tr> is modified.
 */               
int
p7_trace_Doctor(P7_TRACE *tr, int *opt_ndi, int *opt_nid)
{
  int opos;			/* position in old trace             */
  int npos;			/* position in new trace (<= opos)   */
  int ndi = 0;			/* number of DI transitions doctored */
  int nid = 0;			/* number of ID transitions doctored */

  /* overwrite the trace from left to right */
  if (tr) 
    {
      opos = npos = 0;
      while (opos < tr->N) {
	/* fix implied D->I transitions; D transforms to M, I pulled in */
	if (tr->st[opos] == p7T_DL && tr->st[opos+1] == p7T_IL) 
	  {
	    tr->st[npos] = p7T_ML;
	    tr->k[npos]  = tr->k[opos];     /* D transforms to M      */
	    tr->i[npos]  = tr->i[opos+1];   /* insert char moves back */
	    opos += 2;
	    npos += 1;
	    ndi++;
	  }
	/* ditto for glocal */
	else if (tr->st[opos] == p7T_DG && tr->st[opos+1] == p7T_IG) 
	  {
	    tr->st[npos] = p7T_MG;
	    tr->k[npos]  = tr->k[opos];     /* D transforms to M      */
	    tr->i[npos]  = tr->i[opos+1];   /* insert char moves back */
	    opos += 2;
	    npos += 1;
	    ndi++;
	  }
	/* fix implied I->D transitions; D transforms to M, I is pushed in */
	else if (tr->st[opos]== p7T_IL && tr->st[opos+1]== p7T_DL) 
	  {
	    tr->st[npos] = p7T_ML;
	    tr->k[npos]  = tr->k[opos+1];    /* D transforms to M    */
	    tr->i[npos]  = tr->i[opos];      /* insert char moves up */
	    opos += 2;
	    npos += 1;
	    nid++; 
	  } 
	/* ditto for glocal */
	else if (tr->st[opos]== p7T_IG && tr->st[opos+1]== p7T_DG) 
	  {
	    tr->st[npos] = p7T_MG;
	    tr->k[npos]  = tr->k[opos+1];    /* D transforms to M    */
	    tr->i[npos]  = tr->i[opos];      /* insert char moves up */
	    opos += 2;
	    npos += 1;
	    nid++; 
	  } 
	/* everything else is just copied */
	else {
	  tr->st[npos] = tr->st[opos];
	  tr->k[npos]  = tr->k[opos];
	  tr->i[npos]  = tr->i[opos];
	  opos++;
	  npos++;
	}
      }
      tr->N = npos;
    }

  if (opt_ndi != NULL) *opt_ndi = ndi;
  if (opt_nid != NULL) *opt_nid = nid;
  return eslOK;
}
/*-------------- end, faux traces from MSAs ---------------------*/


/*****************************************************************
 * 7. Counting traces into new HMMs.
 *****************************************************************/

/* Function: p7_trace_Count()
 * 
 * Purpose:  Count a traceback into a count-based core HMM structure.
 *           (Usually as part of a model parameter re-estimation.)
 *           
 *           Tracebacks are relative to the profile model, not
 *           the simpler core HMM. Therefore we have to do some
 *           interpretation; for example, ignoring transitions that
 *           aren't parameterized by observed counts, such as
 *           L->Mk entries and Mk->E local exits.
 * 
 *           Because tracebacks might (and often do) come from an
 *           input alignment via <p7_trace_FauxFromMSA()>, we have to
 *           tolerate some "illegal" special cases: I->D and D->I
 *           transitions, "emissions" of the * nonresidue character,
 *           and NULL traces (empty local alignments). We simply
 *           ignore all these.
 *
 * Args:     hmm   - counts-based HMM to count <tr> into
 *           tr    - alignment of seq to HMM
 *           dsq   - digitized sequence that traceback aligns to the HMM (1..L)
 *                   (or can be an ax, aligned digital seq)
 *           wt    - weight on this sequence
 *           
 * Return:   <eslOK> on success.
 *           Weighted count events are accumulated in hmm's mat[][], ins[][],
 *           t[][] fields: the core probability model.
 *           
 * Throws:   <eslEINVAL> if something's corrupt in the trace; effect on hmm
 *           counts is undefined, because it may abort at any point in the trace.
 */
int
p7_trace_Count(P7_HMM *hmm, ESL_DSQ *dsq, float wt, P7_TRACE *tr)
{
  int z;                        /* position in tr         */
  int i;			/* symbol position in seq */
  int st,st2;     		/* state type (cur, nxt)  */
  int k;	        	/* node index (cur, nxt)  */
  
  if (tr == NULL) return eslOK;	/* tolerate empty local alignments from FauxFromMSA() */

  /* Collect emission counts  */
  for (z = 0; z < tr->N-1; z++) 
    {
      /* pull some info into tmp vars for notational clarity below. */
      st  = tr->st[z]; 
      k   = tr->k[z]; 
      i   = tr->i[z];
      
      if (! i) continue;				     /* skip nonemitters */
      if (esl_abc_XIsNonresidue(hmm->abc, dsq[i])) continue; /* ignore '*' symbols in undoctored traces */

      if      (st == p7T_MG || st == p7T_ML) esl_abc_FCount(hmm->abc, hmm->mat[k], dsq[i], wt);
      else if (st == p7T_IG || st == p7T_IL) esl_abc_FCount(hmm->abc, hmm->ins[k], dsq[i], wt);      
    }

  /* Collect transition counts */
  for (z = 0; z < tr->N-1; z++) 
    {
      /* pull some info into tmp vars for notational clarity below. */
      st  = tr->st[z]; 
      st2 = tr->st[z+1];
      k   = tr->k[z]; 

      /* Ignore stuff that could be in a profile trace, esp. one from FauxFromMSA(), but isn't in a core HMM. */
      if (st == p7T_S || st == p7T_N || st == p7T_B || st == p7T_L ||
	  st == p7T_E || st == p7T_J || st == p7T_C)    continue; /* only count G,M,D,I glocal transitions */
      if (st == p7T_DL && st2 == p7T_IL)                continue; /* no D->I transitions */
      if (st == p7T_DG && st2 == p7T_IG)                continue; /* no D->I transitions */
      if (st == p7T_IL && st2 == p7T_DL)                continue; /* no I->D transitions */
      if (st == p7T_IG && st2 == p7T_DG)                continue; /* no I->D transitions */

      switch (st) {

      case p7T_G:
	switch (st2) {
	case p7T_MG: hmm->t[0][p7H_MM] += wt; break;
	case p7T_DG: hmm->t[0][p7H_MD] += wt; break;
	default: ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
	break;

      case p7T_MG: 
      case p7T_ML:
	switch (st2) {
	case p7T_MG: case p7T_ML: hmm->t[k][p7H_MM] += wt; break;
	case p7T_IG: case p7T_IL: hmm->t[k][p7H_MI] += wt; break;
	case p7T_DG: case p7T_DL: hmm->t[k][p7H_MD] += wt; break;
	case p7T_E: if (st == p7T_MG) hmm->t[k][p7H_MM] += wt; break; /* only glocal exits are counted */
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
	break;

      case p7T_IG:
      case p7T_IL:
	switch (st2) {
	case p7T_MG: case p7T_ML: hmm->t[k][p7H_IM] += wt; break;
	case p7T_IG: case p7T_IL: hmm->t[k][p7H_II] += wt; break;
	case p7T_E:                                        break; /* IL->E can arise in model construction as special case */
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
	break;

      case p7T_DG:
      case p7T_DL:
	switch (st2) {
	case p7T_MG: case p7T_ML: hmm->t[k][p7H_DM] += wt; break;
	case p7T_DG: case p7T_DL: hmm->t[k][p7H_DD] += wt; break;
	case p7T_E: if (st == p7T_DG) hmm->t[k][p7H_DM] += wt; break; /* only glocal exits are counted */
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
	break;

      default: ESL_EXCEPTION(eslEINVAL, "bad state in trace");
      }
    } /* end loop over trace position */
  return eslOK;
}
/*--------------------- end, trace counting ---------------------*/




/*****************************************************************
 * 8. Unit tests
 *****************************************************************/			 
#ifdef p7TRACE_TESTDRIVE

#include "esl_sq.h"
#include "esl_msafile.h"

#define P7_TRACE_SPOTCHECK(tr, zidx, stidx, kidx, iidx)	\
  (((tr)->st[zidx] == (stidx) && (tr)->k[zidx] == (kidx) && (tr)->i[zidx] == (iidx)) ? 1 : 0)

static void
utest_create_msa(ESL_ALPHABET **ret_abc, ESL_MSA **ret_msa, int **ret_matassign, ESL_SQ ***ret_sq)
{
  char         *msg       = "p7_trace.c:: create_msa failed";
  int          *matassign = NULL;
  ESL_ALPHABET *abc       = esl_alphabet_Create(eslAMINO);
  ESL_MSA      *msa       = esl_msa_CreateFromString("\
# STOCKHOLM 1.0\n\
#=GC RF ..x.xx.xx.x..\n\
seq0    ..A.CD.EF.G..\n\
seq1    ~~~~CD.EF.G~~\n\
seq2    gg-.-DeE-.-gg\n\
seq3    ..A.C-e-Fy-..\n\
seq4    ~~Ac--e-F~~~~\n\
seq5    ~g-w--eEFy~~~\n\
seq6    ~g-w--e--y~~~\n\
seq7    ..A.C*.EF.G..\n\
//\n", eslMSAFILE_STOCKHOLM);
  ESL_SQ      **sq        = NULL;
  int           idx,apos;

  if (msa == NULL)                                                esl_fatal(msg);
  if (esl_msa_Digitize(abc, msa, /*errbuf=*/ NULL )     != eslOK) esl_fatal(msg);
  if (( matassign = malloc(sizeof(int) * (msa->alen+2))) == NULL) esl_fatal("malloc");

  matassign[0] = matassign[msa->alen+1] = 0;
  for (apos = 1; apos <= msa->alen; apos++)
    matassign[apos] = (esl_abc_CIsGap(msa->abc, msa->rf[apos-1])? FALSE : TRUE);

  sq = malloc(sizeof(ESL_SQ *)   * msa->nseq);
  for (idx = 0; idx < msa->nseq; idx++)
    if (esl_sq_FetchFromMSA(msa, idx, &(sq[idx]))  != eslOK) esl_fatal(msg);

  *ret_abc       = abc;
  *ret_msa       = msa;
  *ret_matassign = matassign;
  *ret_sq        = sq;
}

static void
utest_FauxFromMSA_seqcoords(ESL_MSA *msa, int *matassign, P7_TRACE ***ret_tr)
{
  char      *msg       = "p7_trace.c:: FauxFromMSA, seqcoords unit test failed";
  P7_TRACE **tr        = NULL;
  int        optflags  = p7_DEFAULT;


  if (( tr     = malloc(sizeof(P7_TRACE *) * msa->nseq))  == NULL)  esl_fatal("malloc");
  if ( p7_trace_FauxFromMSA(msa, matassign, optflags, tr) != eslOK) esl_fatal(msg);
  
  /* Spotchecks of the trickiest bits of those traces. They'll also be Validate()'ed later. */
  /* seq0 = a glocal alignment.  */
  if (! P7_TRACE_SPOTCHECK(tr[0], 3, p7T_G,  0, 0))     esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[0], 9, p7T_MG, 6, 6))     esl_fatal(msg);
  if (tr[0]->N != 13 || tr[0]->M != 6 || tr[0]->L != 6) esl_fatal(msg);

  /* seq1 = a local alignment. */
  if (! P7_TRACE_SPOTCHECK(tr[1], 3, p7T_L,  0, 0))     esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[1], 4, p7T_ML, 2, 1))     esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[1], 8, p7T_ML, 6, 5))     esl_fatal(msg);
  if (tr[1]->N != 12 || tr[1]->M != 6 || tr[1]->L != 5) esl_fatal(msg);

  /* seq2 = one reason ~ convention is bad: local ali with flush N,C can't be represented */
  if (! P7_TRACE_SPOTCHECK(tr[2],  3, p7T_N,  0, 2))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[2],  5, p7T_G,  0, 0))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[2],  6, p7T_DG, 1, 0))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[2],  9, p7T_IG, 3, 4))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[2], 12, p7T_DG, 6, 0))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[2], 16, p7T_C,  0, 7))    esl_fatal(msg);
  if (tr[2]->N != 18 || tr[2]->M != 6 || tr[2]->L != 7) esl_fatal(msg);

  /* seq3 = input MSA implies D->I, I->D transitions, illegal in H3 */
  if (! P7_TRACE_SPOTCHECK(tr[3],  6, p7T_DG, 3, 0))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[3],  7, p7T_IG, 3, 3))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[3],  8, p7T_DG, 4, 0))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[3], 10, p7T_IG, 5, 5))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[3], 11, p7T_DG, 6, 0))    esl_fatal(msg);
  if (tr[3]->N != 15 || tr[3]->M != 6 || tr[3]->L != 5) esl_fatal(msg);

  /* seq4 = ... D->I, I->D transitions in a local alignment path */
  if (! P7_TRACE_SPOTCHECK(tr[4],  5, p7T_IL, 1, 2))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[4],  6, p7T_DL, 2, 0))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[4],  7, p7T_DL, 3, 0))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[4],  8, p7T_IL, 3, 3))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[4],  9, p7T_DL, 4, 0))    esl_fatal(msg);
  if (tr[4]->N != 14 || tr[4]->M != 6 || tr[4]->L != 4) esl_fatal(msg);

  /* seq5 = quirk of L->Mk, Mk->E local entry/exit: nonconsensus IL
   * res outside the entry/exit preserve knowledge of what consensus
   * column they follow, rather than getting lumped into NN/CC
   */
  if (! P7_TRACE_SPOTCHECK(tr[5],  2, p7T_N,  0, 1))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[5],  4, p7T_L,  0, 0))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[5],  5, p7T_IL, 1, 2))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[5],  8, p7T_IL, 3, 3))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[5], 11, p7T_IL, 5, 6))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[5], 13, p7T_C,  0, 0))    esl_fatal(msg);
  if (tr[5]->N != 15 || tr[5]->M != 6 || tr[5]->L != 6) esl_fatal(msg);

  /* seq6 = pathological case, local ali can imply no consensus residues; results in a trace that will not Validate() */
  if (! P7_TRACE_SPOTCHECK(tr[6],  2, p7T_N,   0, 1))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[6],  4, p7T_L,   0, 0))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[6],  5, p7T_IL,  1, 2))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[6],  8, p7T_IL,  3, 3))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[6], 11, p7T_IL,  5, 4))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[6], 12, p7T_E,   0, 0))    esl_fatal(msg);
  if (tr[6]->N != 15 || tr[6]->M != 6 || tr[6]->L != 4)  esl_fatal(msg);

  /* seq7 = nonresidue '*' symbol is accepted in input ali as if it's a residue */
  if (! P7_TRACE_SPOTCHECK(tr[7],  6, p7T_MG,  3, 3))   esl_fatal(msg);
  if (tr[7]->N != 13 || tr[7]->M != 6 || tr[7]->L != 6) esl_fatal(msg);

  *ret_tr = tr;
}		      

/* same test as above, but now in MSA coords 1..alen, not unaligned seq coords */
static void
utest_FauxFromMSA_msacoords(ESL_MSA *msa, int *matassign, P7_TRACE ***ret_tr)
{
  char      *msg       = "p7_trace.c:: FauxFromMSA, MSA coords unit test failed";
  P7_TRACE **tr        = NULL;
  int        optflags  = p7_MSA_COORDS;

  if (( tr     = malloc(sizeof(P7_TRACE *) * msa->nseq))  == NULL)  esl_fatal("malloc");
  if ( p7_trace_FauxFromMSA(msa, matassign, optflags, tr) != eslOK) esl_fatal(msg);
  
  if (! P7_TRACE_SPOTCHECK(tr[0], 9, p7T_MG, 6, 11))     esl_fatal(msg);
  if (tr[0]->N != 13 || tr[0]->M != 6 || tr[0]->L != 13) esl_fatal(msg);

  if (! P7_TRACE_SPOTCHECK(tr[1], 4, p7T_ML, 2, 5))      esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[1], 8, p7T_ML, 6, 11))     esl_fatal(msg);
  if (tr[1]->N != 12 || tr[1]->M != 6 || tr[1]->L != 13) esl_fatal(msg);

  if (! P7_TRACE_SPOTCHECK(tr[2],  3, p7T_N,  0, 2))     esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[2],  9, p7T_IG, 3, 7))     esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[2], 16, p7T_C,  0, 13))    esl_fatal(msg);
  if (tr[2]->N != 18 || tr[2]->M != 6 || tr[2]->L != 13) esl_fatal(msg);

  if (! P7_TRACE_SPOTCHECK(tr[3],  7, p7T_IG, 3, 7))     esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[3], 10, p7T_IG, 5, 10))    esl_fatal(msg);
  if (tr[3]->N != 15 || tr[3]->M != 6 || tr[3]->L != 13) esl_fatal(msg);

  if (! P7_TRACE_SPOTCHECK(tr[4],  5, p7T_IL, 1, 4))     esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[4],  8, p7T_IL, 3, 7))     esl_fatal(msg);
  if (tr[4]->N != 14 || tr[4]->M != 6 || tr[4]->L != 13) esl_fatal(msg);

  if (! P7_TRACE_SPOTCHECK(tr[5],  2, p7T_N,  0, 2))     esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[5],  5, p7T_IL, 1, 4))     esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[5],  8, p7T_IL, 3, 7))     esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[5], 11, p7T_IL, 5, 10))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[5], 13, p7T_C,  0, 0))     esl_fatal(msg);
  if (tr[5]->N != 15 || tr[5]->M != 6 || tr[5]->L != 13) esl_fatal(msg);

  if (! P7_TRACE_SPOTCHECK(tr[6],  2, p7T_N,   0, 2))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[6],  4, p7T_L,   0, 0))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[6],  5, p7T_IL,  1, 4))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[6],  8, p7T_IL,  3, 7))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[6], 11, p7T_IL,  5, 10))   esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[6], 12, p7T_E,   0, 0))    esl_fatal(msg);
  if (tr[6]->N != 15 || tr[6]->M != 6 || tr[6]->L != 13) esl_fatal(msg);

  if (! P7_TRACE_SPOTCHECK(tr[7],  6, p7T_MG,  3, 6))    esl_fatal(msg);
  if (tr[7]->N != 13 || tr[7]->M != 6 || tr[7]->L != 13) esl_fatal(msg);

  *ret_tr = tr;
}		      


static void
utest_Doctor(ESL_SQ **sq, P7_TRACE **tr)
{
  char *msg          = "p7_trace.c:: Doctor unit test failed";
  int   preresult[8] = { eslOK, eslOK, eslOK, eslFAIL, eslFAIL, eslFAIL, eslFAIL, eslOK };
  int   ndi, nid;
  int   idx;
  /* traces 3,4,5 have DI,ID transitions and will fail validation until they're doctored */
  /* trace 6 will NOT validate, even after doctoring! */
  for (idx = 0; idx < 8; idx++)
    {
      if (p7_trace_Validate(tr[idx], sq[idx]->abc, sq[idx]->dsq, /*errbuf=*/NULL) != preresult[idx]) esl_fatal(msg);
      if (p7_trace_Doctor(tr[idx], &ndi, &nid) != eslOK) esl_fatal(msg);
      if (preresult[idx] == eslOK   && (ndi != 0 || nid != 0)) esl_fatal(msg);
      if (preresult[idx] == eslFAIL && (ndi == 0 || nid == 0)) esl_fatal(msg);
      if (idx != 6 &&  p7_trace_Validate(tr[idx], sq[idx]->abc, sq[idx]->dsq, /*errbuf=*/NULL) != eslOK) esl_fatal(msg);
    }      

  /* spotcheck the doctoring in 3 and 4 */
  if (! P7_TRACE_SPOTCHECK(tr[3],  6, p7T_MG, 3, 3))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[3],  7, p7T_DG, 4, 0))    esl_fatal(msg);
  if (! P7_TRACE_SPOTCHECK(tr[3],  8, p7T_MG, 5, 4))    esl_fatal(msg);
  if (tr[3]->N != 13 || tr[3]->M != 6 || tr[3]->L != 5) esl_fatal(msg);

  if (! P7_TRACE_SPOTCHECK(tr[4],  5, p7T_ML, 2, 2))    esl_fatal(msg); /* a I->D pushed right, became ML */
  if (! P7_TRACE_SPOTCHECK(tr[4],  6, p7T_ML, 3, 3))    esl_fatal(msg); /* a D->I pulled left, became ML  */
  if (! P7_TRACE_SPOTCHECK(tr[4],  7, p7T_DL, 4, 0))    esl_fatal(msg); /* this D stayed */
  if (tr[4]->N != 12 || tr[4]->M != 6 || tr[4]->L != 4) esl_fatal(msg);
}


static void
utest_check_counts_are_zero(P7_HMM *hmm)
{
 char   *msg = "p7_trace.c:: zero count check failed";
 int k,x,z;

  /* now all emissions should be zero, including the mat[0] state */
  for (k = 0; k <= hmm->M; k++)
    {
      for (x = 0; x < hmm->abc->K; x++)
	{
	  if (hmm->mat[k][x] != 0.0) esl_fatal(msg);
	  if (hmm->ins[k][x] != 0.0) esl_fatal(msg);
	}
      for (z = 0; z < p7H_NTRANSITIONS; z++)
	if (hmm->t[k][z] != 0.0) esl_fatal(msg);
    }
}

static void 
utest_check_counts_undoctored(P7_HMM *hmm)
{
  char   *msg = "p7_trace.c:: count check (undoctored version) failed";
  int     k, sym;
  ESL_DSQ x;
  float   ct;

  /* check nonzero match emissions; zero them as we go */
  k = 1; sym = 'A'; ct = 4.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 2; sym = 'C'; ct = 4.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 3; sym = 'D'; ct = 3.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 4; sym = 'E'; ct = 5.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 5; sym = 'F'; ct = 6.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 6; sym = 'G'; ct = 3.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;

  /* check nonzero insert emissions, zero as we go */
  k = 1; sym = 'C'; ct = 1.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->ins[k][x] != ct) esl_fatal(msg); hmm->ins[k][x] = 0.0;
  k = 1; sym = 'W'; ct = 2.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->ins[k][x] != ct) esl_fatal(msg); hmm->ins[k][x] = 0.0;
  k = 3; sym = 'E'; ct = 5.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->ins[k][x] != ct) esl_fatal(msg); hmm->ins[k][x] = 0.0;
  k = 5; sym = 'Y'; ct = 3.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->ins[k][x] != ct) esl_fatal(msg); hmm->ins[k][x] = 0.0;

  /* check nonzero transitions, zero as we go */
  if (hmm->t[0][p7H_MM] != 3.0) esl_fatal(msg); hmm->t[0][p7H_MM] = 0.0;
  if (hmm->t[0][p7H_MD] != 1.0) esl_fatal(msg); hmm->t[0][p7H_MD] = 0.0;
  if (hmm->t[1][p7H_MM] != 3.0) esl_fatal(msg); hmm->t[1][p7H_MM] = 0.0;
  if (hmm->t[1][p7H_MI] != 1.0) esl_fatal(msg); hmm->t[1][p7H_MI] = 0.0;
  if (hmm->t[1][p7H_DD] != 1.0) esl_fatal(msg); hmm->t[1][p7H_DD] = 0.0;
  if (hmm->t[2][p7H_MM] != 3.0) esl_fatal(msg); hmm->t[2][p7H_MM] = 0.0;
  if (hmm->t[2][p7H_MD] != 1.0) esl_fatal(msg); hmm->t[2][p7H_MD] = 0.0;
  if (hmm->t[2][p7H_DM] != 1.0) esl_fatal(msg); hmm->t[2][p7H_DM] = 0.0;
  if (hmm->t[2][p7H_DD] != 3.0) esl_fatal(msg); hmm->t[2][p7H_DD] = 0.0;
  if (hmm->t[3][p7H_MM] != 3.0) esl_fatal(msg); hmm->t[3][p7H_MM] = 0.0;
  if (hmm->t[3][p7H_MI] != 1.0) esl_fatal(msg); hmm->t[3][p7H_MI] = 0.0;
  if (hmm->t[3][p7H_IM] != 2.0) esl_fatal(msg); hmm->t[3][p7H_IM] = 0.0;
  if (hmm->t[4][p7H_MM] != 4.0) esl_fatal(msg); hmm->t[4][p7H_MM] = 0.0;
  if (hmm->t[4][p7H_MD] != 1.0) esl_fatal(msg); hmm->t[4][p7H_MD] = 0.0;
  if (hmm->t[4][p7H_DM] != 2.0) esl_fatal(msg); hmm->t[4][p7H_DM] = 0.0;
  if (hmm->t[4][p7H_DD] != 1.0) esl_fatal(msg); hmm->t[4][p7H_DD] = 0.0;
  if (hmm->t[5][p7H_MM] != 3.0) esl_fatal(msg); hmm->t[5][p7H_MM] = 0.0;
  if (hmm->t[5][p7H_MI] != 2.0) esl_fatal(msg); hmm->t[5][p7H_MI] = 0.0;
  if (hmm->t[5][p7H_DD] != 1.0) esl_fatal(msg); hmm->t[5][p7H_DD] = 0.0;
  if (hmm->t[6][p7H_MM] != 2.0) esl_fatal(msg); hmm->t[6][p7H_MM] = 0.0;
  if (hmm->t[6][p7H_DM] != 2.0) esl_fatal(msg); hmm->t[6][p7H_DM] = 0.0;

  /* now all counts should be zeroed */
  utest_check_counts_are_zero(hmm);
}

static void 
utest_check_counts_doctored(P7_HMM *hmm)
{
  char   *msg = "p7_trace.c:: count check (doctored version) failed";
  int     k, sym;
  ESL_DSQ x;
  float   ct;

  /* check nonzero match emissions; zero them as we go */
  k = 1; sym = 'A'; ct = 4.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 2; sym = 'C'; ct = 5.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 2; sym = 'W'; ct = 2.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 3; sym = 'D'; ct = 3.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 3; sym = 'E'; ct = 4.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 4; sym = 'E'; ct = 5.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 5; sym = 'F'; ct = 6.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 5; sym = 'Y'; ct = 1.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 6; sym = 'G'; ct = 3.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;
  k = 6; sym = 'Y'; ct = 1.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->mat[k][x] != ct) esl_fatal(msg); hmm->mat[k][x] = 0.0;

  /* check nonzero insert emissions, zero as we go */
  k = 3; sym = 'E'; ct = 1.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->ins[k][x] != ct) esl_fatal(msg); hmm->ins[k][x] = 0.0;
  k = 5; sym = 'Y'; ct = 1.0; x = esl_abc_DigitizeSymbol(hmm->abc, sym);  if (hmm->ins[k][x] != ct) esl_fatal(msg); hmm->ins[k][x] = 0.0;

  /* check nonzero transitions, zero as we go */
  if (hmm->t[0][p7H_MM] != 3.0) esl_fatal(msg); hmm->t[0][p7H_MM] = 0.0;
  if (hmm->t[0][p7H_MD] != 1.0) esl_fatal(msg); hmm->t[0][p7H_MD] = 0.0;
  if (hmm->t[1][p7H_MM] != 4.0) esl_fatal(msg); hmm->t[1][p7H_MM] = 0.0;
  if (hmm->t[1][p7H_DD] != 1.0) esl_fatal(msg); hmm->t[1][p7H_DD] = 0.0;
  if (hmm->t[2][p7H_MM] != 7.0) esl_fatal(msg); hmm->t[2][p7H_MM] = 0.0;
  if (hmm->t[2][p7H_DM] != 1.0) esl_fatal(msg); hmm->t[2][p7H_DM] = 0.0;
  if (hmm->t[3][p7H_MM] != 4.0) esl_fatal(msg); hmm->t[3][p7H_MM] = 0.0;
  if (hmm->t[3][p7H_MI] != 1.0) esl_fatal(msg); hmm->t[3][p7H_MI] = 0.0;
  if (hmm->t[3][p7H_MD] != 3.0) esl_fatal(msg); hmm->t[3][p7H_MD] = 0.0;
  if (hmm->t[3][p7H_IM] != 1.0) esl_fatal(msg); hmm->t[3][p7H_IM] = 0.0;
  if (hmm->t[4][p7H_MM] != 4.0) esl_fatal(msg); hmm->t[4][p7H_MM] = 0.0;
  if (hmm->t[4][p7H_MD] != 1.0) esl_fatal(msg); hmm->t[4][p7H_MD] = 0.0;
  if (hmm->t[4][p7H_DM] != 3.0) esl_fatal(msg); hmm->t[4][p7H_DM] = 0.0;
  if (hmm->t[5][p7H_MM] != 4.0) esl_fatal(msg); hmm->t[5][p7H_MM] = 0.0;
  if (hmm->t[5][p7H_MI] != 1.0) esl_fatal(msg); hmm->t[5][p7H_MI] = 0.0;
  if (hmm->t[5][p7H_DD] != 1.0) esl_fatal(msg); hmm->t[5][p7H_DD] = 0.0;
  if (hmm->t[6][p7H_MM] != 3.0) esl_fatal(msg); hmm->t[6][p7H_MM] = 0.0;
  if (hmm->t[6][p7H_DM] != 1.0) esl_fatal(msg); hmm->t[6][p7H_DM] = 0.0;

  /* now all counts should be zeroed */
  utest_check_counts_are_zero(hmm);
}

/* pass undoctored msa-coord traces into this; they're doctored here */
static void
utest_Count(ESL_MSA *msa, int *matassign, P7_TRACE **trm)
{
  char   *msg = "p7_trace.c:: Count unit test failed";
  P7_HMM *hmm = NULL;
  int     M;
  int     apos,idx;
  
  for (M = 0, apos = 1; apos <= msa->alen; apos++) if (matassign[apos]) M++;
  if ((hmm = p7_hmm_Create(M, msa->abc)) == NULL)  esl_fatal(msg);
  if (p7_hmm_Zero(hmm)                   != eslOK) esl_fatal(msg);

  /* First count and check the undoctored traces */
  for (idx = 0; idx < msa->nseq; idx++) 
    if (p7_trace_Count(hmm, msa->ax[idx], msa->wgt[idx], trm[idx]) != eslOK) esl_fatal(msg);
  utest_check_counts_undoctored(hmm);

  /* Then do it again with doctored traces. */
  for (idx = 0; idx < msa->nseq; idx++) 
    {
      if (p7_trace_Doctor(trm[idx], NULL, NULL)                      != eslOK) esl_fatal(msg);
      if (p7_trace_Validate(trm[idx], msa->abc, msa->ax[idx], NULL)  != eslOK) esl_fatal(msg);
      if (p7_trace_Count(hmm, msa->ax[idx], msa->wgt[idx], trm[idx]) != eslOK) esl_fatal(msg);
    }
  utest_check_counts_doctored(hmm);

  p7_hmm_Destroy(hmm);
}

#if 0
/* convert an MSA to traces; then traces back to MSA; 
 * starting and ending MSA should be the same, provided
 * the msa doesn't have any ambiguously aligned insertions.
 */
static void
utest_faux_tracealign(ESL_MSA *msa, int *matassign, int M)
{
  char      *msg  = "p7_trace.c:: FauxFromMSA unit test failed";
  ESL_MSA   *msa2 = NULL;
  ESL_SQ   **sq   = malloc(sizeof(ESL_SQ)   * msa->nseq);
  P7_TRACE **tr   = malloc(sizeof(P7_TRACE) * msa->nseq);
  int        i;
  int        optflags = p7_DIGITIZE;

  for (i = 0; i < msa->nseq; i++)
    if (esl_sq_FetchFromMSA(msa, i, &(sq[i]))                   != eslOK) esl_fatal(msg);

  if (p7_trace_FauxFromMSA(msa, matassign, p7_MSA_COORDS, tr)   != eslOK) esl_fatal(msg);
  if (p7_tracealign_MSA(msa, tr, M, optflags, &msa2)            != eslOK) esl_fatal(msg);
  if (esl_msa_Compare(msa, msa2)                                != eslOK) esl_fatal(msg);
  esl_msa_Destroy(msa2);
  for (i = 0; i < msa->nseq; i++) p7_trace_Destroy(tr[i]);

  if (p7_trace_FauxFromMSA(msa, matassign, p7_DEFAULT, tr)            != eslOK) esl_fatal(msg);
  if (p7_tracealign_Seqs(sq, tr, msa->nseq, M, optflags, NULL, &msa2) != eslOK) esl_fatal(msg);
  if (esl_msa_Compare(msa, msa2)                                      != eslOK) esl_fatal(msg);

  esl_msa_Destroy(msa2);
  for (i = 0; i < msa->nseq; i++) p7_trace_Destroy(tr[i]);
  for (i = 0; i < msa->nseq; i++) esl_sq_Destroy(sq[i]);
  free(tr);
  free(sq);
  return;
}
#endif

#endif /*p7TRACE_TESTDRIVE*/
/*--------------------- end, unit tests -------------------------*/



/*****************************************************************
 * 9. Test driver
 *****************************************************************/			 
#ifdef p7TRACE_TESTDRIVE
#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for P7_TRACE";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go        = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_ALPHABET   *abc       = NULL;
  ESL_MSA        *msa       = NULL;
  int            *matassign = NULL;
  ESL_SQ        **sq        = NULL;
  P7_TRACE      **trs       = NULL; /* traces in unaligned seq coords */
  P7_TRACE      **trm       = NULL; /* traces in MSA coords */
  int             idx;

  fprintf(stderr, "## %s\n", argv[0]);

  utest_create_msa(&abc, &msa, &matassign, &sq);
  utest_FauxFromMSA_seqcoords(msa, matassign, &trs);
  utest_FauxFromMSA_msacoords(msa, matassign, &trm);
  utest_Doctor(sq, trs);	    /* now trs's are doctored on return */
  utest_Count(msa, matassign, trm); /* trm's are doctored here too      */

  p7_trace_DestroyArray(trs, msa->nseq);
  p7_trace_DestroyArray(trm, msa->nseq);
  for (idx = 0; idx < msa->nseq; idx++) esl_sq_Destroy(sq[idx]);
  free(sq);
  free(matassign);
  esl_msa_Destroy(msa);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);

  fprintf(stderr, "#  status = ok\n");
  return eslOK;
}
#endif /*p7TRACE_TESTDRIVE*/
/*--------------------- end, test driver ------------------------*/


/*****************************************************************
 * 10. Example
 *****************************************************************/
#ifdef p7TRACE_EXAMPLE

/* To run:     ./p7_trace_example foo.sto
 * The foo.sto alignment must have a RF annotation line marking consensus columns.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type         default   env  range   toggles   reqs   incomp   help                                    docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL,   "show brief help on version and usage",        0 },
  { "-d",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL,   "doctor traces",                               0 },
  { "-m",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL,   "use msa coords, not unaligned seq coords",    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msafile> <seqfile>";
static char banner[] = "example driver for P7_TRACE";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go        = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *msafile   = esl_opt_GetArg(go, 1);
  int             infmt     = eslMSAFILE_UNKNOWN;
  ESL_ALPHABET   *abc       = NULL;
  ESLX_MSAFILE   *afp       = NULL;
  ESL_MSA        *msa       = NULL;
  int            *matassign = NULL;
  int             optflags  = (esl_opt_GetBoolean(go, "-m") ? p7_MSA_COORDS : p7_DEFAULT);
  P7_TRACE      **trarr     = NULL;
  int             apos,idx;
  int             status;

  if ( (status = eslx_msafile_Open(&abc, msafile, NULL, infmt, NULL, &afp)) != eslOK)
    eslx_msafile_OpenFailure(afp, status);

  if ( (status = eslx_msafile_Read(afp, &msa)) != eslOK) 
    eslx_msafile_ReadFailure(afp, status);

  if (! msa->rf)
    esl_fatal("MSA must have #=GC RF consensus column annotation line for this example to work.");

  if (( matassign = malloc(sizeof(int) *        (msa->alen+2))) == NULL) esl_fatal("malloc failed");
  if (( trarr     = malloc(sizeof(P7_TRACE *) * msa->nseq))     == NULL) esl_fatal("malloc failed");

  matassign[0] = matassign[msa->alen+1] = 0;
  for (apos = 1; apos <= msa->alen; apos++)
    matassign[apos] = (esl_abc_CIsGap(msa->abc, msa->rf[apos-1])? FALSE : TRUE);

  if ( (status = p7_trace_FauxFromMSA(msa, matassign, optflags, trarr)) != eslOK)
    esl_fatal("FauxFromMSA() failed");

  if (esl_opt_GetBoolean(go, "-d")) {
    for (idx = 0; idx < msa->nseq; idx++)
      p7_trace_Doctor(trarr[idx], NULL, NULL);
  }

  for (idx = 0; idx < msa->nseq; idx++)
    if (p7_trace_Validate(trarr[idx], abc, msa->ax[idx], NULL) == eslFAIL)
      printf("warning: trace %d failed validation\n", idx);

  for (idx = 0; idx < msa->nseq; idx++)
    {
      printf("\n### trace %d:\n", idx);
      p7_trace_Dump(stdout, trarr[idx]);
    }

  for (idx = 0; idx < msa->nseq; idx++) p7_trace_Destroy(trarr[idx]);
  free(trarr);
  free(matassign);
  esl_msa_Destroy(msa);
  eslx_msafile_Close(afp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7TRACE_EXAMPLE*/



/************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id: p7_trace.c 3474 2011-01-17 13:25:32Z eddys $
 ************************************************************/

