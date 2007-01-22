/* Emitting (sampling) sequences from an HMM, in either core or
 * profile form.
 * 
 *    1. Exported API: sequence emission routines.
 *    2. Copyright and license.
 * 
 * SRE, Tue Jan  9 08:55:53 2007 [Janelia] [The Crystal Method, Vegas]
 * SVN $Id$
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"
#include "esl_sqio.h"

#include "p7_hmm.h"
#include "p7_trace.h"

/*****************************************************************
 * 1. Exported API: sequence emission routines.
 *****************************************************************/

/* Function:  p7_CoreEmit()
 * Incept:    SRE, Tue Jan  9 10:20:51 2007 [Janelia]
 *
 * Purpose:   Generate (sample) a sequence from a core profile HMM <hmm>.
 *            
 *            Optionally return the sequence and/or its trace in <sq>
 *            and <tr>, respectively, which the caller has
 *            allocated. Having the caller provide these reusable
 *            objects allows re-use of both <sq> and <tr> in repeated
 *            calls, saving malloc/free wastage. Either can be passed
 *            as <NULL> if it isn't needed.
 *            
 *            This does not set any fields in the <sq> except for the
 *            sequence itself. Caller must set the name, and any other
 *            annotation it wants to add.
 *
 * Note:      Traces are always relative to the wing-retracted profile 
 *            form model, but core models can traverse D_1 and D_M states;
 *            in fact, core model can generate an empty B->DDDD->E 
 *            path that the profile form model can't accommodate. 
 *            So, we do some special stuff to retract wings of the trace,
 *            and to catch the empty case.
 *            
 * Args:      r     -  source of randomness
 *            hmm   -  core HMM to generate from
 *            sq    -  opt: digital sequence sampled (or NULL)
 *            tr    -  opt: trace sampled            (or NULL)
 *
 * Returns:   <eslOK> on success; 
 *            optionally return the digital sequence through <ret_sq>,
 *            and optionally return its trace in <ret_tr>.
 *
 * Throws:    <eslECORRUPT> if emission gets us into an illegal state, 
 *            probably indicating that a probability that should have
 *            been zero wasn't. 
 *
 *            Throws <eslEMEM> on a reallocation error.
 * 
 *            In these cases, the contents of <sq> and <tr> may be
 *            corrupted. Caller should not trust their data, but may
 *            safely reuse them.
 *
 * Xref:      STL11/124.
 */
int
p7_CoreEmit(ESL_RANDOMNESS *r, P7_HMM *hmm, ESL_SQ *sq, P7_TRACE *tr)
{
  char      st;			/* state type */
  int       k;			/* position in model nodes 1..M */
  int       i;			/* position in sequence 1..L */
  int       x;			/* sampled residue */
  int       status;

  /* The entire routine is wrapped in a do/while, in order to reject
   * one pathological case of an L=0 sampled sequence.
   */
  do {
    if (sq != NULL) esl_sq_Reuse(sq);    
    if (tr != NULL) {
      if ((status = p7_trace_Reuse(tr))                != eslOK) goto ERROR;
      if ((status = p7_trace_Append(tr, p7_STS, 0, 0)) != eslOK) goto ERROR;
      if ((status = p7_trace_Append(tr, p7_STN, 0, 0)) != eslOK) goto ERROR;
      if ((status = p7_trace_Append(tr, p7_STB, 0, 0)) != eslOK) goto ERROR;
    }
    st    = p7_STB;
    k     = 0;
    i     = 0;
    while (st != p7_STE)
      {
	/* Sample next state type, given current state type (and current k) */
	switch (st) {
	case p7_STB:
	  switch (esl_rnd_FChoose(r, hmm->t[0], 3)) {
	  case 0:  st = p7_STM; break;
	  case 1:  ESL_XEXCEPTION(eslECORRUPT, "Whoa, shouldn't be able to sample a B->I path."); 
	  case 2:  st = p7_STD; break;
          default: ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
	  }
	  break;

	case p7_STM:
	  if (k == hmm->M) st = p7_STE;
	  else {
	    switch (esl_rnd_FChoose(r, hmm->t[k], 3)) {
	    case 0:  st = p7_STM; break;
	    case 1:  st = p7_STI; break;
	    case 2:  st = p7_STD; break;
	    default: ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
	    }
	  }
	  break;

	case p7_STI:
	  switch (esl_rnd_FChoose(r, hmm->t[k]+3, 2)) {
          case 0: st = p7_STM; break;
	  case 1: st = p7_STI; break;
	  default: ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
	  }
	  break;

	case p7_STD:
	  if (k == hmm->M) st = p7_STE;
	  else {
	    switch (esl_rnd_FChoose(r, hmm->t[k]+5, 2)) {
	    case 0: st = p7_STM; break;
	    case 1: st = p7_STD; break;
	    default: ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
	    }
	  }
	  break;

	default: ESL_XEXCEPTION(eslECORRUPT, "impossible state reached during emission");
	}

	/* Bump k,i if needed, depending on new state type */
	if (st == p7_STM || st == p7_STD) k++;
	if (st == p7_STM || st == p7_STI) i++;

	/* Sample new residue x if in match or insert */
	if      (st == p7_STM) x = esl_rnd_FChoose(r, hmm->mat[k], hmm->abc->K);
	else if (st == p7_STI) x = esl_rnd_FChoose(r, hmm->ins[k], hmm->abc->K);
	else    x = eslDSQ_SENTINEL;

	/* Add state to trace */
	if (tr != NULL) {
	  if (st == p7_STE)		/* Handle right wing retraction: rewind over D's in a Mk->DDD->E path */
	    while (tr->st[tr->N-1] == p7_STD) tr->N--; /* unwarranted chumminess with trace structure     */
	  if (! (tr->st[tr->N-1] == p7_STB && st == p7_STD))   /* the case handles left wing retraction */
	    if ((status = p7_trace_Append(tr, st, k, i)) != eslOK) goto ERROR;
	}
	/* Note we might be B->E here! */      

	/* Add x to sequence */
	if (sq != NULL && x != eslDSQ_SENTINEL) 
	  if ((status = esl_sq_XAddResidue(sq, x)) != eslOK) goto ERROR;
      }

    /* Last state reached was E; now finish the trace: */
    if (tr != NULL) {
      if ((status = p7_trace_Append(tr, p7_STC, 0, 0)) != eslOK) goto ERROR;
      if ((status = p7_trace_Append(tr, p7_STT, 0, 0)) != eslOK) goto ERROR;
    }
  } while (i == 0);

  /* Terminate the sequence and return */
  if ((status = esl_sq_XAddResidue(sq, eslDSQ_SENTINEL)) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  return status;
}



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
