/* Emitting (sampling) sequences from an HMM, in either core or
 * profile form.
 * 
 *    1. Exported API: sequence emission routines.
 *    2. Private functions.
 *    3. Routines emulating HMMER2, for testing.
 *    4. Stats driver.
 *    5. Copyright and license.
 * 
 * SRE, Tue Jan  9 08:55:53 2007 [Janelia] [The Crystal Method, Vegas]
 * SVN $Id$
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "hmmer.h"

static int sample_endpoints(ESL_RANDOMNESS *r, P7_PROFILE *gm, int *ret_kstart, int *ret_kend);


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
 *            Trace is relative to the core model: it may include
 *            I_0 and I_M states, B->DD->M entry is explicit, and a
 *            0 length generated sequence is possible.
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

  if (sq != NULL) esl_sq_Reuse(sq);    
  if (tr != NULL) {
    if ((status = p7_trace_Reuse(tr))                != eslOK) goto ERROR;
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
      case p7_STM:
	switch (esl_rnd_FChoose(r, hmm->t[k], 3)) {
	case 0:  st = p7_STM; break;
	case 1:  st = p7_STI; break;
	case 2:  st = p7_STD; break;
	default: ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
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
	switch (esl_rnd_FChoose(r, hmm->t[k]+5, 2)) {
	case 0: st = p7_STM; break;
	case 1: st = p7_STD; break;
	default: ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
	}
	break;

      default: ESL_XEXCEPTION(eslECORRUPT, "impossible state reached during emission");
      }

      /* Bump k,i if needed, depending on new state type */
      if (st == p7_STM || st == p7_STD) k++;
      if (st == p7_STM || st == p7_STI) i++;

      /* a transit to M_M+1 is a transit to the E state */
      if (k == hmm->M+1) {
	if   (st == p7_STM) { st = p7_STE; k = 0; }
	else ESL_XEXCEPTION(eslECORRUPT, "failed to reach E state properly");
      }

      /* Sample new residue x if in match or insert */
      if      (st == p7_STM) x = esl_rnd_FChoose(r, hmm->mat[k], hmm->abc->K);
      else if (st == p7_STI) x = esl_rnd_FChoose(r, hmm->ins[k], hmm->abc->K);
      else                   x = eslDSQ_SENTINEL;

      /* Add state to trace */
      if (tr != NULL) {
	if (x == eslDSQ_SENTINEL) {
	  if ((status = p7_trace_Append(tr, st, k, 0)) != eslOK) goto ERROR;
	} else {
	  if ((status = p7_trace_Append(tr, st, k, i)) != eslOK) goto ERROR;
	}
      }

      /* Add x to sequence */
      if (sq != NULL && x != eslDSQ_SENTINEL) 
	if ((status = esl_sq_XAddResidue(sq, x)) != eslOK) goto ERROR;
    }

  /* Terminate the sequence and return */
  if ((status = esl_sq_XAddResidue(sq, eslDSQ_SENTINEL)) != eslOK) goto ERROR;
  return eslOK;

ERROR:
  return status;
}


/* Function:  p7_ProfileEmit()
 * Incept:    SRE, Mon Jan 22 10:23:28 2007 [Janelia]
 *
 * Purpose:   Sample a sequence from the implicit 
 *            probabilistic model of a Plan7 profile <gm>. 
 *            
 *            Optionally return the sequence and/or its trace in <sq>
 *            and <tr>, respectively. Caller has allocated space for
 *            both of these (though they may need to be
 *            reallocated/grown here). Either can be passed as <NULL>
 *            if unneeded. 
 *            
 *            Only the sequence field is set in the <sq>. Caller must
 *            set the name, plus any other fields it wants to set.
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
p7_ProfileEmit(ESL_RANDOMNESS *r, P7_PROFILE *gm, ESL_SQ *sq, P7_TRACE *tr)
{
  char      prv, st;		/* prev, current state type */
  int       k;		        /* position in model nodes 1..M */
  int       i;			/* position in sequence 1..L */
  int       x;			/* sampled residue */
  int       kend;	        /* predestined end node */
  int       status;

  if (sq != NULL) esl_sq_Reuse(sq);    
  if (tr != NULL) {
    if ((status = p7_trace_Reuse(tr))                != eslOK) goto ERROR;
    if ((status = p7_trace_Append(tr, p7_STS, 0, 0)) != eslOK) goto ERROR;
    if ((status = p7_trace_Append(tr, p7_STN, 0, 0)) != eslOK) goto ERROR;
  }
  st    = p7_STN;
  k     = 0;
  i     = 0;
  while (st != p7_STT)
    {
      /* Sample a state transition. After this section, prv and st (prev->current state) are set;
       * k also gets set if we make a B->Mk entry transition.
       */
      prv = st;
      switch (st) {
      case p7_STB:  /* Enter the implicit profile: choose our entry and our predestined exit */
	if ((status = sample_endpoints(r, gm, &k, &kend)) != eslOK) goto ERROR;
	st = p7_STM;		/* must be, because left wing is retracted */
	break;
	
      case p7_STM:
	if (k == kend) st = p7_STE; /* check our preordained fate in the implicit model */
	else {
	  switch (esl_rnd_FChoose(r, gm->hmm->t[k], 3)) {
	  case 0:  st = p7_STM; break;
	  case 1:  st = p7_STI; break;
	  case 2:  st = p7_STD; break;
	  default: ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
	  }
	}
	break;

      case p7_STD:
	if (k == kend) st = p7_STE; 
	else           st = (esl_rnd_FChoose(r, gm->hmm->t[k]+5, 2) == 0) ? p7_STM : p7_STD; 
	break;

      case p7_STI: st = (esl_rnd_FChoose(r, gm->hmm->t[k]+3,2) == 0) ? p7_STM : p7_STI;  break;
      case p7_STN: st = (esl_rnd_FChoose(r, gm->xt[p7_XTN], 2) == 0) ? p7_STB : p7_STN;  break;
      case p7_STE: st = (esl_rnd_FChoose(r, gm->xt[p7_XTE], 2) == 0) ? p7_STC : p7_STJ;  break;
      case p7_STJ: st = (esl_rnd_FChoose(r, gm->xt[p7_XTJ], 2) == 0) ? p7_STB : p7_STJ;  break;
      case p7_STC: st = (esl_rnd_FChoose(r, gm->xt[p7_XTC], 2) == 0) ? p7_STT : p7_STC;  break;
      default:     ESL_XEXCEPTION(eslECORRUPT, "impossible state reached during emission");
      }
     
      /* Based on the transition we just sampled, update k. */
      if      (st == p7_STE)                  k = 0;
      else if (st == p7_STM && prv != p7_STB) k++;    /* be careful about B->Mk, where we already set k */
      else if (st == p7_STD)                  k++;

      /* Based on the transition we just sampled, generate a residue. */
      if      (st == p7_STM)                                              x = esl_rnd_FChoose(r, gm->hmm->mat[k], gm->abc->K);
      else if (st == p7_STI)                                              x = esl_rnd_FChoose(r, gm->hmm->ins[k], gm->abc->K);
      else if ((st == p7_STN || st == p7_STC || st == p7_STJ) && prv==st) x = esl_rnd_FChoose(r, gm->bg->f,       gm->abc->K);
      else    x = eslDSQ_SENTINEL;

      if (x != eslDSQ_SENTINEL) i++;

      /* Add residue (if any) to sequence */
      if (sq != NULL && x != eslDSQ_SENTINEL && (status = esl_sq_XAddResidue(sq, x)) != eslOK) goto ERROR;

      /* Add state to trace; distinguish emitting position (pass i=i) from non (pass i=0) */
      if (tr != NULL) {
	if (x == eslDSQ_SENTINEL) {
	  if ((status = p7_trace_Append(tr, st, k, 0)) != eslOK) goto ERROR;
	} else {
	  if ((status = p7_trace_Append(tr, st, k, i)) != eslOK) goto ERROR;
	}
      }
    }
  /* Terminate the sequence (if we're generating one) */
  if (sq != NULL && (status = esl_sq_XAddResidue(sq, eslDSQ_SENTINEL)) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  return status;
}

/*****************************************************************
 * 2. Private functions.
 *****************************************************************/

/* sample_endpoints()
 * Incept:    SRE, Mon Jan 22 10:43:20 2007 [Janelia]
 *
 * Purpose:   Given a profile <gm> and random number source <r>, sample
 *            a begin transition from the implicit probabilistic profile
 *            model, yielding a sampled start and end node; return these
 *            via <ret_kstart> and <ret_kend>.
 *            
 *            By construction, the entry at node <kstart> is into a
 *            match state, but the exit from node <kend> might turn
 *            out to be from either a match or delete state.
 *            
 *            We assume that exits j are uniformly distributed for a
 *            particular entry point i: $a_{ij} =$ constant $\forall
 *            j$.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            
 * Xref:      STL11/138           
 */
static int
sample_endpoints(ESL_RANDOMNESS *r, P7_PROFILE *gm, int *ret_kstart, int *ret_kend)
{
  float *pstart = NULL;
  int    k;
  int    kstart, kend;
  int    status;

  ESL_ALLOC(pstart, sizeof(float) * (gm->M+1));
  pstart[0] = 0.;
  for (k = 1; k <= gm->M; k++)
    pstart[k] = gm->begin[k] * (gm->M - k + 1);	/* multiply p_ij by the number of exits j */
  kstart = esl_rnd_FChoose(r, pstart, gm->M+1);	/* sample the starting position from that distribution */
  kend   = kstart + esl_rnd_Choose(r, gm->M-kstart+1); /* and the exit uniformly from possible exits for it */

  free(pstart);
  *ret_kstart = kstart;
  *ret_kend   = kend;
  return eslOK;
  
 ERROR:
  if (pstart != NULL) free(pstart);
  *ret_kstart = 0;
  *ret_kend   = 0;
  return status;
}


/*****************************************************************
 * 3. Routines emulating HMMER2, for testing.
 *****************************************************************/

/* In a HMMER2 profile, the DP model is probabilistic (there is no implicit
 * model), with begin[] as a probability distribution and end[k] being treated
 * as a fourth match transition; sample_endpoints (for the H3 style) is
 * not used.
 */
int
p7_H2_ProfileEmit(ESL_RANDOMNESS *r, P7_PROFILE *gm, ESL_SQ *sq, P7_TRACE *tr)
{
  char      prv, st;		/* prev, current state type */
  int       k;		        /* position in model nodes 1..M */
  int       i;			/* position in sequence 1..L */
  int       x;			/* sampled residue */
  int       status;
  float     p[4];

  if (sq != NULL) esl_sq_Reuse(sq);    
  if (tr != NULL) {
    if ((status = p7_trace_Reuse(tr))                != eslOK) goto ERROR;
    if ((status = p7_trace_Append(tr, p7_STS, 0, 0)) != eslOK) goto ERROR;
    if ((status = p7_trace_Append(tr, p7_STN, 0, 0)) != eslOK) goto ERROR;
  }
  st    = p7_STN;
  k     = 0;
  i     = 0;
  while (st != p7_STT)
    {
      /* Sample a state transition. After this section, prv and st (prev->current state) are set;
       * k also gets set if we make a B->Mk entry transition.
       */
      prv = st;
      switch (st) {
      case p7_STB:  /* Enter the implicit profile: choose our entry and our predestined exit */
	k  = 1 + esl_rnd_FChoose(r, gm->begin+1, gm->M);
	st = p7_STM;		/* must be, because left wing is retracted */
	break;
	
      case p7_STM:
	if (k == gm->M) st = p7_STE;
	else {
	  esl_vec_FCopy(gm->hmm->t[k], 3, p);
	  p[3] = gm->end[k];
	  switch (esl_rnd_FChoose(r, p, 4)) {
	  case 0:  st = p7_STM; break;
	  case 1:  st = p7_STI; break;
	  case 2:  st = p7_STD; break;
	  case 3:  st = p7_STE; break;
	  default: ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
	  }
	}
	break;

      case p7_STD:
	if (k == gm->M) st = p7_STE;
	else            st = (esl_rnd_FChoose(r, gm->hmm->t[k]+5,2) == 0) ? p7_STM : p7_STD;  
	break;

      case p7_STI: st = (esl_rnd_FChoose(r, gm->hmm->t[k]+3,2) == 0) ? p7_STM : p7_STI;  break;
      case p7_STN: st = (esl_rnd_FChoose(r, gm->xt[p7_XTN], 2) == 0) ? p7_STB : p7_STN;  break;
      case p7_STE: st = (esl_rnd_FChoose(r, gm->xt[p7_XTE], 2) == 0) ? p7_STC : p7_STJ;  break;
      case p7_STJ: st = (esl_rnd_FChoose(r, gm->xt[p7_XTJ], 2) == 0) ? p7_STB : p7_STJ;  break;
      case p7_STC: st = (esl_rnd_FChoose(r, gm->xt[p7_XTC], 2) == 0) ? p7_STT : p7_STC;  break;
      default:     ESL_XEXCEPTION(eslECORRUPT, "impossible state reached during emission");
      }
     
      /* Based on the transition we just sampled, update k. */
      if      (st == p7_STE)                  k = 0;
      else if (st == p7_STM && prv != p7_STB) k++;    /* be careful about B->Mk, where we already set k */
      else if (st == p7_STD)                  k++;

      /* Based on the transition we just sampled, generate a residue. */
      if      (st == p7_STM)                                              x = esl_rnd_FChoose(r, gm->hmm->mat[k], gm->abc->K);
      else if (st == p7_STI)                                              x = esl_rnd_FChoose(r, gm->hmm->ins[k], gm->abc->K);
      else if ((st == p7_STN || st == p7_STC || st == p7_STJ) && prv==st) x = esl_rnd_FChoose(r, gm->bg->f,       gm->abc->K);
      else    x = eslDSQ_SENTINEL;

      if (x != eslDSQ_SENTINEL) i++;

      /* Add residue (if any) to sequence */
      if (sq != NULL && x != eslDSQ_SENTINEL && (status = esl_sq_XAddResidue(sq, x)) != eslOK) goto ERROR;

      /* Add state to trace; distinguish emitting position (pass i=i) from non (pass i=0) */
      if (tr != NULL) {
	if (x == eslDSQ_SENTINEL) {
	  if ((status = p7_trace_Append(tr, st, k, 0)) != eslOK) goto ERROR;
	} else {
	  if ((status = p7_trace_Append(tr, st, k, i)) != eslOK) goto ERROR;
	}
      }
    }
  /* Terminate the sequence (if we're generating one) */
  if (sq != NULL && (status = esl_sq_XAddResidue(sq, eslDSQ_SENTINEL)) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  return status;
}


/*****************************************************************
 * 4. Stats driver
 *****************************************************************/

/* A small driver providing a testbed for sequence-emission related development testing.
 * 
 * gcc -g -Wall -o stats -L. -I. -L../easel -I../easel -Dp7EMIT_STATS emit.c -lhmmer -leasel -lm
 */
#ifdef p7EMIT_STATS
#include "p7_config.h"
#include <stdio.h>
#include <easel.h>
#include <esl_alphabet.h>
#include <esl_random.h>
#include <esl_sqio.h>
#include "hmmer.h"

int
main(int argc, char **argv)
{
  char            *hmmfile = argv[1];  /* name of HMM file to read one HMM from   */
  ESL_ALPHABET    *abc     = NULL;     /* sequence alphabet                       */
  ESL_RANDOMNESS  *r       = NULL;     /* source of randomness                    */
  P7_HMMFILE      *hfp     = NULL;     /* open hmmfile                            */
  P7_HMM          *hmm     = NULL;     /* HMM to emit from                        */
  P7_PROFILE      *gm      = NULL;     /* profile HMM (scores)                    */
  P7_BG           *bg      = NULL;     /* null model                              */
  P7_TRACE        *tr      = NULL;     /* sampled trace                           */
  ESL_SQ          *sq      = NULL;     /* sampled digital sequence                */
  int              n       = 1000;
  int              counts[p7_NSTATETYPES];
  int              i;
  int              sc;
  int              nullsc;
  double           bitscore;

  r = esl_randomness_CreateTimeseeded();
  p7_trace_Create(256, &tr);
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("failed to open %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read HMM");
  sq = esl_sq_CreateDigital(abc);
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);

  p7_ProfileConfig(hmm, bg, gm, p7_LOCAL);

  for (i = 0; i < n; i++) 
    {
      p7_ProfileEmit(r, gm, sq, tr);
      p7_trace_StateUseCounts(tr, counts);

      p7_ReconfigLength(gm, sq->n);
      p7_bg_SetLength(bg, sq->n);
      p7_trace_Score(tr, sq->dsq, gm, &sc);
      p7_bg_NullOne (bg, sq->dsq, sq->n, &nullsc);
      bitscore = p7_SILO2Bitscore(sc - nullsc);

      printf("%d  %8.2f\n",
	     counts[p7_STM] + (counts[p7_STI] + counts[p7_STD])/2,
	     bitscore);
    }

  p7_profile_Destroy(gm);
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  p7_hmmfile_Close(hfp);
  p7_hmm_Destroy(hmm);
  return eslOK;
} 


#endif /*p7EMIT_STATS*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
