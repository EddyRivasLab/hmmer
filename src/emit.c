/************************************************************
 * @LICENSE@
 ************************************************************/

/* emit.c
 * SRE, Sun Mar  8 12:26:58 1998
 * CVS $Id$
 * 
 * Generation of sequences/traces from an HMM.
 */

#include "config.h"
#include "squidconf.h"

#include <ctype.h>

#include "structs.h"
#include "funcs.h"
#include "squid.h"



/* Function: EmitSequence()
 * Date:     SRE, Sun Mar  8 12:28:03 1998 [St. Louis]
 *
 * Purpose:  Given a model, sample a sequence and/or traceback.
 *
 * Args:     hmm     - the model
 *           ret_dsq - RETURN: generated digitized sequence (pass NULL if unwanted)
 *           ret_L   - RETURN: length of generated sequence 
 *           ret_tr  - RETURN: generated trace (pass NULL if unwanted)
 *
 * Returns:  void
 */
void
EmitSequence(struct plan7_s *hmm, char **ret_dsq, int *ret_L, struct p7trace_s **ret_tr)
{
  struct p7trace_s *tr;
  char  type;               	/* current state type */
  int   k;			/* current node index */
  char *dsq;                    /* generated sequence, digitized */
  int   L;			/* length of sequence */
  int   alloc_tlen;		/* allocated space for traceback */
  int   alloc_L;		/* allocated space for sequence  */
  int   tpos;			/* position in traceback */
  int   sym;			/* a generated symbol index */
  float t[4];			/* little array for choosing M transition from */
  
  /* Initialize; allocations
   */
  P7AllocTrace(64, &tr);
  alloc_tlen = 64;
  dsq = MallocOrDie(sizeof(char) * 64);
  alloc_L = 64;

  TraceSet(tr, 0, STS, 0, 0);
  TraceSet(tr, 1, STN, 0, 0);
  dsq[0] = (char) Alphabet_iupac;
  L      = 1;
  k      = 0;
  type   = STN;
  tpos   = 2;

  while (type != STT) 
    {
      /* Deal with state transition
       */
      switch (type) {
      case STB:
	hmm->begin[0] = hmm->tbd1; /* begin[0] hack (documented in structs.h) */
	k = FChoose(hmm->begin, hmm->M+1);
	if (k == 0) { type = STD; k = 1; } else {type = STM; }
	break;

      case STI:	type = (FChoose(hmm->t[k]+TIM, 2) == 0)    ? STM : STI; if (type == STM) k++; break;
      case STN: type = (FChoose(hmm->xt[XTN], 2)  == LOOP) ? STN : STB; k = 0; break;
      case STE:	type = (FChoose(hmm->xt[XTE], 2)  == LOOP) ? STJ : STC; k = 0; break;
      case STC:	type = (FChoose(hmm->xt[XTC], 2)  == LOOP) ? STC : STT; k = 0; break;
      case STJ:	type = (FChoose(hmm->xt[XTJ], 2)  == LOOP) ? STJ : STB; k = 0; break;

      case STD:	
	if (k < hmm->M) {
	  type = (FChoose(hmm->t[k]+TDM, 2) == 0) ? STM : STD; 
	  k++;   
	} else {
	  type = STE;
	  k = 0;
	}
	break;

      case STM:
	if (k < hmm->M) {
	  FCopy(t, hmm->t[k], 3);
	  t[3] = hmm->end[k];
	  switch (FChoose(t,4)) {
	  case 0: k++;  type = STM; break;
	  case 1:       type = STI; break;
	  case 2: k++;  type = STD; break;
	  case 3: k=0;  type = STE; break;
	  default: Die("never happens");
	  }
	} else {
	  k    = 0;
	  type = STE;
	}
	break;

      case STT:
      case STBOGUS:
      default:
	Die("can't happen.");
      }
  
      /* Choose a symbol emission, if necessary
       */
      sym = -1;
      if      (type == STM) sym = FChoose(hmm->mat[k], Alphabet_size);
      else if (type == STI) sym = FChoose(hmm->ins[k], Alphabet_size); 
      else if ((type == STN && tr->statetype[tpos-1] == STN) ||
	       (type == STC && tr->statetype[tpos-1] == STC) ||
	       (type == STJ && tr->statetype[tpos-1] == STJ))
	sym = FChoose(hmm->null, Alphabet_size);
	
      /* Add to the traceback; deal with realloc if necessary
       */
      TraceSet(tr, tpos, type, k, (sym != -1) ? L : 0);
      tpos++;
      if (tpos == alloc_tlen) {
	alloc_tlen += 64; 
	P7ReallocTrace(tr, alloc_tlen);
      }

      /* Add to the digitized seq; deal with realloc, if necessary
       */
      if (sym != -1) {
	dsq[L] = (char) sym;
	L++;
	if (L+1 == alloc_L) {	/* L+1 leaves room for sentinel byte + \0 */
	  alloc_L += 64;
	  dsq = ReallocOrDie(dsq, sizeof(char) * alloc_L);
	}
      }
    }
  
  /* Finish off the trace
   */ 
  tr->tlen = tpos;

  /* Finish off the dsq with sentinel byte and null terminator.
   * Emitted Sequence length is L-1.
   */
  dsq[L]   = (char) Alphabet_iupac;
  dsq[L+1] = '\0';
  L--;

  /* Return
   */
  if (ret_dsq != NULL) *ret_dsq = dsq; else free(dsq);
  if (ret_L   != NULL) *ret_L   = L;
  if (ret_tr  != NULL) *ret_tr  = tr;  else P7FreeTrace(tr);
  return;
}

#ifdef SRE_REMOVED
/* Function: EmitBestSequence()
 * Date:     SRE, Tue Nov 10 16:21:59 1998 [St. Louis]
 *
 * Purpose:  Given a model, emit the maximum probability sequence
 *           from it: argmax_{seq} P(seq | model).
 *           This is a sensible HMM equivalent to a "consensus"
 *           sequence.
 *           The model should be Plan7NakedConfig()'ed; 
 *           in particular, if we allowed B->M and M->E,
 *           the highest probability sequence would be
 *           artifactually short. (We could do the highest
 *           scoring sequence instead, to get around this problem,
 *           but the highest scoring sequence is prone to
 *           other artifacts -- any looping state N,C,J, or I
 *           with a positively scoring residue leads to
 *           an infinitely long "best scoring" sequence.)
 *
 * Args:     hmm     - the model
 *           ret_seq - RETURN: best sequence
 *           ret_L   - RETURN: length of sequence
 *           ret_tr  - RETURN: traceback of the model/seq alignment; or NULL.
 *
 * Returns:  void
 */
void
EmitBestSequence(struct plan7_s *hmm, char **ret_dsq, int *ret_L, struct p7trace_s **ret_tr)
{
  char              *seq;                  /* RETURN: best seq */
  struct p7trace_s  *tr;                   /* RETURN: traceback */
  float             *mmx, *imx, *dmx;      /* log P forward scores for M,D,I */
  char              *mtb, *itb, *dtb;      /* traceback ptrs for M,D,I */
  int  x;		        /* counter for symbols */
  int  k;			/* counter for nodes   */
  float sc;			/* tmp var for a log P */
  int  bestsym;
  int  rpos;			/* position in a sequence */
  int  tpos;			/* position in a trace */
  int  tlen;			/* length of the traceback */

  /* Initial allocations. We only need a 1D matrix and its shadow;
   * it's overkill to use the Plan7Matrix structures, so don't.
   */
  mmx = MallocOrDie(sizeof(float) * (hmm->M+1));
  imx = MallocOrDie(sizeof(float) * (hmm->M));
  dmx = MallocOrDie(sizeof(float) * (hmm->M));
  mtb = MallocOrDie(sizeof(char)  * (hmm->M+1));
  itb = MallocOrDie(sizeof(char)  * (hmm->M));
  dtb = MallocOrDie(sizeof(char)  * (hmm->M));

  /* Initialization. 
   * We can safely assume a max probability path of S->N->B->(M1 or D1),
   * so just init M1 and D1.
   */
  mmx[1] = log(hmm->xt[XTN][MOVE]) + log(1.F - hmm->tbd1);
  dmx[1] = 


  /* Main recursion, done as a push.
   * The model is used in probability form; no wing folding needed.
   */
  for (k = 1; k < hmm->M; k++)
    {
      /* Transits out of match state (init with these)
       */
      mmx[k+1] = mmx[k] + log(hmm->t[k][TMM]); mtb[k+1] = STM;
      dmx[k+1] = mmx[k] + log(hmm->t[k][TMD]); dtb[k+1] = STM;
      if (k < hmm->M-1) 
	imx[k]   = mmx[k] + log(hmm->t[k][TMI]); itb[k]   = STM;
      
      /* Transits out of delete state
       */
      if ((sc = dmx[k] + log(hmm->t[k][TDM])) > mmx[k+1]) 
	{ mmx[k+1] = sc; mtb[k+1] = STD; }
      if ((sc = dmx[k] + log(hmm->t[k][TDD])) > dmx[k+1])
	{ dmx[k+1] = sc; dtb[k+1] = STD; }

      /* Transits out of insert state (self-loops are never good)
       */
      if ((sc = imx[k] + log(hmm->t[k][TIM])) > mmx[k+1])
	{ mmx[k+1] = sc; mtb[k+1] = STI; }
      
      /* Best emissions
       */
      x = FArgMax(hmm->mat[k+1], Alphabet_size);
      mmx[k+1] += log(hmm->mat[k+1][x]);

      if (k < hmm->M-1) {
	x = FArgMax(hmm->ins[k+1], Alphabet_size);
	imx[k+1] += log(hmm->ins[k+1][x]);
      }
    }
}
#endif /* SRE_REMOVED */


/* Function: EmitConsensusSequence()
 * Date:     SRE, Wed Nov 11 11:08:59 1998 [St. Louis]
 *
 * Purpose:  Generate a "consensus sequence". For the purposes
 *           of a profile HMM, this is defined as:
 *              - for each node:
 *                 - if StateOccupancy() says that M is used 
 *                     with probability >= 0.5, this M is "consensus".
 *                     Then, choose maximally likely residue.
 *                     if P>0.5 (protein) or P>0.9 (DNA), make
 *                     it upper case; else make it lower case. 
 *                 - if StateOccupancy() says that I
 *                     is used with P >= 0.5, this I is "consensus";
 *                     use it 1/(1-TII) times (its expectation value).
 *                     Generate an "x" from each I.
 *                     
 *           The function expects that the model is config'ed
 *           by Plan7NakedConfig(): that is, for a single global pass
 *           with no N,C,J involvement.
 *                     
 *
 * Args:     hmm     - the model
 *           ret_seq - RETURN: consensus sequence (pass NULL if unwanted)
 *           ret_dsq - RETURN: digitized consensus sequence (pass NULL if unwanted)
 *           ret_L   - RETURN: length of generated sequence 
 *           ret_tr  - RETURN: generated trace (pass NULL if unwanted)
 *
 * Returns:  void        
 */
void
EmitConsensusSequence(struct plan7_s *hmm, char **ret_seq, char **ret_dsq, int *ret_L, struct p7trace_s **ret_tr)
{
  struct p7trace_s *tr;         /* RETURN: traceback */
  char *dsq, *seq;              /* sequence in digitized and undigitized form */
  float *mp, *ip, *dp;          /* state occupancies from StateOccupancy() */
  int    nmat, ndel, nins;	/* number of matches, deletes, inserts used */
  int    k;			/* counter for nodes */
  int    tpos;			/* position in trace */
  int    i;                     /* position in seq (equiv pos in dsq is i+1 */
  int    x;			/* symbol choice (M) or # symbols (I) */
  float  mthresh;		/* >= this, show symbol as upper case */

  if (Alphabet_type == hmmAMINO) mthresh = 0.5;
  else                           mthresh = 0.9;

  StateOccupancy(hmm, &mp, &ip, &dp);

  /* First pass: how many states do we need in the trace?
   *             how long will the sequence be?
   */
  nmat = ndel = nins = 0;
  for (k = 1; k <= hmm->M; k++)
    {
      if (mp[k] >= 0.5) nmat++; else ndel++;
      if (k < hmm->M && ip[k] >= 0.5) 
	nins += (int) (1.f / (1.f - hmm->t[k][TII]));
    }
  
  /* Allocations
   */
  P7AllocTrace(6 + nmat + ndel + nins, &tr);
  dsq = MallocOrDie(sizeof(char) * (nmat+nins+3));
  seq = MallocOrDie(sizeof(char) * (nmat+nins+1));

  /* Main pass.
   * Construct consensus trace, seq, and dsq.
   */
  TraceSet(tr, 0, STS, 0, 0);
  TraceSet(tr, 1, STN, 0, 0);
  TraceSet(tr, 2, STB, 0, 0);
  dsq[0] = Alphabet_iupac;	/* guard byte */
  tpos = 3;
  i    = 0;
  for (k = 1; k <= hmm->M; k++)
    {
      if (mp[k] >= 0.5)
	{
	  x = FArgMax(hmm->mat[k], Alphabet_size);
	  TraceSet(tr, tpos, STM, k, i+1);
	  seq[i]   = Alphabet[x];
	  dsq[i+1] = x;
	  if (hmm->mat[k][x] < mthresh)
	    seq[i] = tolower((int) seq[i]);
	  i++;
	  tpos++;
	}
      else
	{
	  TraceSet(tr, tpos, STD, k, 0);
	  tpos++;
	}

      if (k < hmm->M && ip[k] >= 0.5)
	{
	  x = (int) (1.f / (1.f - hmm->t[k][TII]));
	  while (x--) 
	    {
	      TraceSet(tr, tpos, STI, k, i+1);
	      seq[i]   = 'x';
	      dsq[i+1] = Alphabet_iupac - 1;
	      i++; 
	      tpos++;
	    }
	}
    }
  TraceSet(tr, tpos, STE, 0, 0); tpos++;
  TraceSet(tr, tpos, STC, 0, 0); tpos++;
  TraceSet(tr, tpos, STT, 0, 0); tpos++;
  dsq[i+1] = Alphabet_iupac;
    
  free(mp);
  free(ip);
  free(dp);
  if (ret_seq != NULL) *ret_seq = seq; else free(seq);
  if (ret_dsq != NULL) *ret_dsq = dsq; else free(dsq);
  if (ret_L   != NULL) *ret_L   = i;   
  if (ret_tr  != NULL) *ret_tr  = tr;  else P7FreeTrace(tr);
}



/* Function: StateOccupancy()
 * Date:     SRE, Wed Nov 11 09:46:15 1998 [St. Louis]
 *
 * Purpose:  Calculate the expected state occupancy for
 *           a given HMM in generated traces.
 *           
 *           Note that expected prob of getting into
 *           any special state in a trace is trivial:
 *              S,N,B,E,C,T = 1.0
 *              J = E->J transition prob 
 *
 * Args:     hmm    - the model
 *           ret_mp - RETURN: [1..M] prob's of occupying M
 *           ret_ip - RETURN: [1..M-1] prob's of occupying I
 *           ret_dp - RETURN: [1..M] prob's of occupying D
 *
 * Returns:  void
 *           mp, ip, dp are malloc'ed here. Caller must free().
 */
void
StateOccupancy(struct plan7_s *hmm, float **ret_mp, float **ret_ip, float **ret_dp)
{
  float *fmp, *fip, *fdp;       /* forward probabilities  */
  int k;			/* counter for nodes      */

  /* Initial allocations
   */
  fmp = MallocOrDie (sizeof(float) * (hmm->M+1));
  fip = MallocOrDie (sizeof(float) * (hmm->M));
  fdp = MallocOrDie (sizeof(float) * (hmm->M+1));
  
  /* Forward pass. 
   */
  fdp[1] = hmm->tbd1;
  fmp[1] = hmm->begin[1];
  fip[1] = fmp[1] * hmm->t[1][TMI];
  for (k = 2; k <= hmm->M; k++)
    {
			/* M: from M,D,I at k-1, or B; count t_II as 1.0 */
      fmp[k] = fmp[k-1] * hmm->t[k-1][TMM] +
	       fip[k-1] +
               fdp[k-1] * hmm->t[k-1][TDM] +
	       hmm->begin[k];
			/* D: from M,D at k-1 */
      fdp[k] = fmp[k-1] * hmm->t[k-1][TMD] +
	       fdp[k-1] * hmm->t[k-1][TDD];
			/* I: from M at k; don't count II */
      if (k < hmm->M) {
	fip[k] = fmp[k] * hmm->t[k][TMI];
      }

      SQD_DASSERT2((fabs(1.0f - fmp[k] - fdp[k]) < 1e-6f));
      fmp[k] /= fmp[k]+fdp[k];	/* prevent propagating fp errors */
      fdp[k] /= fmp[k]+fdp[k];
    }
  /* We don't need a backward pass; all backwards P's are 1.0
   * by definition (you can always get out of a state with P=1).
   * The only situation where this might not be true is for
   * a TII of 1.0, when TIM = 0 -- but in that case, if there's
   * a finite chance of getting into that insert state, the model
   * generates infinitely long sequences, so we can consider this
   * situation "perverse" and disallow it elsewhere in building
   * profile HMMs.
   */

  /* Return.
   */
  *ret_mp = fmp;
  *ret_dp = fdp;
  *ret_ip = fip;
}
