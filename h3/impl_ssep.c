/* SSE implementation of Forward filter
 * 
 * SRE, Thu Dec 13 08:42:56 2007 [Janelia]
 * SVN $Id$
 */

#include "p7_config.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_alphabet.h"

#include "hmmer.h"
#include "impl_ssep.h"


/*****************************************************************
 * 1. The P7_OPROFILE structure: a score profile.
 *****************************************************************/

/* Function:  p7_oprofile_Create()
 * Synopsis:  Allocate an optimized profile structure.
 * Incept:    SRE, Thu Dec 13 08:49:04 2007 [Janelia]
 *
 * Purpose:   Create a profile of <M> nodes for digital alphabet <abc>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_OPROFILE *
p7_oprofile_Create(int M, const ESL_ALPHABET *abc)
{
  int          status;
  P7_OPROFILE *om = NULL;
  int          x;
  int          nq = p7O_NQ(M);	/* # of "blocks" (aligned, striped 128-bit units) needed for query */

  /* level 0 */
  ESL_ALLOC(om, sizeof(P7_OPROFILE));
  om->tp     = NULL;
  om->rp     = NULL;

  /* level 1 */
  ESL_ALLOC(om->tp, sizeof(__m128) * nq  * p7O_NTRANS);    
  ESL_ALLOC(om->rp, sizeof(__m128 *) * abc->Kp); 
  om->rp[0] = NULL;

  /* level 2 */
  ESL_ALLOC(om->rp[0], sizeof(__m128) * nq  * p7O_NR * abc->Kp);                     
  for (x = 1; x < abc->Kp; x++)
    om->rp[x] = om->rp[0] + (x * nq * p7O_NR);

  /* remaining initializations */
  om->mode   = p7_NO_MODE;
  om->M      = M;
  om->abc    = abc;
  om->allocQ = nq;		/* remember how big our alloc is: allows future reuse */
  return om;

 ERROR:
  p7_oprofile_Destroy(om);
  return NULL;
}


/* Function:  p7_oprofile_Destroy()
 * Synopsis:  Frees an optimized profile structure.
 * Incept:    SRE, Thu Dec 13 08:49:13 2007 [Janelia]
 */
void
p7_oprofile_Destroy(P7_OPROFILE *om)
{
  if (om == NULL) return;

  if (om->tp != NULL) free(om->tp);
  if (om->rp != NULL) 
    {
      if (om->rp[0] != NULL) free(om->rp[0]);
      free(om->rp);
    }
  free(om);
}


/* Function:  p7_oprofile_Convert()
 * Synopsis:  Converts standard profile to an optimized one.
 * Incept:    SRE, Thu Dec 13 08:49:21 2007 [Janelia]
 */
int
p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om)
{
  int x;			/* counter over residues */
  int q;			/* q counts over total # of striped quads, always 1.. p7O_NQ */
  int j;			/* j counts over position of quads in some particular memory arrangement */
  int k;			/* the usual counter over model nodes 1..M */
  int kb;			/* base k for loading om's TSC quads */
  int t;			/* counter over transitions 0..7 = p7O_{BM, MM, IM, DM, MD, MI, II, DD} */
  int tg;			/* transition index in gm */
  int M  = gm->M;		/* length of the query */
  int nq = p7O_NQ(M);	        /* segment length; total # of striped quads */


  if (gm->abc_r->type != om->abc->type)  ESL_EXCEPTION(eslEINVAL, "alphabets of the two profiles don't match");
  if (nq > om->allocQ)                   ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  /* match, insert emission scores: start at k=1 */
  for (x = 0; x < gm->abc_r->Kp; x++)
    for (j = 0, k = 1, q = 0; q < nq; q++, k++)
      {
	om->rp[x][j++] = _mm_setr_ps((k+0*nq <= M) ? exp(p7P_MSC(gm, k,      x)) : 0.,
				     (k+1*nq <= M) ? exp(p7P_MSC(gm, k+nq,   x)) : 0.,
				     (k+2*nq <= M) ? exp(p7P_MSC(gm, k+2*nq, x)) : 0.,
				     (k+3*nq <= M) ? exp(p7P_MSC(gm, k+3*nq, x)) : 0.);

	om->rp[x][j++] = _mm_setr_ps((k+0*nq <= M) ? exp(p7P_ISC(gm, k,      x)) : 0.,
				     (k+1*nq <= M) ? exp(p7P_ISC(gm, k+nq,   x)) : 0.,
				     (k+2*nq <= M) ? exp(p7P_ISC(gm, k+2*nq, x)) : 0.,
				     (k+3*nq <= M) ? exp(p7P_ISC(gm, k+3*nq, x)) : 0.);
      }

  /* Transition scores, all but the DD's. */
  for (j = 0, k = 1, q = 0; q < nq; q++, k++)
    {
      for (t = p7O_BM; t <= p7O_II; t++) /* this loop of 7 transitions depends on the order in the definition of p7o_tsc_e */
	{
	  switch (t) {
	  case p7O_BM: tg = p7P_BM;  kb = k-1; break; /* gm has tBMk stored off by one! start from k=0 not 1 */
	  case p7O_MM: tg = p7P_MM;  kb = k-1; break; /* MM, DM, IM quads are rotated by -1, start from k=0  */
	  case p7O_IM: tg = p7P_IM;  kb = k-1; break;
	  case p7O_DM: tg = p7P_DM;  kb = k-1; break;
	  case p7O_MD: tg = p7P_MD;  kb = k;   break; /* the remaining ones are straight up  */
	  case p7O_MI: tg = p7P_MI;  kb = k;   break; 
	  case p7O_II: tg = p7P_II;  kb = k;   break; 
	  }

	  om->tp[j++] = _mm_setr_ps((kb      < M) ? exp(p7P_TSC(gm, kb,      tg)) : 0.,
				    (kb+1*nq < M) ? exp(p7P_TSC(gm, kb+nq,   tg)) : 0.,
				    (kb+2*nq < M) ? exp(p7P_TSC(gm, kb+2*nq, tg)) : 0.,
				    (kb+3*nq < M) ? exp(p7P_TSC(gm, kb+3*nq, tg)) : 0.);
	}
    }

  /* And finally the DD's, which are at the end of the optimized tsc vector; (j is already there) */
  for (k = 1, q = 0; q < nq; q++, k++)
    om->tp[j++] = _mm_setr_ps((k+0*nq < M) ? exp(p7P_TSC(gm, k,      p7P_DD)) : 0.,
			      (k+1*nq < M) ? exp(p7P_TSC(gm, k+nq,   p7P_DD)) : 0.,
			      (k+2*nq < M) ? exp(p7P_TSC(gm, k+2*nq, p7P_DD)) : 0.,
			      (k+3*nq < M) ? exp(p7P_TSC(gm, k+3*nq, p7P_DD)) : 0.);

  /* Specials. (These are actually in exactly the same order in om and
   *  gm, but we copy in general form anyway.)
   */
  om->xp[p7O_E][p7O_LOOP] = exp(gm->xsc[p7P_E][p7P_LOOP]);  
  om->xp[p7O_E][p7O_MOVE] = exp(gm->xsc[p7P_E][p7P_MOVE]);
  om->xp[p7O_N][p7O_LOOP] = exp(gm->xsc[p7P_N][p7P_LOOP]);
  om->xp[p7O_N][p7O_MOVE] = exp(gm->xsc[p7P_N][p7P_MOVE]);
  om->xp[p7O_C][p7O_LOOP] = exp(gm->xsc[p7P_C][p7P_LOOP]);
  om->xp[p7O_C][p7O_MOVE] = exp(gm->xsc[p7P_C][p7P_MOVE]);
  om->xp[p7O_J][p7O_LOOP] = exp(gm->xsc[p7P_J][p7P_LOOP]);
  om->xp[p7O_J][p7O_MOVE] = exp(gm->xsc[p7P_J][p7P_MOVE]);

  om->mode = gm->mode;
  om->M    = M;
  return eslOK;
}

/* Function:  p7_oprofile_Dump()
 * Synopsis:  Dump internals of a <P7_OPROFILE>
 * Incept:    SRE, Thu Dec 13 08:49:30 2007 [Janelia]
 *
 * Purpose:   Dump the internals of <P7_OPROFILE> structure <om>
 *            to stream <fp>; generally for testing or debugging
 *            purposes.
 *
 * Args:      fp   - output stream (often stdout)
 *            om   - optimized profile to dump
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_oprofile_Dump(FILE *fp, P7_OPROFILE *om)
{
  int x;			/* counter over residues */
  int q;			/* q counts over total # of striped quads, always 1.. p7O_NQ */
  int j;			/* j counts over position of quads in some particular memory arrangement */
  int z;
  int k;			/* the usual counter over model nodes 1..M */
  int kb;
  int t;
  int M  = om->M;		/* length of the profile */
  int nq = p7O_NQ(M);	        /* segment length; total # of striped quads */
  union { __m128 v; float x[p7O_QWIDTH]; } tmp;

  /* Residue emissions
   */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 
      for (k =1, q = 0; q < nq; q++, k++)
	{
	  fprintf(fp, "[ ");
	  for (z = 0; z < p7O_QWIDTH; z++) 
	    if (k+z*nq <= M) fprintf(fp, "%8d ", k+z*nq);
	    else             fprintf(fp, "%8s ", "xx");
	  fprintf(fp, "]");
	}

      fprintf(fp, "\nmat: ");
      for (j = 0, q = 0; q < nq; q++, j+=2)
	{
	  fprintf(fp, "[ ");
	  tmp.v = om->rp[x][j];
	  for (z = 0; z < p7O_QWIDTH; z++) fprintf(fp, "%8.5f ", tmp.x[z]);
	  fprintf(fp, "]");
	}

      fprintf(fp, "\nins: ");
      for (j = 1, q = 0; q < nq; q++, j+=2)
	{
	  fprintf(fp, "[ ");
	  tmp.v = om->rp[x][j];
	  for (z = 0; z < p7O_QWIDTH; z++) fprintf(fp, "%8.5f ", tmp.x[z]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n\n");
    }

  /* Transitions */
  for (t = p7O_BM; t <= p7O_II; t++)
    {
      switch (t) {
      case p7O_BM: fprintf(fp, "\ntBM: "); break;
      case p7O_MM: fprintf(fp, "\ntMM: "); break;
      case p7O_IM: fprintf(fp, "\ntIM: "); break;
      case p7O_DM: fprintf(fp, "\ntDM: "); break;
      case p7O_MD: fprintf(fp, "\ntMD: "); break;
      case p7O_MI: fprintf(fp, "\ntMI: "); break;
      case p7O_II: fprintf(fp, "\ntII: "); break;
      }

      for (k = 1, q = 0; q < nq; q++, k++)
	{
	  switch (t) {
	  case p7O_BM: kb = k;                 break; 
	  case p7O_MM: kb = 1 + (nq+k-2) % nq; break; /* MM, DM, IM quads rotated by +1  */
	  case p7O_IM: kb = 1 + (nq+k-2) % nq; break;  
	  case p7O_DM: kb = 1 + (nq+k-2) % nq; break;  
	  case p7O_MD: kb = k;                 break; /* the remaining ones are straight up  */
	  case p7O_MI: kb = k;                 break; 
	  case p7O_II: kb = k;                 break; 
	  }
	  fprintf(fp, "[ ");
	  for (z = 0; z < p7O_QWIDTH; z++) 
	    if (kb+z*nq <= M) fprintf(fp, "%8d ", kb+z*nq);
	    else              fprintf(fp, "%8s ", "xx");
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n     ");	  
      for (q = 0; q < nq; q++)
	{
	  fprintf(fp, "[ ");
	  tmp.v = om->tp[q*7 + t];
	  for (z = 0; z < p7O_QWIDTH; z++) fprintf(fp, "%8.5f ", tmp.x[z]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");	  
    }

  /* DD transitions */
  fprintf(fp, "\ntDD: ");
  for (k =1, q = 0; q < nq; q++, k++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < p7O_QWIDTH; z++) 
	if (k+z*nq <= M) fprintf(fp, "%8d ", k+z*nq);
	else             fprintf(fp, "%8s ", "xx");
      fprintf(fp, "]");
    }
  fprintf(fp, "\n     ");	  
  for (j = nq*7, q = 0; q < nq; q++, j++)
    {
      fprintf(fp, "[ ");
      tmp.v = om->tp[j];
      for (z = 0; z < p7O_QWIDTH; z++) fprintf(fp, "%8.5f ", tmp.x[z]);
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");	  
  
  /* Specials */
  fprintf(fp, "E->C: %8.5f    E->J: %8.5f\n", om->xp[p7O_E][p7O_MOVE], om->xp[p7O_E][p7O_LOOP]);
  fprintf(fp, "N->B: %8.5f    N->N: %8.5f\n", om->xp[p7O_N][p7O_MOVE], om->xp[p7O_N][p7O_LOOP]);
  fprintf(fp, "J->B: %8.5f    J->J: %8.5f\n", om->xp[p7O_J][p7O_MOVE], om->xp[p7O_J][p7O_LOOP]);
  fprintf(fp, "C->T: %8.5f    C->C: %8.5f\n", om->xp[p7O_C][p7O_MOVE], om->xp[p7O_C][p7O_LOOP]);

  fprintf(fp, "Q:     %d\n",   nq);  
  fprintf(fp, "M:     %d\n",   M);  
  return eslOK;
}




/*****************************************************************
 * 2. The P7_OMX structure: a dynamic programming matrix
 *****************************************************************/

/* Function:  p7_omx_Create()
 * Synopsis:  Create an optimized dynamic programming matrix.
 * Incept:    SRE, Thu Dec 13 08:52:46 2007 [Janelia]
 *
 * Purpose:   Allocates a reusable, resizeable <P7_OMX> for models up to
 *            size <allocM>.
 *
 * Returns:   a pointer to the new <P7_OMX>.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_OMX *
p7_omx_Create(int allocM)
{
  int      status;
  P7_OMX  *ox;
  int      nq = p7O_NQ(allocM);	       /* segment length; total # of striped quads */

  ESL_ALLOC(ox, sizeof(P7_OMX));
  ox->dp = NULL;

  ESL_ALLOC(ox->dp,  sizeof(__m128) * p7X_NSCELLS * nq);
  ox->M      = 0;
  ox->Q      = 0;
  ox->allocQ = nq;
#ifdef p7_DEBUGGING
  ox->debugging = FALSE;
  ox->dfp       = NULL;
#endif
  return ox;

 ERROR:
  p7_omx_Destroy(ox);
  return NULL;
}


/* Function:  p7_omx_Destroy()
 * Synopsis:  Frees an optimized DP matrix.
 * Incept:    SRE, Thu Dec 13 08:52:53 2007 [Janelia]
 *
 * Purpose:   Frees optimized DP matrix <ox>.
 *
 * Returns:   (void)
 */
void
p7_omx_Destroy(P7_OMX *ox)
{
  if (ox == NULL) return;
  if (ox->dp != NULL) free(ox->dp);
  free(ox);
  return;
}

/* Function:  p7_omx_SetDumpMode()
 * Synopsis:  Set an optimized DP matrix to be dumped for debugging.
 * Incept:    SRE, Thu Dec 13 10:24:38 2007 [Janelia]
 *
 * Purpose:   Whenever a dynamic programming calculation is run, 
 *            dump DP matrix <ox> to stream <fp> for diagnostics. 
 *            
 *            When the dump mode is on, the DP routine itself actually
 *            does the dumping, because it has to dump after every row
 *            is calculated. (We're doing an optimized one-row
 *            calculation.)
 *            
 *            If the code has not been compiled with the
 *            <p7_DEBUGGING> flag up, this function is a no-op.
 *
 * Args:      fp   - output stream for diagnostics (stdout, perhaps)
 *            ox   - DP matrix to set to debugging mode
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      J2/62.
 */
int
p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox)
{
#if p7_DEBUGGING
  ox->debugging = TRUE;
  ox->dfp       = fp;
#endif
  return eslOK;
}

#ifdef p7_DEBUGGING
/* omx_dump_row()
 * Dump current row of optimized DP matrix for debugging.
 * SRE, Thu Dec 13 08:53:01 2007 [Janelia]
 *
 * Dump DP matrix <ox> current row for diagnostics The matrix contains
 * values for one current row <rowi>, and this is used as a row label.
 * 
 * If <rowi> is 0, it prints a header first too.
 * 
 * The output format is coordinated with <p7_gmx_Dump()> to
 * facilitate comparison to a known answer.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
static int
omx_dump_row(P7_OMX *ox, int rowi, float xE, float xN, float xJ, float xB, float xC)
{
  int     q,z,k;
  float  *v;
  float   tmp[4];
  __m128 *dp = ox->dp;
  int     Q  = ox->Q;
  int     M  = ox->M;
  int     status;

  ESL_ALLOC(v, sizeof(float) * ((Q*4)+1));
  v[0] = 0.;

  if (rowi == 0)
    {
      fprintf(ox->dfp, "      ");
      for (k = 0; k <= M;  k++) fprintf(ox->dfp, "%8d ", k);
      fprintf(ox->dfp, "%8s %8s %8s %8s %8s\n", "E", "N", "J", "B", "C");
      fprintf(ox->dfp, "      ");
      for (k = 0; k <= M+5;  k++) fprintf(ox->dfp, "%8s ", "--------");
      fprintf(ox->dfp, "\n");
    }

  /* Unpack, unstripe, then print M's. */
  for (q = 0; q < Q; q++) {
    _mm_store_ps(tmp, MMX(q));
    for (z = 0; z < 4; z++) v[q+Q*z+1] = tmp[z];
  }
  fprintf(ox->dfp, "%3d M ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%8.5f ", v[k] == 0. ? -eslINFINITY : log(v[k]));

 /* The specials */
  fprintf(ox->dfp, "%8.5f %8.5f %8.5f %8.5f %8.5f\n",
	  xE == 0. ? -eslINFINITY : log(xE),
	  xN == 0. ? -eslINFINITY : log(xN),
	  xJ == 0. ? -eslINFINITY : log(xJ),
	  xB == 0. ? -eslINFINITY : log(xB), 
	  xC == 0. ? -eslINFINITY : log(xC));

  /* Unpack, unstripe, then print I's. */
  for (q = 0; q < Q; q++) {
    _mm_store_ps(tmp, IMX(q));
    for (z = 0; z < 4; z++) v[q+Q*z+1] = tmp[z];
  }
  fprintf(ox->dfp, "%3d I ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%8.5f ", v[k] == 0. ? -eslINFINITY : log(v[k]));
  fprintf(ox->dfp, "\n");

  /* Unpack, unstripe, then print D's. */
  for (q = 0; q < Q; q++) {
    _mm_store_ps(tmp, DMX(q));
    for (z = 0; z < 4; z++) v[q+Q*z+1] = tmp[z];
  }
  fprintf(ox->dfp, "%3d D ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%8.5f ", v[k] == 0. ? -eslINFINITY : log(v[k]));
  fprintf(ox->dfp, "\n\n");

  free(v);
  return eslOK;

ERROR:
  free(v);
  return status;
}
#endif /*p7_DEBUGGING*/


/*****************************************************************
 * 3. The p7_ForwardFilter() DP implementation.
 *****************************************************************/



/* Function:  p7_ForwardFilter()
 * Synopsis:  Calculates Forward score, insanely fast.
 * Incept:    SRE, Thu Dec 13 08:54:07 2007 [Janelia]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      [Farrar07] for ideas behind striped SIMD DP.
 *            J2/46-47 for layout of HMMER's striped SIMD DP.
 *            J2/50 for single row DP.
 *            
 * Notes:     We use 7 xmm registers (I think). There are 16 xmm registers
 *            on 64-bit AMD and Intel (and 8 on IA-32). So on IA-64, we might
 *            try unrolling this and doing two calculations (double the chunk
 *            width) instead of one, as Lindahl did in the H2 Altivec code.  
 */
int
p7_ForwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register __m128 mpv, dpv, ipv;   /* previous row values */
  register __m128 sv;		   /* temp storage of 1 curr row value in progress */
  register __m128 dcv;		   /* delayed storage of D(i,q+1) */
  register __m128 xEv;		   /* E state: keeps max for Mk->E as we go */
  register __m128 xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  __m128   zerov;
  float    xN, xE, xB, xC, xJ;	   /* special states */
  int i;			   /* counter over sequence positions 1..L */
  int q;			   /* counter over quads 1..nq */
  int z;
  int Q       = p7O_NQ(om->M);	   /* segment length: # of vectors */
  __m128 *dp  = ox->dp;
  __m128 *rp;			   /* will point at om->rp[x] for residue x[i] */
  __m128 *tp;			   /* will point into (and step thru) om->tp   */


  /* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ) ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M = om->M;
  ox->Q = Q;

  /* Initialization.
   */
  zerov = _mm_setzero_ps();
  for (q = 0; q < Q; q++)
    MMX(q) = IMX(q) = DMX(q) = zerov;
  xE    = 0.;
  xN    = 1.;
  xJ    = 0.;
  xB    = om->xp[p7O_N][p7O_MOVE];
  xC    = 0.;


#if p7_DEBUGGING
  if (ox->debugging) omx_dump_row(ox, 0, xE, xN, xJ, xB, xC);   
#endif

  for (i = 1; i <= L; i++)
    {
      rp    = om->rp[dsq[i]];
      tp    = om->tp;
      dcv   = _mm_setzero_ps();
      xEv   = _mm_setzero_ps();
      xBv   = _mm_set1_ps(xB);

      /* Right shifts by 4 bytes. 4,8,12,x becomes x,4,8,12.  Shift zeros on.
       */
      mpv = MMX(Q-1);  mpv = _mm_shuffle_ps(mpv, mpv, _MM_SHUFFLE(2, 1, 0, 0));   mpv = _mm_move_ss(mpv, zerov);
      dpv = DMX(Q-1);  dpv = _mm_shuffle_ps(dpv, dpv, _MM_SHUFFLE(2, 1, 0, 0));   dpv = _mm_move_ss(dpv, zerov);
      ipv = IMX(Q-1);  ipv = _mm_shuffle_ps(ipv, ipv, _MM_SHUFFLE(2, 1, 0, 0));   ipv = _mm_move_ss(ipv, zerov);
      
      for (q = 0; q < Q; q++)
	{
	  /* Calculate new MMX(i,q); don't store it yet, hold it in sv. */
	  sv   =                _mm_mul_ps(xBv, *tp);  tp++;
	  sv   = _mm_add_ps(sv, _mm_mul_ps(mpv, *tp)); tp++;
	  sv   = _mm_add_ps(sv, _mm_mul_ps(ipv, *tp)); tp++;
	  sv   = _mm_add_ps(sv, _mm_mul_ps(dpv, *tp)); tp++;
	  sv   = _mm_mul_ps(sv, *rp);                  rp++;
	  xEv  = _mm_add_ps(xEv, sv);
	  
	  /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
	   * {MDI}MX(q) is then the current, not the prev row
	   */
	  mpv = MMX(q);
	  dpv = DMX(q);
	  ipv = IMX(q);

	  /* Do the delayed stores of {MD}(i,q) now that memory is usable */
	  MMX(q) = sv;
	  DMX(q) = dcv;

	  /* Calculate the next D(i,q+1) partially: M->D only;
           * delay storage, holding it in dcv
	   */
	  dcv   = _mm_mul_ps(sv, *tp); tp++;

	  /* Calculate and store I(i,q) */
	  sv     =                _mm_mul_ps(mpv, *tp);  tp++;
	  sv     = _mm_add_ps(sv, _mm_mul_ps(ipv, *tp)); tp++;
	  IMX(q) = _mm_mul_ps(sv, *rp);                  rp++;
	}	  

      /* Now the DD paths. We would rather not serialize them but 
       * in an accurate Forward calculation, we have few options.
       */
      /* dcv has carried through from end of q loop above; store it 
       * in first pass, we add M->D and D->D path into DMX
       */
      /* We're almost certainly're obligated to do at least one complete 
       * DD path to be sure: 
       */
      dcv    = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
      dcv    = _mm_move_ss(dcv, zerov);
      DMX(0) = zerov;
      tp     = om->tp + 7*Q;	/* set tp to start of the DD's */
      for (q = 0; q < Q; q++) 
	{
	  DMX(q) = _mm_add_ps(dcv, DMX(q));	
	  dcv    = _mm_mul_ps(DMX(q), *tp); tp++; /* extend DMX(q), so we include M->D and D->D paths */
	}

      /* now. on small models, it seems best (empirically) to just go
       * ahead and serialize. on large models, we can do a bit better,
       * by testing for when dcv (DD path) accrued to DMX(q) is below
       * machine epsilon for all q, in which case we know DMX(q) are all
       * at their final values. The tradeoff point is (empirically) somewhere around M=100,
       * at least on my desktop. We don't worry about the conditional here;
       * it's outside any inner loops.
       */
      if (om->M < 100)
	{			/* Fully serialized version */
	  for (z = 1; z < 4; z++)
	    {
	      dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
	      dcv = _mm_move_ss(dcv, zerov);
	      tp = om->tp + 7*Q;	/* set tp to start of the DD's */
	      for (q = 0; q < Q; q++) 
		{
		  DMX(q) = _mm_add_ps(dcv, DMX(q));	
		  dcv    = _mm_mul_ps(dcv, *tp);   tp++; /* note, extend dcv, not DMX(q); only adding DD paths now */
		}	    
	    }
	} 
      else
	{			/* Slightly parallelized version, but which incurs some overhead */
	  for (z = 1; z < 4; z++)
	    {
	      register __m128 cv;	/* keeps track of whether any DD's change DMX(q) */

	      dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
	      dcv = _mm_move_ss(dcv, zerov);
	      tp  = om->tp + 7*Q;	/* set tp to start of the DD's */
	      cv  = zerov;
	      for (q = 0; q < Q; q++) 
		{
		  sv     = _mm_add_ps(dcv, DMX(q));	
		  cv     = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMX(q))); /* remember if DD paths changed any DMX(q): *without* conditional branch */
		  DMX(q) = sv;	                            /* store new DMX(q) */
		  dcv    = _mm_mul_ps(dcv, *tp);   tp++;        /* note, extend dcv, not DMX(q); only adding DD paths now */
		}	    
	      if (! _mm_movemask_ps(cv)) break; /* DD's didn't change any DMX(q)? Then we're done, break out. */
	    }
	}

      /* Add D's to xEv */
      for (q = 0; q < Q; q++) 
	xEv = _mm_add_ps(DMX(q), xEv);

      /* Finally the "special" states, which start from Mk->E (->C, ->J->B) */
      /* The following incantation is a horizontal sum of xEv's elements  */
      /* These must follow DD calculations, because D's contribute to E in Forward
       * (as opposed to Viterbi)
       */
      xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(0, 3, 2, 1)));
      xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(1, 0, 3, 2)));
      _mm_store_ss(&xE, xEv);

      xN = xN * om->xp[p7O_N][p7O_LOOP];
      xC = (xC * om->xp[p7O_C][p7O_LOOP]) +  (xE * om->xp[p7O_E][p7O_MOVE]);
      xJ = (xJ * om->xp[p7O_J][p7O_LOOP]) +  (xE * om->xp[p7O_E][p7O_LOOP]);
      xB = (xJ * om->xp[p7O_J][p7O_MOVE]) +  (xN * om->xp[p7O_N][p7O_MOVE]);
      /* and now xB will carry over into next i, and xC carries over after i=L */


#if p7_DEBUGGING
      if (ox->debugging) omx_dump_row(ox, i, xE, xN, xJ, xB, xC);   
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T, and flip back to log space */
  *ret_sc  = log(xC * om->xp[p7O_C][p7O_MOVE]);
  return eslOK;
}


/*****************************************************************
 * Exam drivers.
 *****************************************************************/

#ifdef p7IMPL_SSEP_EXAM
/* Examining internals of a converted optimized profile
   gcc -msse2 -g -Wall -I. -L. -I../easel -L../easel -o exam-impl-ssep -Dp7IMPL_SSEP_EXAM impl_ssep.c -lhmmer -leasel -lm
   ./exam <hmmfile>
   Pfam 22.0 Bombesin is a good small example (M=14)
 */ 
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "examine internals of a P7_OPROFILE after conversion";


int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  int             L       = 400;

  /* Read a test HMM in, make a local profile out of it */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM from %s", hmmfile);
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);

  /* Convert profile to an optimized profile, and dump it */
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  p7_oprofile_Dump(stdout, om);
  
  /* Clean up and return */
  p7_oprofile_Destroy(om);
  p7_hmm_Destroy(hmm);      
  p7_profile_Destroy(gm);
  p7_hmmfile_Close(hfp);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7IMPL_SSEP_EXAM*/


#ifdef p7IMPL_SSEP_EXAM2
/* Comparing internals of small DP matrices
   gcc -msse2 -g -Wall -I. -L. -I../easel -L../easel -o exam -Dp7IMPL_SSEP_EXAM2 impl_ssep.c -lhmmer -leasel -lm
   ./exam <hmmfile> <seqfile>

   Example alignment: 
# STOCKHOLM 1.0
foo1  GAATTC
foo2  GAATTC
//
   Example sequence:
>foo3
GAATTC
 */ 
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",          0 },
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "output internals of dp_generic vs. SSE DP matrices for debugging";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_OMX         *ox      = NULL;
  P7_GMX         *gx      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           sc1, sc2;
  int             status;

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
  if  (esl_sqio_Read(sqfp, sq) != eslOK) p7_Fail("Failed to read sequence");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  p7_oprofile_Dump(stdout, om);
 
  ox = p7_omx_Create(gm->M);
  gx = p7_gmx_Create(gm->M, sq->n);
  p7_omx_SetDumpMode(stdout, ox);

  p7_ForwardFilter(sq->dsq, sq->n, om, ox, &sc1); 
  p7_GForward     (sq->dsq, sq->n, gm, gx, &sc2);
  p7_gmx_Dump(stdout, gx);

  printf("raw score from SSE:        %f\n", sc1);
  printf("raw score from dp_generic: %f\n", sc2);

  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7IMPL_SSEP_EXAM2*/



/*****************************************************************
 * Benchmark driver.
 *****************************************************************/
#ifdef p7IMPL_SSEP_BENCHMARK
/* gcc -o benchmark-ssep -g -O3 -msse2 -I. -L. -I../easel -L../easel -Dp7IMPL_SSEP_BENCHMARK impl_ssep.c -lhmmer -leasel -lm
 * icc -o benchmark-ssep -O3 -static -I. -L. -I../easel -L../easel -Dp7IMPL_SSEP_BENCHMARK impl_ssep.c -lhmmer -leasel -lm 
 * icc -o benchmark-ssep -g -Wall -static -I. -L. -I../easel -L../easel -Dp7IMPL_SSEP_BENCHMARK impl_ssep.c -lhmmer -leasel -lm 
 * 
 *   ./benchmark-ssep <hmmfile>         runs benchmark
 *   ./benchmark-ssep -b <hmmfile>      gets baseline time to subtract: just random seq generation
 *   ./benchmark-ssep -c <hmmfile>      compare scores of SSE to generic impl
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                  0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose: show individual scores",               0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "compare scores of generic, SSE DP        (debug)", 0 }, 
  { "-b",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "baseline timing: don't run DP at all",             0 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "time the generic (serial) version",                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for the generic implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_OMX         *ox      = NULL;
  P7_GMX         *gx      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc1, sc2;

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  ox = p7_omx_Create(gm->M);
  gx = p7_gmx_Create(gm->M, L);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rnd_xfIID(r, bg->f, abc->K, L, dsq);

      if (! esl_opt_GetBoolean(go, "-b")) {
	if (esl_opt_GetBoolean(go, "-g")) p7_GForward     (dsq, L, gm, gx, &sc1);   
	else                              p7_ForwardFilter(dsq, L, om, ox, &sc1);   

	if (esl_opt_GetBoolean(go, "-c")) {
	  p7_GForward     (dsq, L, gm, gx, &sc2); 
	  printf("%.4f %.4f\n", sc1, sc2);  
	}
      }
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M = %d\n", gm->M);

  free(dsq);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
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
#endif 
