/* SSE implementation of Viterbi filter in compressed precision.
 * 
 *   1. The P7_OPROFILE structure: a score profile.
 *   2. The P7_OMX structure: a dynamic programming matrix.
 *   3. Viterbi filter implementation.
 * 
 * SRE, Sun Dec  9 12:04:30 2007 [Janelia]
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
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_sse8.h"


/*****************************************************************
 * 1. The P7_OPROFILE structure: a score profile.
 *****************************************************************/

/* Function:  p7_oprofile_Create()
 * Synopsis:  Allocate an optimized profile structure.
 * Incept:    SRE, Sun Dec  9 12:05:29 2007 [Janelia]
 *
 * Purpose:   Allocate a profile of <M> nodes for digital alphabet <abc>.
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
  om->tsc     = NULL;
  om->rsc     = NULL;

  /* level 1 */
  ESL_ALLOC(om->tsc, sizeof(__m128i) * nq  * p7O_NTRANS);    
  ESL_ALLOC(om->rsc, sizeof(__m128i *) * abc->Kp); 
  om->rsc[0] = NULL;

  /* level 2 */
  ESL_ALLOC(om->rsc[0], sizeof(__m128i) * nq  * p7O_NR * abc->Kp);                     
  for (x = 1; x < abc->Kp; x++)
    om->rsc[x] = om->rsc[0] + (x * nq * p7O_NR);

  /* remaining initializations */
  om->abc      = abc;
  om->mode     = p7_NO_MODE;
  om->M        = M;
  om->allocQ   = nq;		/* remember how big our alloc is: allows future reuse */
  om->dd_bound = 0;
  om->scale    = 0.0;
  om->base     = 0;
  om->bias     = 0;
  return om;

 ERROR:
  p7_oprofile_Destroy(om);
  return NULL;
}

/* Function:  p7_oprofile_Destroy()
 * Synopsis:  Frees an optimized profile structure.
 * Incept:    SRE, Sun Dec  9 12:05:37 2007 [Janelia]
 */
void
p7_oprofile_Destroy(P7_OPROFILE *om)
{
  if (om == NULL) return;

  if (om->tsc != NULL) free(om->tsc);
  if (om->rsc != NULL) 
    {
      if (om->rsc[0] != NULL) free(om->rsc[0]);
      free(om->rsc);
    }
  free(om);
}

static unsigned char
biased_charify(P7_OPROFILE *om, float sc)
{
  sc  = -1. * roundf(om->scale * sc);          	/* ugh. sc is now a nonnegative integer cost represented in a float... */
  sc += om->bias;		                /* and now we add an unsigned char to it... */
  return (sc > 255 ? 255 : (unsigned char) sc);	/* and now we cast it to an unsigned char cost. */
}
 

static unsigned char 
unbiased_charify(P7_OPROFILE *om, float sc)
{
  sc  = -1.0 * roundf(om->scale * sc);          
  return (sc > 255 ? 255 : (unsigned char) sc);	
}

/* Function:  p7_oprofile_Convert()
 * Synopsis:  Converts standard profile to an optimized one.
 * Incept:    SRE, Sun Dec  9 12:06:02 2007 [Janelia]
 */
int
p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om)
{
  int x;			/* counter over residues */
  int q;			/* q counts over total # of striped vectors, always 1.. p7O_NQ */
  int j;			/* j counts over position of vectors in some particular memory arrangement */
  int k;			/* the usual counter over model nodes 1..M */
  int kb;			/* base k for loading om's TSC quads */
  int t;			/* counter over transitions 0..7 = p7O_{BM, MM, IM, DM, MD, MI, II, DD} */
  int tg;			/* transition index in gm */
  int M  = gm->M;		/* length of the query */
  int nq = p7O_NQ(M);	        /* segment length; total # of striped quads */
  float max;
  int   ddtmp;
  unsigned char val[16];
  int   z;


  if (gm->abc_r->type != om->abc->type)  ESL_EXCEPTION(eslEINVAL, "alphabets of the two profiles don't match");
  if (nq > om->allocQ)                   ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  /* Determine the arithmetic basis for the reduced precision in unsigned ints */
  max = 0.;
  for (x = 0; x < gm->abc_r->K; x++)
    max = ESL_MAX(max, esl_vec_FMax(gm->rsc[x], (M+1)*2));

  /* We could implement a more smooth func for determining scale below,  
   * and squeaking maximal precision out; [xref J2/66] for considerations.
   */
  if (gm->M <=800)  om->scale = 3.0 / eslCONST_LOG2;    /* scores in units of third-bits for most models [xref J2/66]     */
  else              om->scale = 2.0 / eslCONST_LOG2;    /* scores in units of half-bits for large models                  */
  om->base  = 195;		                        /* bias in 0..255 scores in DP matrix: scores range -65..20 bits  */
  om->bias  = -1 * unbiased_charify(om, max);           /* match, insert emission costs all negative, by subtracting bias */

  /* match, insert emission costs: start at k=1.  */
  for (x = 0; x < gm->abc_r->Kp; x++)
    for (j = 0, k = 1, q = 0; q < nq; q++, k++)
      {
	for (z = 0; z < 16; z++) val[z] = ((k+ z*nq <= M) ? biased_charify(om, p7P_MSC(gm, k+z*nq, x)) : 255);
	om->rsc[x][j++] = _mm_load_si128((__m128i *) val);

	for (z = 0; z < 16; z++) val[z] = ((k+ z*nq <= M) ? biased_charify(om, p7P_ISC(gm, k+z*nq, x)) : 255);
	om->rsc[x][j++] = _mm_load_si128((__m128i *) val);
      }

  /* Transition costs, all but the DD's. */
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

	  for (z = 0; z < 16; z++) val[z] = ((kb+ z*nq < M) ? unbiased_charify(om, p7P_TSC(gm, kb+ z*nq, tg)) : 255);
	  om->tsc[j++] = _mm_load_si128((__m128i *) val);
	}
    }

  /* And finally the DD's, which are at the end of the optimized tsc vector; (j is already there) */
  for (k = 1, q = 0; q < nq; q++, k++)
    {
      for (z = 0; z < 16; z++) val[z] = ((k+ z*nq < M) ? unbiased_charify(om, p7P_TSC(gm, k+ z*nq, p7P_DD)) : 255);
      om->tsc[j++] = _mm_load_si128((__m128i *) val);
    }

  /* Specials. (These are actually in exactly the same order in om and
   *  gm, but we copy in general form anyway.)
   *  NN, CC, JJ are all hardcoded 0; part of precision-maximizing strategy [xref J2/66]
   */
  om->xsc[p7O_E][p7O_LOOP] = unbiased_charify(om, gm->xsc[p7P_E][p7P_LOOP]);  
  om->xsc[p7O_E][p7O_MOVE] = unbiased_charify(om, gm->xsc[p7P_E][p7P_MOVE]);
  om->xsc[p7O_N][p7O_LOOP] = 0;
  om->xsc[p7O_N][p7O_MOVE] = unbiased_charify(om, gm->xsc[p7P_N][p7P_MOVE]);
  om->xsc[p7O_C][p7O_LOOP] = 0;
  om->xsc[p7O_C][p7O_MOVE] = unbiased_charify(om, gm->xsc[p7P_C][p7P_MOVE]);
  om->xsc[p7O_J][p7O_LOOP] = 0;
  om->xsc[p7O_J][p7O_MOVE] = unbiased_charify(om, gm->xsc[p7P_J][p7P_MOVE]);

  /* Transition score bound for "lazy F" DD path evaluation (xref J2/52) */
  om->dd_bound = 255;	
  for (k = 2; k < M-1; k++) 
    {
      ddtmp  =  (int) unbiased_charify(om, p7P_TSC(gm, k+1, p7P_BM));
      ddtmp -=  (int) unbiased_charify(om, p7P_TSC(gm, k,   p7P_DD));
      ddtmp -=  (int) unbiased_charify(om, p7P_TSC(gm, k+1, p7P_DM));
      om->dd_bound = ESL_MIN(om->dd_bound, ddtmp);
    }

  om->M    = M;
  om->mode = gm->mode;
  return eslOK;
}


/* Function:  p7_oprofile_Dump()
 * Synopsis:  Dump internals of a <P7_OPROFILE>
 * Incept:    SRE, Sun Dec  9 14:58:52 2007 [Janelia]
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
  int k;			/* the usual counter over model nodes 1..M */
  int kb;
  int t;
  int i;			/* counter in elements within a vector */
  int M  = om->M;		/* length of the profile */
  int nq = p7O_NQ(M);	        /* segment length; total # of striped quads */
  union {
    __m128i v;
    unsigned char i[16];
  } tmp;

  /* Residue emissions
   */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 

      /* Header (rearranged column numbers, in the vectors)  */
      for (k =1, q = 0; q < nq; q++, k++)
	{
	  fprintf(fp, "[ ");
	  for (i = 0; i < p7O_QWIDTH; i++) 
	    if (k+i*nq <= M) fprintf(fp, "%4d ", k+i*nq);
	    else             fprintf(fp, "%4s ", "xx");
	  fprintf(fp, "]");
	}

      /* Match emission scores */
      fprintf(fp, "\nmat: ");
      for (j = 0, q = 0; q < nq; q++, j+=2)
	{
	  fprintf(fp, "[ ");
	  _mm_store_si128(&tmp.v, om->rsc[x][j]);
	  for (i = 0; i < p7O_QWIDTH; i++) fprintf(fp, "%4d ", tmp.i[i]);
	  fprintf(fp, "]");
	}

      fprintf(fp, "\nins: ");
      for (j = 1, q = 0; q < nq; q++, j+=2)
	{
	  fprintf(fp, "[ ");
	  _mm_store_si128(&tmp.v, om->rsc[x][j]);
	  for (i = 0; i < p7O_QWIDTH; i++) fprintf(fp, "%4d ", tmp.i[i]);
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
	  for (i = 0; i < p7O_QWIDTH; i++) 
	    if (kb+i*nq <= M) fprintf(fp, "%4d ", kb+i*nq);
	    else              fprintf(fp, "%4s ", "xx");
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n     ");	  
      for (q = 0; q < nq; q++)
	{
	  fprintf(fp, "[ ");
	  _mm_store_si128(&tmp.v, om->tsc[q*7 + t]);
	  for (i = 0; i < p7O_QWIDTH; i++) fprintf(fp, "%4d ", tmp.i[i]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");	  
    }

  /* DD transitions */
  fprintf(fp, "\ntDD: ");
  for (k =1, q = 0; q < nq; q++, k++)
    {
      fprintf(fp, "[ ");
      for (i = 0; i < p7O_QWIDTH; i++) 
	if (k+i*nq <= M) fprintf(fp, "%4d ", k+i*nq);
	else             fprintf(fp, "%4s ", "xx");
      fprintf(fp, "]");
    }
  fprintf(fp, "\n     ");	  
  for (j = nq*7, q = 0; q < nq; q++, j++)
    {
      fprintf(fp, "[ ");
      _mm_store_si128(&tmp.v, om->tsc[j]);
      for (i = 0; i < p7O_QWIDTH; i++) fprintf(fp, "%4d ", tmp.i[i]);
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");	  
  

  fprintf(fp, "E->C: %4d    E->J: %4d\n", om->xsc[p7O_E][p7O_MOVE], om->xsc[p7O_E][p7O_LOOP]);
  fprintf(fp, "N->B: %4d    N->N: %4d\n", om->xsc[p7O_N][p7O_MOVE], om->xsc[p7O_N][p7O_LOOP]);
  fprintf(fp, "J->B: %4d    J->J: %4d\n", om->xsc[p7O_J][p7O_MOVE], om->xsc[p7O_J][p7O_LOOP]);
  fprintf(fp, "C->T: %4d    C->C: %4d\n", om->xsc[p7O_C][p7O_MOVE], om->xsc[p7O_C][p7O_LOOP]);

  fprintf(fp, "bound: %4d\n",  om->dd_bound);
  fprintf(fp, "scale: %.2f\n", om->scale);
  fprintf(fp, "base:  %4d\n",  om->base);
  fprintf(fp, "bias:  %4d\n",  om->bias);
  fprintf(fp, "Q:     %d\n",   nq);  
  fprintf(fp, "M:     %d\n",   M);  


  return eslOK;
}


/*****************************************************************
 * 2. The P7_OMX structure: a dynamic programming matrix
 *****************************************************************/

/* Function:  p7_omx_Create()
 * Synopsis:  Create an optimized dynamic programming matrix.
 * Incept:    SRE, Sun Dec  9 15:35:49 2007 [Janelia]
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

  ESL_ALLOC(ox->dp,  sizeof(__m128i) * p7X_NSCELLS * nq);
  ox->M      = 0;
  ox->Q      = 0;
  ox->allocQ = nq;
#if p7_DEBUGGING
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
 * Incept:    SRE, Sun Dec  9 15:35:59 2007 [Janelia]
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
 * Incept:    SRE, Tue Dec 11 10:09:48 2007 [Janelia]
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
 * SRE, Sun Dec  9 15:36:04 2007 [Janelia]
 *
 * Dump DP matrix row <rowi> in <ox> to stream <fp> for diagnostics. 
 *
 * If <rowi> is 0, it prints a header first too.
 *            
 * The output format is coordinated with <p7_gmx_Dump()> to
 * facilitate comparison to a known answer.
 */
static int
omx_dump_row(P7_OMX *ox, int rowi, unsigned char xE, unsigned char xJ, unsigned char xB, unsigned char xC)
{
  int      q,z,k;
  unsigned char    *v;
  __m128i *dp = ox->dp;
  int      Q  = ox->Q;
  int      M  = ox->M;
  int      status;
  union { __m128i v; unsigned char i[16]; } tmp;

  ESL_ALLOC(v, sizeof(unsigned char) * ((Q*16)+1));
  v[0] = 0;

  /* Header (if we're on the 0th row)
   */
  if (rowi == 0)
    {
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M;  k++) fprintf(ox->dfp, "%3d ", k);
      fprintf(ox->dfp, "%3s %3s %3s %3s %3s\n", "E", "N", "J", "B", "C");
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M+5;  k++) fprintf(ox->dfp, "%3s ", "---");
      fprintf(ox->dfp, "\n");
    }

  /* Unpack, unstripe, then print M's. */
  for (q = 0; q < Q; q++) {
    _mm_store_si128(&tmp.v, MMX(q));
    for (z = 0; z < 16; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d M ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", v[k]);

  /* The specials */
  fprintf(ox->dfp, "%3d %3d %3d %3d %3d\n", xE, 0, xJ, xB, xC);

  /* Unpack, unstripe, then print I's. */
  for (q = 0; q < Q; q++) {
    _mm_store_si128(&tmp.v, IMX(q));
    for (z = 0; z < 16; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d I ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", v[k]);
  fprintf(ox->dfp, "\n");

  /* Unpack, unstripe, then print D's. */
  for (q = 0; q < Q; q++) {
    _mm_store_si128(&tmp.v, DMX(q));
    for (z = 0; z < 16; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d D ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", v[k]);
  fprintf(ox->dfp, "\n\n");

  free(v);
  return eslOK;

ERROR:
  free(v);
  return status;

}
#endif /*p7_DEBUGGING*/

/*****************************************************************
 * 3. The p7_ViterbiFilter() DP implementation.
 *****************************************************************/

/* Returns TRUE if any a[z] > b[z] .
 * This is an incantation. 
 * SSE provides no cmpgt_epu8 instruction!
 * Note that cmpeq_epi8 works fine for unsigned ints (there is no 
 * cmpeq_epu8 instruction either). 
 */
static int 
sse_any_gt_epu8(__m128i a, __m128i b)
{
  __m128i mask    = _mm_cmpeq_epi8(_mm_max_epu8(a,b), b); /* anywhere a>b, mask[z] = 0x0; elsewhere 0xff */
  int   maskbits  = _mm_movemask_epi8(_mm_xor_si128(mask,  _mm_cmpeq_epi8(mask, mask)));
  return maskbits != 0;
}

/* Returns maximum element \max_z a[z] in epu8 vector */
static unsigned char
sse_hmax_epu8(__m128i a)
{
  union { __m128i v; unsigned char i[16]; } tmp;
  
  tmp.v = _mm_max_epu8(a,     _mm_slli_si128(a,     1));
  tmp.v = _mm_max_epu8(tmp.v, _mm_slli_si128(tmp.v, 2));
  tmp.v = _mm_max_epu8(tmp.v, _mm_slli_si128(tmp.v, 4));
  tmp.v = _mm_max_epu8(tmp.v, _mm_slli_si128(tmp.v, 8));
  return tmp.i[15];
}



/* Function:  p7_ViterbiFilter()
 * Synopsis:  Calculates Viterbi score, insanely fast, in limited precision.
 * Incept:    SRE, Tue Nov 27 09:15:24 2007 [Janelia]
 *
 * Purpose:   Calculates an approximation of the Viterbi score for sequence
 *            <dsq> of length <L> residues, using optimized profile <om>,
 *            and a preallocated one-row DP matrix <ox>. Return the 
 *            Viterbi score (in nats) in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow.
 *            
 *            This is a striped SIMD Viterbi implementation in Intel
 *            SSE/SSE2 integer intrinsics \citep{Farrar07}, in reduced
 *            precision (unsigned chars).
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - DP matrix
 *            ret_sc  - RETURN: Viterbi score (in nats)          
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small, or if
 *            profile isn't in a local alignment mode. (Must be in local
 *            alignment mode because that's what helps us guarantee 
 *            limited dynamic range.)


 *
 * Xref:      [Farrar07] for ideas behind striped SIMD DP.
 *            J2/46-47 for layout of HMMER's striped SIMD DP.
 *            J2/50 for single row DP.
 *            J2/60 for reduced precision (epu8)
 *            J2/65 for initial benchmarking
 *            J2/66 for precision maximization
 */
int
p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register __m128i mpv, dpv, ipv;  /* previous row values */
  register __m128i sv;		   /* temp storage of 1 curr row value in progress */
  register __m128i dcv;		   /* delayed storage of D(i,q+1) */
  register __m128i xEv;		   /* E state: keeps max for Mk->E as we go */
  register __m128i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i Dmaxv;          /* keeps track of maximum D cell on row */
  __m128i  biasv;		   /* emission bias in a vector */
  unsigned char xE, xB, xC, xJ;	   /* special states */
  unsigned char Dmax;		   /* maximum D cell on row */
  int i;			   /* counter over sequence positions 1..L */
  int q;			   /* counter over quads 1..nq */
  int Q        = p7O_NQ(om->M);	   /* segment length: # of quads */
  __m128i *dp  = ox->dp;
  __m128i *rsc;			   /* will point at om->rsc[x] for residue x[i] */
  __m128i *tsc;			   /* will point into (and step thru) om->tsc   */

  /* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ) ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  if (om->mode != p7_LOCAL && om->mode != p7_UNILOCAL) ESL_EXCEPTION(eslEINVAL, "Fast filter only works for local alignment");
  ox->M = om->M;
  ox->Q = Q;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base.
   */
  biasv = _mm_set1_epi8((char) om->bias); /* yes, you can set1() an unsigned char vector this way */
  for (q = 0; q < Q; q++)
    MMX(q) = IMX(q) = DMX(q) = _mm_setzero_si128();
  xB   = om->base - om->xsc[p7O_N][p7O_MOVE]; /* remember, all values are costs to be subtracted. */
  xJ   = 0;
  xC   = 0;
  xE   = 0;

#if p7_DEBUGGING
  if (ox->debugging) omx_dump_row(ox, 0, xE, xJ, xB, xC);   
#endif

  for (i = 1; i <= L; i++)
    {
      rsc   = om->rsc[dsq[i]];
      tsc   = om->tsc;
      dcv   = _mm_setzero_si128();      /* "-infinity" */
      xEv   = _mm_setzero_si128();     
      Dmaxv = _mm_setzero_si128();     
      xBv   = _mm_set1_epi8((char) xB);

      /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically, which is our -infinity.
       */
      mpv = MMX(Q-1);  mpv = _mm_slli_si128(mpv, 1);  
      dpv = DMX(Q-1);  dpv = _mm_slli_si128(dpv, 1);  
      ipv = IMX(Q-1);  ipv = _mm_slli_si128(ipv, 1);  

      for (q = 0; q < Q; q++)
	{
	  /* Calculate new MMX(i,q); don't store it yet, hold it in sv. */
	  sv   =                   _mm_subs_epu8(xBv, *tsc);  tsc++;
	  sv   = _mm_max_epu8 (sv, _mm_subs_epu8(mpv, *tsc)); tsc++;
	  sv   = _mm_max_epu8 (sv, _mm_subs_epu8(ipv, *tsc)); tsc++;
	  sv   = _mm_max_epu8 (sv, _mm_subs_epu8(dpv, *tsc)); tsc++;
	  sv   = _mm_adds_epu8(sv, biasv);     
	  sv   = _mm_subs_epu8(sv, *rsc);                     rsc++;
	  xEv  = _mm_max_epu8(xEv, sv);
	  
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
	  dcv   = _mm_subs_epu8(sv, *tsc);  tsc++;
	  Dmaxv = _mm_max_epu8(dcv, Dmaxv);

	  /* Calculate and store I(i,q) */
	  sv     =                   _mm_subs_epu8(mpv, *tsc);  tsc++;
	  sv     = _mm_max_epu8 (sv, _mm_subs_epu8(ipv, *tsc)); tsc++;
	  sv     = _mm_adds_epu8(sv, biasv);
	  IMX(q) = _mm_subs_epu8(sv, *rsc);                     rsc++;
	}	  

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      xE = sse_hmax_epu8(xEv);
      xC = ESL_MAX(xC - om->xsc[p7O_C][p7O_LOOP],        xE - om->xsc[p7O_E][p7O_MOVE]);
      xJ = ESL_MAX(xJ - om->xsc[p7O_J][p7O_LOOP],        xE - om->xsc[p7O_E][p7O_LOOP]);
      xB = ESL_MAX(xJ - om->xsc[p7O_J][p7O_MOVE],  om->base - om->xsc[p7O_N][p7O_MOVE]);
      /* and now xB will carry over into next i, and xC carries over after i=L */

      /* Finally the "lazy F" loop (sensu [Farrar07]). We can often
       * prove that we don't need to evaluate any D->D paths at all.
       *
       * The observation is that if we can show that on the next row,
       * B->M(i+1,k) paths always dominate M->D->...->D->M(i+1,k) paths
       * for all k, then we don't need any D->D calculations.
       * 
       * The test condition is:
       *      max_k D(i,k) + max_k ( TDD(k-2) + TDM(k-1) - TBM(k) ) < xB(i)
       * So:
       *   max_k (TDD(k-2) + TDM(k-1) - TBM(k)) is precalc'ed in om->dd_bound;
       *   max_k D(i,k) is why we tracked Dmaxv;
       *   xB(i) was just calculated above.
       */
      Dmax = sse_hmax_epu8(Dmaxv);
      if ((int) Dmax + om->dd_bound > (int) xB) 
	{
	  /* Now we're obligated to do at least one complete DD path to be sure. */
	  /* dcv has carried through from end of q loop above */
	  dcv = _mm_slli_si128(dcv, 1);
	  tsc = om->tsc + 7*Q;	/* set tsc to start of the DD's */
	  for (q = 0; q < Q; q++) 
	    {
	      DMX(q) = _mm_max_epu8(dcv, DMX(q));	
	      dcv    = _mm_subs_epu8(DMX(q), *tsc); tsc++;
	    }

	  /* We may have to do up to three more passes; the check
	   * is for whether crossing a segment boundary can improve
	   * our score. 
	   */
	  do {
	    dcv = _mm_slli_si128(dcv, 1);
	    tsc = om->tsc + 7*Q;	/* set tsc to start of the DD's */
	    for (q = 0; q < Q; q++) 
	      {
		if (! sse_any_gt_epu8(dcv, DMX(q))) break;
		DMX(q) = _mm_max_epu8(dcv, DMX(q));	
		dcv    = _mm_subs_epu8(DMX(q), *tsc);   tsc++;
	      }	    
	  } while (q == Q);
	}
      else  /* not calculating DD? then just store the last M->D vector calc'ed.*/
	DMX(0) = _mm_slli_si128(dcv, 1);

	  
#if p7_DEBUGGING
      if (ox->debugging) omx_dump_row(ox, i, xE, xJ, xB, xC);   
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T, and add our missing precision on the NN,CC,JJ back */
  *ret_sc = ((float) (xC - om->xsc[p7O_C][p7O_MOVE]) - (float) om->base);
  *ret_sc /= om->scale;
  if      (om->mode == p7_UNILOCAL) *ret_sc -= 2.0; /* that's ~ L \log \frac{L}{L+2}, for our NN,CC,JJ */
  else if (om->mode == p7_LOCAL)    *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
  return eslOK;
}

/*****************************************************************
 * 3. Private (static) functions for debugging.
 *****************************************************************/

/* Round all the scores in a generic P7_PROFILE in exactly the same
 * way that the scores in an sse8 P7_OPROFILE were rounded.  Then the
 * two profiles should give identical internal scores.
 */
static void
round_profile(P7_OPROFILE *om, P7_PROFILE *gm)
{
  int k;
  int x;

  /* Transitions */
  /* <= -eslINFINITY test is used solely to silence compiler. really testing == -eslINFINITY */
  for (x = 0; x < gm->M*p7P_NTRANS; x++)
      gm->tsc[x] = (gm->tsc[x] <= -eslINFINITY) ? -255 : roundf(om->scale * gm->tsc[x]);
  
  /* Emissions */
  for (x = 0; x < gm->abc_r->Kp; x++)
    for (k = 0; k <= p7P_NR*gm->M; k++)
      gm->rsc[x][k] = (gm->rsc[x][k] <= -eslINFINITY) ? -255 : roundf(om->scale * gm->rsc[x][k]);

  /* Specials */
  for (k = 0; k < p7P_NXSTATES; k++)
    for (x = 0; x < p7P_NXTRANS; x++)
      gm->xsc[k][x] = (gm->xsc[k][x] <= -eslINFINITY) ? -255 : roundf(om->scale * gm->xsc[k][x]);
}



/*****************************************************************
 * Exam drivers.
 *****************************************************************/

#ifdef p7IMPL_SSE8_EXAM
/* Examining internals of a converted optimized profile
   gcc -msse2 -g -Wall -I. -L. -I../easel -L../easel -o exam-impl-sse8 -Dp7IMPL_SSE8_EXAM impl_sse8.c -lhmmer -leasel -lm
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
#endif /*p7IMPL_SSE8_EXAM*/



#ifdef p7IMPL_SSE8_EXAM2
/* Comparing internals of small DP matrices
   gcc -msse2 -g -Wall -I. -L. -I../easel -L../easel -o exam -Dp7IMPL_SSE8_EXAM2 impl_sse8.c -lhmmer -leasel -lm
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
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "round generic profile so scores match (debug)", 0 },
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

  if (esl_opt_GetBoolean(go, "-x")) round_profile(om, gm);

 
  ox = p7_omx_Create(gm->M);
  gx = p7_gmx_Create(gm->M, sq->n);
  p7_omx_SetDumpMode(stdout, ox);

  p7_ViterbiFilter(sq->dsq, sq->n, om, ox, &sc1); 

  p7_GViterbi     (sq->dsq, sq->n, gm, gx, &sc2);
  p7_gmx_Dump(stdout, gx);
  if (esl_opt_GetBoolean(go, "-x")) sc2 /= om->scale;

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
#endif /*p7IMPL_SSE8_EXAM2*/


/*****************************************************************
 * Benchmark driver.
 *****************************************************************/
#ifdef p7IMPL_SSE8_BENCHMARK
/* gcc -o benchmark-sse8 -g -O3 -msse2 -I. -L. -I../easel -L../easel -Dp7IMPL_SSE8_BENCHMARK impl_sse8.c -lhmmer -leasel -lm
 * icc -o benchmark-sse8 -O3 -static -I. -L. -I../easel -L../easel -Dp7IMPL_SSE8_BENCHMARK impl_sse8.c -lhmmer -leasel -lm 
 * icc -o benchmark-sse8 -g -Wall -static -I. -L. -I../easel -L../easel -Dp7IMPL_SSE8_BENCHMARK impl_sse8.c -lhmmer -leasel -lm 
 * 
 *   ./benchmark-sse8 <hmmfile>         runs benchmark
 *   ./benchmark-sse8 -b <hmmfile>      gets baseline time to subtract: just random seq generation
 *   ./benchmark-sse8 -cx <hmmfile>     compare scores of SSE to generic impl rounded the same way
 *   ./benchmark-sse8 -c <hmmfile>      compare scores of SSE to generic impl in full precision 
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
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "round generic profile, make scores match (debug)", 0 },
  { "-b",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "baseline timing: don't run DP at all",             0 },
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

  if (esl_opt_GetBoolean(go, "-x")) round_profile(om, gm);

  ox = p7_omx_Create(gm->M);
  gx = p7_gmx_Create(gm->M, L);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rnd_xfIID(r, bg->f, abc->K, L, dsq);

      if (! esl_opt_GetBoolean(go, "-b")) {
	p7_ViterbiFilter(dsq, L, om, ox, &sc1);   

	if (esl_opt_GetBoolean(go, "-c")) {
	  p7_GViterbi     (dsq, L, gm, gx, &sc2); 
	  if (esl_opt_GetBoolean(go, "-x")) sc2 /= om->scale;
	  
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
