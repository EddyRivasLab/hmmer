/* SSE implementation of Viterbi filter.
 * 
 * A filter module provides a standard API to the p7_ViterbiFilter()
 * call. The API requires that the module implement support for two
 * objects, a P7_OPROFILE optimized score profile and a P7_OMX
 * optimized dynamic programming matrix.
 * 
 *   1. The P7_OPROFILE structure: a score profile.
 *   2. The P7_OMX structure: a dynamic programming matrix.
 *   
 *   
 *-----------------------------------------------------------------------
 * 
 * Note on pointer alignment: many SSE calls require __m128 values to
 * to be aligned on 16-byte boundaries. I believe (but am not sure)
 * that ISO C malloc() must return a pointer properly aligned for the
 * largest legal type, and since __m128t is a legal type, malloc()
 * should return us pointers allocated on 16-byte boundaries. (This is
 * certainly true on my development machine.) Moreover, because we
 * malloc() vectors for SSE ops as arrays of __m128 (rather than, say,
 * loading our __m128 registers from float arrays that might be
 * unaligned), malloc() should *definitely* be returning us pointers
 * aligned on 16-byte boundaries. Nonetheless, to be Extra Double
 * Sure, we manually align all our malloc()'ed pointers anyway: you'll
 * see the idiom of
 *     p_mem = malloc();
 *     p     = (MYTYPE *) (((size_t) p + 15) & (~0xf))
 *     ...
 *     free(p_mem);
 * even though this is almost certainly overkill.
 *-----------------------------------------------------------------------
 * 
 * SRE, Sun Nov 25 11:26:48 2007 [Casa de Gatos]
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
#include "impl_sse.h"



/*****************************************************************
 * 1. The P7_OPROFILE structure: a score profile.
 *****************************************************************/

/* Function:  p7_oprofile_Create()
 * Synopsis:  Allocate an optimized profile structure.
 * Incept:    SRE, Sun Nov 25 12:03:19 2007 [Casa de Gatos]
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
  om->tsc     = NULL;
  om->rsc     = NULL;

  /* level 1 */
  ESL_ALLOC(om->tsc, sizeof(__m128) * nq  * p7O_NTRANS);    
  ESL_ALLOC(om->rsc, sizeof(__m128 *) * abc->Kp); 
  om->rsc[0] = NULL;

  /* level 2 */
  ESL_ALLOC(om->rsc[0], sizeof(__m128) * nq  * p7O_NR * abc->Kp);                     
  for (x = 1; x < abc->Kp; x++)
    om->rsc[x] = om->rsc[0] + (x * nq * p7O_NR);

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
 * Incept:    SRE, Sun Nov 25 12:22:21 2007 [Casa de Gatos]
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


/* Function:  p7_oprofile_Convert()
 * Synopsis:  Converts standard profile to an optimized one.
 * Incept:    SRE, Mon Nov 26 07:38:57 2007 [Janelia]
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
	om->rsc[x][j++] = _mm_setr_ps((k+0*nq <= M) ? p7P_MSC(gm, k,      x) : -eslINFINITY,
				      (k+1*nq <= M) ? p7P_MSC(gm, k+nq,   x) : -eslINFINITY,
				      (k+2*nq <= M) ? p7P_MSC(gm, k+2*nq, x) : -eslINFINITY,
				      (k+3*nq <= M) ? p7P_MSC(gm, k+3*nq, x) : -eslINFINITY);

	om->rsc[x][j++] = _mm_setr_ps((k+0*nq <= M) ? p7P_ISC(gm, k,      x) : -eslINFINITY,
				      (k+1*nq <= M) ? p7P_ISC(gm, k+nq,   x) : -eslINFINITY,
				      (k+2*nq <= M) ? p7P_ISC(gm, k+2*nq, x) : -eslINFINITY,
				      (k+3*nq <= M) ? p7P_ISC(gm, k+3*nq, x) : -eslINFINITY);
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

	  om->tsc[j++] = _mm_setr_ps((kb      < M) ? p7P_TSC(gm, kb,      tg) : -eslINFINITY,
	                             (kb+1*nq < M) ? p7P_TSC(gm, kb+nq,   tg) : -eslINFINITY,
				     (kb+2*nq < M) ? p7P_TSC(gm, kb+2*nq, tg) : -eslINFINITY,
				     (kb+3*nq < M) ? p7P_TSC(gm, kb+3*nq, tg) : -eslINFINITY);
	}
    }

  /* And finally the DD's, which are at the end of the optimized tsc vector; (j is already there) */
  for (k = 1, q = 0; q < nq; q++, k++)
    om->tsc[j++] = _mm_setr_ps((k+0*nq < M) ? p7P_TSC(gm, k,      p7P_DD) : -eslINFINITY,
			       (k+1*nq < M) ? p7P_TSC(gm, k+nq,   p7P_DD) : -eslINFINITY,
			       (k+2*nq < M) ? p7P_TSC(gm, k+2*nq, p7P_DD) : -eslINFINITY,
			       (k+3*nq < M) ? p7P_TSC(gm, k+3*nq, p7P_DD) : -eslINFINITY);

  /* Specials. (These are actually in exactly the same order in om and
   *  gm, but we copy in general form anyway.)
   */
  om->xsc[p7O_E][p7O_LOOP] = gm->xsc[p7P_E][p7P_LOOP];  
  om->xsc[p7O_E][p7O_MOVE] = gm->xsc[p7P_E][p7P_MOVE];
  om->xsc[p7O_N][p7O_LOOP] = gm->xsc[p7P_N][p7P_LOOP];
  om->xsc[p7O_N][p7O_MOVE] = gm->xsc[p7P_N][p7P_MOVE];
  om->xsc[p7O_C][p7O_LOOP] = gm->xsc[p7P_C][p7P_LOOP];
  om->xsc[p7O_C][p7O_MOVE] = gm->xsc[p7P_C][p7P_MOVE];
  om->xsc[p7O_J][p7O_LOOP] = gm->xsc[p7P_J][p7P_LOOP];
  om->xsc[p7O_J][p7O_MOVE] = gm->xsc[p7P_J][p7P_MOVE];

  om->M    = M;
  om->mode = gm->mode;
  return eslOK;
}


/* Function:  p7_oprofile_Dump()
 * Synopsis:  Dump the internals of a P7_OPROFILE
 * Incept:    SRE, Mon Nov 26 09:09:30 2007 [Janelia]
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
  int M  = om->M;		/* length of the profile */
  int nq = p7O_NQ(M);	        /* segment length; total # of striped quads */
  float tmp[4];

  /* Residue emissions
   */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 
      for (k =1, q = 0; q < nq; q++, k++)
	{
	  fprintf(fp, "[ %5d %5d %5d ", k, k+nq, k+2*nq);
	  if (k+3*nq <= M) fprintf(fp, "%5d ] ", k+3*nq);
	  else             fprintf(fp, "%5s ] ", "xx");
	}

      fprintf(fp, "\nmat: ");
      for (j = 0, q = 0; q < nq; q++, j+=2)
	{
	  _mm_store_ps(tmp, om->rsc[x][j]);
	  fprintf(fp, "[ %5.2f %5.2f %5.2f %5.2f ] ", tmp[0], tmp[1], tmp[2], tmp[3]);
	}

      fprintf(fp, "\nins: ");
      for (j = 1, q = 0; q < nq; q++, j+=2)
	{
	  _mm_store_ps(tmp, om->rsc[x][j]);
	  fprintf(fp, "[ %5.2f %5.2f %5.2f %5.2f ] ", tmp[0], tmp[1], tmp[2], tmp[3]);
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
	  fprintf(fp, "[ %5d %5d %5d ", kb, kb+nq, kb+2*nq);
	  if (kb+3*nq <= M) fprintf(fp, "%5d ] ", kb+3*nq);
	  else              fprintf(fp, "%5s ] ", "xx");
	}
      fprintf(fp, "\n     ");	  
      for (q = 0; q < nq; q++)
	{
	  _mm_store_ps(tmp, om->tsc[q*7 + t]);
	  fprintf(fp, "[ %5.2f %5.2f %5.2f %5.2f ] ", tmp[0], tmp[1], tmp[2], tmp[3]);
	}
      fprintf(fp, "\n");	  
    }

  /* DD transitions */
  fprintf(fp, "\ntDD: ");
  for (k =1, q = 0; q < nq; q++, k++)
    {
      fprintf(fp, "[ %5d %5d %5d ", k, k+nq, k+2*nq);
      if (k+3*nq <= M) fprintf(fp, "%5d ] ", k+3*nq);
      else             fprintf(fp, "%5s ] ", "xx");
    }
  fprintf(fp, "\n     ");	  
  for (j = nq*7, q = 0; q < nq; q++, j++)
    {
      _mm_store_ps(tmp, om->tsc[j]);
      fprintf(fp, "[ %5.2f %5.2f %5.2f %5.2f ] ", tmp[0], tmp[1], tmp[2], tmp[3]);
    }
  fprintf(fp, "\n");	  
  
  return eslOK;
}


/*****************************************************************
 * 2. The P7_OMX structure: a dynamic programming matrix
 *****************************************************************/

/* Function:  p7_omx_Create()
 * Synopsis:  Create an optimized dynamic programming matrix.
 * Incept:    SRE, Tue Nov 27 08:48:20 2007 [Janelia]
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
  return ox;

 ERROR:
  p7_omx_Destroy(ox);
  return NULL;
}


/* Function:  p7_omx_Destroy()
 * Synopsis:  Frees an optimized DP matrix.
 * Incept:    SRE, Tue Nov 27 09:11:42 2007 [Janelia]
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


/* Function:  p7_omx_Dump()
 * Synopsis:  Dump an optimized DP matrix for debugging.
 * Incept:    SRE, Wed Nov 28 09:00:38 2007 [Janelia]
 *
 * Purpose:   Dump DP matrix <ox> to stream <fp> for diagnostics.  The
 *            matrix contains values for one current row <rowi>.
 *
 *            If <rowi> is 0, it prints a header first too.
 * 
 *            Because the SSE implementation uses a one-row matrix,
 *            this routine is generally going to be called (during
 *            code development) from within <p7_ViterbiFilter()>,
 *            after each new row <rowi> is completed.
 *            
 *            The output format is coordinated with <p7_gmx_Dump()> to
 *            facilitate comparison to a known answer.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_omx_Dump(FILE *ofp, P7_OMX *ox, int rowi)
{
  int     q,z,k;
  float  *v;
  float   tmp[4];
  __m128 *dp = ox->dp;
  int     Q  = ox->Q;
  int     M  = ox->M;
  int     status;

  ESL_ALLOC(v, sizeof(float) * (M+1));
  v[0] = -eslINFINITY;

  if (rowi == 0)
    {
      fprintf(ofp, "     ");
      for (k = 0; k <= M;  k++) fprintf(ofp, "%8d ", k);
      fprintf(ofp, "\n");
    }

  /* Unpack, unstripe, then print M's. */
  for (q = 0; q < Q; q++) {
    _mm_store_ps(tmp, MMX(q));
    for (z = 0; z < 4; z++) v[q+Q*z+1] = tmp[z];
  }
  fprintf(ofp, "%3d M ", rowi);
  for (k = 0; k <= M; k++) fprintf(ofp, "%8.4f ", v[k]);
  fprintf(ofp, "\n");

  /* Unpack, unstripe, then print I's. */
  for (q = 0; q < Q; q++) {
    _mm_store_ps(tmp, IMX(q));
    for (z = 0; z < 4; z++) v[q+Q*z+1] = tmp[z];
  }
  fprintf(ofp, "%3d I ", rowi);
  for (k = 0; k <= M; k++) fprintf(ofp, "%8.4f ", v[k]);
  fprintf(ofp, "\n");

  /* Unpack, unstripe, then print D's. */
  for (q = 0; q < Q; q++) {
    _mm_store_ps(tmp, DMX(q));
    for (z = 0; z < 4; z++) v[q+Q*z+1] = tmp[z];
  }
  fprintf(ofp, "%3d D ", rowi);
  for (k = 0; k <= M; k++) fprintf(ofp, "%8.4f ", v[k]);
  fprintf(ofp, "\n\n");

  free(v);
  return eslOK;

ERROR:
  return status;

}




/*****************************************************************
 * 3. The p7_ViterbiFilter() DP implementation.
 *****************************************************************/

/* Return TRUE if any a[i] != b[i]; from Apple's Altivec/SSE migration guide */
static int 
sse_any_neq(__m128 a, __m128 b)
{
  __m128 mask     = _mm_cmpneq_ps(a,b);
  int   maskbits = _mm_movemask_ps( mask );
  return maskbits != 0;
}


/* Function:  p7_ViterbiFilter()
 * Synopsis:  Calculates Viterbi score, insanely fast.
 * Incept:    SRE, Tue Nov 27 09:15:24 2007 [Janelia]
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
p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register __m128 mpv, dpv, ipv;   /* previous row values */
  register __m128 sv;		   /* temp storage of 1 curr row value in progress */
  register __m128 dcv;		   /* delayed storage of D(i,q+1) */
  register __m128 xEv;		   /* E state: keeps max for Mk->E as we go */
  register __m128 xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  __m128  infv;			   /* -eslINFINITY in a vector */
  float    xE, xB, xC, xJ;	   /* special states */
  int i;			   /* counter over sequence positions 1..L */
  int q;			   /* counter over quads 1..nq */
  int z;
  int nq      = p7O_NQ(om->M);	   /* segment length: # of quads */
  __m128 *dp  = ox->dp;
  __m128 *rsc;			   /* will point at om->rsc[x] for residue x[i] */
  __m128 *tsc;			   /* will point into (and step thru) om->tsc   */

  /* Check that the DP matrix is ok for us. */
  if (nq > ox->allocQ) ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M = om->M;
  ox->Q = nq;

  /* Initialization.
   */
  infv = _mm_set1_ps(-eslINFINITY);
  for (q = 0; q < nq; q++)
    MMX(q) = IMX(q) = DMX(q) = infv;
  xB   = om->xsc[p7O_N][p7O_MOVE];
  xJ   = -eslINFINITY;
  xC   = -eslINFINITY;
  /* p7_omx_Dump(stdout, ox, i); */

  for (i = 1; i <= L; i++)
    {
      rsc  = om->rsc[dsq[i]];
      tsc  = om->tsc;
      dcv  = infv;
      xEv  = infv;
      xBv  = _mm_set1_ps(xB);

      /* Right shifts by 4 bytes. 4,8,12,x becomes x,4,8,12. 
       */
      mpv = MMX(nq-1);  mpv = _mm_shuffle_ps(mpv, mpv, _MM_SHUFFLE(2, 1, 0, 0));   mpv = _mm_move_ss(mpv, infv);
      dpv = DMX(nq-1);  dpv = _mm_shuffle_ps(dpv, dpv, _MM_SHUFFLE(2, 1, 0, 0));   dpv = _mm_move_ss(dpv, infv);
      ipv = IMX(nq-1);  ipv = _mm_shuffle_ps(ipv, ipv, _MM_SHUFFLE(2, 1, 0, 0));   ipv = _mm_move_ss(ipv, infv);
      
      for (q = 0; q < nq; q++)
	{
	  /* Calculate new MMX(i,q); don't store it yet, hold it in sv. */
	  sv   =                _mm_add_ps(xBv, *tsc);  tsc++;
	  sv   = _mm_max_ps(sv, _mm_add_ps(mpv, *tsc)); tsc++;
	  sv   = _mm_max_ps(sv, _mm_add_ps(ipv, *tsc)); tsc++;
	  sv   = _mm_max_ps(sv, _mm_add_ps(dpv, *tsc)); tsc++;
	  sv   = _mm_add_ps(sv, *rsc);                  rsc++;
	  xEv  = _mm_max_ps(xEv, sv);
	  
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
	  dcv = _mm_add_ps(sv, *tsc); tsc++;

	  /* Calculate and store I(i,q) */
	  sv     =                _mm_add_ps(mpv, *tsc);  tsc++;
	  sv     = _mm_max_ps(sv, _mm_add_ps(ipv, *tsc)); tsc++;
	  IMX(q) = _mm_add_ps(sv, *rsc);                  rsc++;
	}	  

      /* Now the "lazy F" loop (sensu [Farrar07]), evaluating the D->D
       * path.  The dependencies within the current row make this a
       * pain to parallelize. The hope is that lazy F only requires one
       * pass; worst case, it requires 4, amounting to serialization.
       */
      dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
      dcv = _mm_move_ss(dcv, infv);
      for (z = 0; z < 4; z++)
      /*      while (sse_any_neq(dcv, DMX(0))) */
	{
	  tsc = om->tsc + 7*nq;	/* set tsc to start of the DD's */
	  for (q = 0; q < nq; q++)
	    {
	      DMX(q) = _mm_max_ps(dcv, DMX(q));	/* delayed max! */
	      dcv    = _mm_add_ps(DMX(q), *tsc);   tsc++;
	    }
	  dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
	  dcv = _mm_move_ss(dcv, infv);
	}

      
      /* Now the "special" states */
      /* The following incantation takes the max of xEv's elements  */
      xEv = _mm_max_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(0, 3, 2, 1)));
      xEv = _mm_max_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(1, 0, 3, 2)));
      _mm_store_ss(&xE, xEv);

      xJ = ESL_MAX(xJ + om->xsc[p7O_J][p7O_LOOP],  xE + om->xsc[p7O_E][p7O_LOOP]);
      xC = ESL_MAX(xC + om->xsc[p7O_C][p7O_LOOP],  xE + om->xsc[p7O_E][p7O_MOVE]);
      xB = ESL_MAX(xJ + om->xsc[p7O_J][p7O_MOVE],
		   i * om->xsc[p7O_N][p7O_LOOP] + om->xsc[p7O_N][p7O_MOVE]);
      /* and now xB carries over into next i, and xC carries over after i=L */

      /* p7_omx_Dump(stdout, ox, i); */
    } /* end loop over sequence residues 1..L */

  /* finally C->T */
  *ret_sc = xC + om->xsc[p7O_C][p7O_MOVE];
  return eslOK;
}


/*****************************************************************
 * Exam drivers.
 *****************************************************************/

#ifdef p7IMPL_SSE_EXAM
/* Examining internals of a converted optimized profile
   gcc -msse2 -g -Wall -I. -L. -I../easel -L../easel -o exam -Dp7IMPL_SSE_EXAM impl_sse.c -lhmmer -leasel -lm
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
#endif /*p7IMPL_SSE_EXAM*/


#ifdef p7IMPL_SSE_EXAM2
/* Comparing internals of small DP matrices
   gcc -msse2 -g -Wall -I. -L. -I../easel -L../easel -o exam -Dp7IMPL_SSE_EXAM2 impl_sse.c -lhmmer -leasel -lm
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
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
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

  ox = p7_omx_Create(gm->M);
  gx = p7_gmx_Create(gm->M, sq->n);

  p7_ViterbiFilter(sq->dsq, sq->n, om, ox, &sc1); 

  p7_GViterbi     (sq->dsq, sq->n, gm, gx, &sc2);
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
#endif /*p7IMPL_SSE_EXAM2*/


/*****************************************************************
 * Benchmark driver.
 *****************************************************************/
#ifdef p7IMPL_SSE_BENCHMARK
/* gcc -o benchmark-sse -g -O3 -msse2 -I. -L. -I../easel -L../easel -Dp7IMPL_SSE_BENCHMARK impl_sse.c -lhmmer -leasel -lm
 * ./benchmark-sse <hmmfile>
 */
/* As of Tue Jul 17 13:14:43 2007
 * 61 Mc/s for Viterbi, 8.6 Mc/s for Forward.
 * (gcc -g -O2, 3.2GHz Xeon, N=50K, L=400, M=72 RRM_1 model)
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
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose: show individual scores",             0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
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
      
      /* p7_ViterbiFilter(dsq, L, om, ox, &sc1);  */
      /* p7_GViterbi     (dsq, L, gm, gx, &sc2); */

      /* printf("sc1= %.4f  sc2 = %.4f\n", sc1, sc2); */

    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

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
#endif /*p7DP_GENERIC_BENCHMARK*/

