/*  SSE versions of Routines for the P7_OPROFILE structure:  
 * a search profile in an optimized implementation.
 * 
 * Contents:
 *   1. The P7_OPROFILE object: allocation, initialization, destruction.
 *   2. Conversion from generic P7_PROFILE to optimized P7_OPROFILE
 *   3. Conversion from optimized P7_OPROFILE to compact score arrays
 *   4. Debugging and development utilities.
 *   5. Benchmark driver.
 *   6. Example.
 *   7. Copyright and license information.
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* roundf() */

#if p7_CPU_ARCH == intel
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
#endif /* intel arch */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "base/p7_bg.h"
#include "base/p7_hmm.h"
#include "base/p7_profile.h"

#include "build/modelsample.h"
#include "search/modelconfig.h"

#include "dp_vector/p7_oprofile.h"

/*****************************************************************
 * 1. The P7_OPROFILE structure: a score profile.
 *****************************************************************/

/* Function:  p7_oprofile_Create_sse()
 * Synopsis:  Allocate an optimized profile structure, configured to use the SSE vector ISA..
 * Incept:    SRE, Sun Nov 25 12:03:19 2007 [Casa de Gatos]
 *
 * Purpose:   Allocate for profiles of up to <allocM> nodes for digital alphabet <abc>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_OPROFILE *
p7_oprofile_Create_sse(int allocM, const ESL_ALPHABET *abc)
{
#ifdef HAVE_SSE2

  int          status;
  P7_OPROFILE *om  = NULL;
  int          nqb = P7_NVB(allocM); /* # of uchar vectors needed for query */
  int          nqw = P7_NVW(allocM); /* # of sword vectors needed for query */
  int          nqf = P7_NVF(allocM); /* # of float vectors needed for query */
  int          nqs = nqb + p7O_EXTRA_SB;


  int          x;

  /* level 0 */
  ESL_ALLOC(om, sizeof(P7_OPROFILE));
  om->simd = SSE;
  om->rbv_mem   = NULL; 
  om->rbv       = NULL;
  om->rfv_mem   = NULL;
  om->tfv_mem   = NULL;
  om->rfv       = NULL;
  om->tfv       = NULL;
  om->sbv_mem   = NULL;
  om->sbv       = NULL;
  om->rwv_mem   = NULL;
  om->twv_mem   = NULL;
  om->rwv       = NULL;
  om->twv       = NULL;

  om->is_shadow = FALSE;

  om->name    = NULL;
  om->acc     = NULL;
  om->desc    = NULL;

  om->rf        = NULL;
  om->mm        = NULL;
  om->cs        = NULL;
  om->consensus = NULL;

  /* level 1 */
  ESL_ALLOC(om->rbv_mem, sizeof(__m128i) * nqb  * abc->Kp          +15); /* +15 is for manual 16-byte alignment */
    ESL_ALLOC(om->sbv_mem, sizeof(__m128i) * nqs  * abc->Kp          +15); 
  ESL_ALLOC(om->rwv_mem, sizeof(__m128i) * nqw  * abc->Kp          +15);                     
  ESL_ALLOC(om->twv_mem, sizeof(__m128i) * nqw  * p7O_NTRANS       +15);   
  ESL_ALLOC(om->rfv_mem, sizeof(__m128)  * nqf  * abc->Kp          +15);                     
  ESL_ALLOC(om->tfv_mem, sizeof(__m128)  * nqf  * p7O_NTRANS       +15);    
  ESL_ALLOC(om->rbv, sizeof(__m128i *) * abc->Kp); 
  ESL_ALLOC(om->sbv, sizeof(__m128i *) * abc->Kp); 
  ESL_ALLOC(om->rwv, sizeof(__m128i *) * abc->Kp); 
  ESL_ALLOC(om->rfv, sizeof(__m128  *) * abc->Kp); 

  /* align vector memory on 16-byte boundaries */
  om->rbv[0] = (__m128i *) (((unsigned long int) om->rbv_mem + 15) & (~0xf));
  om->sbv[0] = (__m128i *) (((unsigned long int) om->sbv_mem + 15) & (~0xf));
  om->rwv[0] = (__m128i *) (((unsigned long int) om->rwv_mem + 15) & (~0xf));
  om->twv    = (__m128i *) (((unsigned long int) om->twv_mem + 15) & (~0xf));
  om->rfv[0] = (__m128  *) (((unsigned long int) om->rfv_mem + 15) & (~0xf));
  om->tfv    = (__m128  *) (((unsigned long int) om->tfv_mem + 15) & (~0xf));

  /* set the rest of the row pointers for match emissions */
  for (x = 1; x < abc->Kp; x++) {
    om->rbv[x] = om->rbv[0] + (x * nqb);
    om->sbv[x] = om->sbv[0] + (x * nqs);
    om->rwv[x] = om->rwv[0] + (x * nqw);
    om->rfv[x] = om->rfv[0] + (x * nqf);

  }

  om->allocQ16  = nqb;
  om->allocQ8   = nqw;
  om->allocQ4   = nqf;

  /* Remaining initializations */
  om->tbm_b     = 0;
  om->tec_b     = 0;
  om->tjb_b     = 0;
  om->scale_b   = 0.0f;
  om->base_b    = 0;
  om->bias_b    = 0;

  om->scale_w      = 0.0f;
  om->base_w       = 0;
  om->ddbound_w    = 0;
  om->ncj_roundoff = 0.0f;	

  for (x = 0; x < p7_NOFFSETS; x++) om->offs[x]    = -1;
  for (x = 0; x < p7_NEVPARAM; x++) om->evparam[x] = p7_EVPARAM_UNSET;
  for (x = 0; x < p7_NCUTOFFS; x++) om->cutoff[x]  = p7_CUTOFF_UNSET;
  for (x = 0; x < p7_MAXABET;  x++) om->compo[x]   = p7_COMPO_UNSET;

  /* in a P7_OPROFILE, we always allocate for the optional RF, CS annotation.  
   * we only rely on the leading \0 to signal that it's unused, but 
   * we initialize all this memory to zeros to shut valgrind up about 
   * fwrite'ing uninitialized memory in the io functions.
   */
  ESL_ALLOC(om->rf,          sizeof(char) * (allocM+2));
  ESL_ALLOC(om->mm,          sizeof(char) * (allocM+2));
  ESL_ALLOC(om->cs,          sizeof(char) * (allocM+2));
  ESL_ALLOC(om->consensus,   sizeof(char) * (allocM+2));
  memset(om->rf,       '\0', sizeof(char) * (allocM+2));
  memset(om->mm,       '\0', sizeof(char) * (allocM+2));
  memset(om->cs,       '\0', sizeof(char) * (allocM+2));
  memset(om->consensus,'\0', sizeof(char) * (allocM+2));

  om->abc        = abc;
  om->L          = 0;
  om->M          = 0;
  om->max_length = -1;
  om->allocM     = allocM;
  om->mode       = p7_NO_MODE;
  om->nj         = 0.0f;
  return om;
ERROR:
  p7_oprofile_Destroy_sse(om);
  return NULL;
#endif /* HAVE_SSE2 */
#ifndef HAVE_SSE2
  return NULL; // Stub so we have something to link when we don't have SSE support
#endif
 
}

/* Function:  p7_oprofile_Destroy()
 * Synopsis:  Frees an optimized profile structure.
 * Incept:    SRE, Sun Nov 25 12:22:21 2007 [Casa de Gatos]
 */
void
p7_oprofile_Destroy_sse(P7_OPROFILE *om)
{
#ifdef HAVE_SSE2
  if (om == NULL) return;

  if (! om->is_shadow)
    {    
      if (om->rbv_mem)   free(om->rbv_mem);
      if (om->sbv_mem)   free(om->sbv_mem);    
      if (om->rwv_mem)   free(om->rwv_mem);
      if (om->twv_mem)   free(om->twv_mem);
     if (om->rfv_mem)   free(om->rfv_mem);
      if (om->tfv_mem)   free(om->tfv_mem);
      if (om->rbv)       free(om->rbv);
      if (om->sbv)       free(om->sbv);
      if (om->rwv)       free(om->rwv);
      if (om->rfv)       free(om->rfv);

      if (om->name)      free(om->name);
      if (om->acc)       free(om->acc);
      if (om->desc)      free(om->desc);
      if (om->rf)        free(om->rf);
      if (om->mm)        free(om->mm);
      if (om->cs)        free(om->cs);
      if (om->consensus) free(om->consensus);
    }
 #endif   
  free(om);  // This will leave dangling allocated memory if we don't have SSE support, but trying to run the SSE
  // code without SSE support will fail spectacularly in so many ways that that's the least of our problems  
}

/* Function:  p7_oprofile_Sizeof()
 * Synopsis:  Return the allocated size of a <P7_OPROFILE>.
 * Incept:    SRE, Wed Mar  2 10:09:21 2011 [Janelia]
 *
 * Purpose:   Returns the allocated size of a <P7_OPROFILE>,
 *            in bytes.
 *            
 *            Very roughly, M*284 bytes, for a model of length M; 60KB
 *            for a typical model; 30MB for a design limit M=100K
 *            model.
 */
size_t
p7_oprofile_Sizeof_sse(const P7_OPROFILE *om)
{
 
  size_t n   = 0;
#ifdef HAVE_SSE2 
  int    nqb = om->allocQ16;  /* # of uchar vectors needed for query */
  int    nqw = om->allocQ8;     /* # of sword vectors needed for query */
  int    nqf = om->allocQ4;     /* # of float vectors needed for query */
  int    nqs = nqb + p7O_EXTRA_SB;

  /* Stuff below exactly mirrors the malloc()'s in
   * p7_oprofile_Create(); so even though we could
   * write this more compactly, leave it like this
   * w/ one:one correspondence to _Create(), for
   * maintainability and clarity.
   */
  n  += sizeof(P7_OPROFILE);
  n  += sizeof(__m128i) * nqb  * om->abc->Kp +15; /* om->rbv_mem   */
  n  += sizeof(__m128i) * nqs  * om->abc->Kp +15; /* om->sbv_mem   */
  n  += sizeof(__m128i) * nqw  * om->abc->Kp +15; /* om->rwv_mem   */
  n  += sizeof(__m128i) * nqw  * p7O_NTRANS  +15; /* om->twv_mem   */
  n  += sizeof(__m128)  * nqf  * om->abc->Kp +15; /* om->rfv_mem   */
  n  += sizeof(__m128)  * nqf  * p7O_NTRANS  +15; /* om->tfv_mem   */
 
  n  += sizeof(__m128i *) * om->abc->Kp;          /* om->rbv       */
  n  += sizeof(__m128i *) * om->abc->Kp;          /* om->sbv       */
  n  += sizeof(__m128i *) * om->abc->Kp;          /* om->rwv       */
 n  += sizeof(__m128  *) * om->abc->Kp;          /* om->rfv       */

  n  += sizeof(char) * (om->allocM+2);            /* om->rf        */
  n  += sizeof(char) * (om->allocM+2);            /* om->mm        */
  n  += sizeof(char) * (om->allocM+2);            /* om->cs        */
  n  += sizeof(char) * (om->allocM+2);            /* om->consensus */
#endif
  return n;
}


/* Function:  p7_oprofile_Clone()
 * Synopsis:  Create a new copy of an optimized profile structure.
 * Incept:    SRE, Sun Nov 25 12:03:19 2007 [Casa de Gatos]
 *
 * Purpose:   Create a newly allocated copy of <om1> and return a ptr
 *            to it.
 *            
 * Throws:    <NULL> on allocation error.
 */
P7_OPROFILE *
p7_oprofile_Clone_sse(const P7_OPROFILE *om1)
{
  const ESL_ALPHABET *abc = om1->abc;
  P7_OPROFILE  *om2  = NULL;
#ifdef HAVE_SSE2
 int           nqb  = P7_NVB(om1->allocM); /* # of uchar vectors needed for query */
  int           nqw  = P7_NVW(om1->allocM); /* # of sword vectors needed for query */
  int           nqf  = P7_NVF(om1->allocM); /* # of float vectors needed for query */
  int           nqs  = nqb + p7O_EXTRA_SB;

  size_t        size = sizeof(char) * (om1->allocM+2);
  int           x, y;
  int           status;

  /* level 0 */
  ESL_ALLOC(om2, sizeof(P7_OPROFILE));

  om2->rbv_mem   = NULL;
  om2->sbv_mem   = NULL;
  om2->rwv_mem   = NULL;
  om2->twv_mem   = NULL;
   om2->rfv_mem   = NULL;
  om2->tfv_mem   = NULL;
 
  om2->rbv       = NULL;
  om2->sbv       = NULL;
  om2->rwv       = NULL;
  om2->twv       = NULL;
    om2->rfv       = NULL;
  om2->tfv       = NULL;

  om2->is_shadow = FALSE;  // om1 can be a shadow, but the resulting copy is a full-fledged profile
  
  om2->name      = NULL;
  om2->acc       = NULL;
  om2->desc      = NULL;
  om2->rf        = NULL;
  om2->mm        = NULL;
  om2->cs        = NULL;
  om2->consensus = NULL;

  /* level 1 */
  ESL_ALLOC(om2->rbv_mem, sizeof(__m128i) * nqb  * abc->Kp    +15); /* +15 is for manual 16-byte alignment */
  ESL_ALLOC(om2->sbv_mem, sizeof(__m128i) * nqs  * abc->Kp    +15);
  ESL_ALLOC(om2->rwv_mem, sizeof(__m128i) * nqw  * abc->Kp    +15);                     
  ESL_ALLOC(om2->twv_mem, sizeof(__m128i) * nqw  * p7O_NTRANS +15); 
  ESL_ALLOC(om2->rfv_mem, sizeof(__m128)  * nqf  * abc->Kp    +15);                     
  ESL_ALLOC(om2->tfv_mem, sizeof(__m128)  * nqf  * p7O_NTRANS +15);      

  ESL_ALLOC(om2->rbv, sizeof(__m128i *) * abc->Kp); 
  ESL_ALLOC(om2->sbv, sizeof(__m128i *) * abc->Kp); 
  ESL_ALLOC(om2->rwv, sizeof(__m128i *) * abc->Kp); 
  ESL_ALLOC(om2->rfv, sizeof(__m128  *) * abc->Kp); 


  /* align vector memory on vector size boundaries */
  om2->rbv[0] = (__m128i *) (((unsigned long int) om2->rbv_mem + 15) & (~0xf));
  om2->sbv[0] = (__m128i *) (((unsigned long int) om2->sbv_mem + 15) & (~0xf));
  memcpy(om2->rbv[0], om1->rbv[0], sizeof(__m128i) * nqb  * abc->Kp);
  memcpy(om2->sbv[0], om1->sbv[0], sizeof(__m128i) * nqs  * abc->Kp);
  om2->rwv[0] = (__m128i *) (((unsigned long int) om2->rwv_mem + 15) & (~0xf));
  om2->twv    = (__m128i *) (((unsigned long int) om2->twv_mem + 15) & (~0xf));
   memcpy(om2->rwv[0], om1->rwv[0], sizeof(__m128i) * nqw  * abc->Kp);
     om2->rfv[0] = (__m128  *) (((unsigned long int) om2->rfv_mem + 15) & (~0xf));
  om2->tfv    = (__m128  *) (((unsigned long int) om2->tfv_mem + 15) & (~0xf));

  /* copy the vector data */

 
  memcpy(om2->rfv[0], om1->rfv[0], sizeof(__m128i) * nqf  * abc->Kp);

  /* set the rest of the row pointers for match emissions */
  for (x = 1; x < abc->Kp; x++) {
    om2->rbv[x] = om2->rbv[0] + (x * nqb);
    om2->sbv[x] = om2->sbv[0] + (x * nqs);
    om2->rwv[x] = om2->rwv[0] + (x * nqw);
    om2->rfv[x] = om2->rfv[0] + (x * nqf);
 
  }
  om2->allocQ16  = nqb;
  om2->allocQ8   = nqw;
  om2->allocQ4   = nqf;

  /* Remaining initializations */
  om2->tbm_b     = om1->tbm_b;
  om2->tec_b     = om1->tec_b;
  om2->tjb_b     = om1->tjb_b;
  om2->scale_b   = om1->scale_b;
  om2->base_b    = om1->base_b;
  om2->bias_b    = om1->bias_b;

  om2->scale_w      = om1->scale_w;
  om2->base_w       = om1->base_w;
  om2->ddbound_w    = om1->ddbound_w;
  om2->ncj_roundoff = om1->ncj_roundoff;	

  for (x = 0; x < p7_NOFFSETS; x++) om2->offs[x]    = om1->offs[x];
  for (x = 0; x < p7_NEVPARAM; x++) om2->evparam[x] = om1->evparam[x];
  for (x = 0; x < p7_NCUTOFFS; x++) om2->cutoff[x]  = om1->cutoff[x];
  for (x = 0; x < p7_MAXABET;  x++) om2->compo[x]   = om1->compo[x];
  for (x = 0; x < nqw  * p7O_NTRANS; ++x) om2->twv[x] = om1->twv[x];
  for (x = 0; x < nqf  * p7O_NTRANS; ++x) om2->tfv[x] = om1->tfv[x];

  for (x = 0; x < p7O_NXSTATES; x++)
    for (y = 0; y < p7O_NXTRANS; y++)
      {
	om2->xw[x][y] = om1->xw[x][y];
	om2->xf[x][y] = om1->xf[x][y];
      }

  if ((status = esl_strdup(om1->name, -1, &om2->name)) != eslOK) goto ERROR;
  if ((status = esl_strdup(om1->acc,  -1, &om2->acc))  != eslOK) goto ERROR;
  if ((status = esl_strdup(om1->desc, -1, &om2->desc)) != eslOK) goto ERROR;

  /* in a P7_OPROFILE, we always allocate for the optional RF, CS annotation.  
   * we only rely on the leading \0 to signal that it's unused, but 
   * we initialize all this memory to zeros to shut valgrind up about 
   * fwrite'ing uninitialized memory in the io functions.
   */
  ESL_ALLOC(om2->rf,          size);
  ESL_ALLOC(om2->mm,          size);
  ESL_ALLOC(om2->cs,          size);
  ESL_ALLOC(om2->consensus,   size);

  memcpy(om2->rf,        om1->rf,        size);
  memcpy(om2->mm,        om1->mm,        size);
  memcpy(om2->cs,        om1->cs,        size);
  memcpy(om2->consensus, om1->consensus, size);

  om2->abc        = om1->abc;
  om2->L          = om1->L;
  om2->M          = om1->M;
  om2->allocM     = om1->allocM;
  om2->mode       = om1->mode;
  om2->nj         = om1->nj;
  om2->max_length = om1->max_length;
  return om2;
 #endif // HAVE_SSE2
 #ifndef HAVE_SSE2
 return NULL;
 #endif 

 ERROR:
  p7_oprofile_Destroy(om2);
  return NULL;
}

/*----------------- end, P7_OPROFILE structure ------------------*/



/*****************************************************************
 * 2. Conversion from generic P7_PROFILE to optimized P7_OPROFILE
 *****************************************************************/

/* sf_conversion():
 * Author: Bjarne Knudsen
 * 
 * Generates the SSVFilter() parts of the profile <om> scores
 * from the completed MSV score.  This includes calculating 
 * special versions of the match scores for using the the
 * ssv filter.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
sf_conversion_sse(P7_OPROFILE *om)
{
   int     x;     /* counter over residues                                        */
  int     q;      /* q counts over total # of striped vectors, 0..nq-1           */
  int     M   = om->M;		/* length of the query                                          */


#ifdef HAVE_SSE2
  __m128i tmp;
  __m128i tmp2;
  int     nq  = P7_NVB(M);     /* segment length; total # of striped vectors needed            */

  /* We now want to fill out om->sbv with om->rbv - bias for use in the
   * SSV filter. The only challenge is that the om->rbv values are
   * unsigned and generally use the whole scale while the om->sbv
   * values are signed. To solve that problem we perform the following
   * calculation:
   *
   *   ((127 + bias) - rbv) ^ 127
   *
   * where the subtraction is unsigned saturated and the addition is
   * unsigned (it will not overflow, since bias is a small positive
   * number). The f(x) = x ^ 127 combined with a change from unsigned
   * to signed numbers have the same effect as f(x) = -x + 127. So if
   * we regard the above as signed instead of unsigned it is equal to:
   *
   *   -((127 + bias) - rbv) + 127 = rbv - bias
   *
   * which is what we want. The reason for this slightly complex idea
   * is that we wish the transformation to be fast, especially for
   * hmmscan where many models are loaded.
   */

  tmp = _mm_set1_epi8((int8_t) (om->bias_b + 127));
  tmp2  = _mm_set1_epi8(127);

  for (x = 0; x < om->abc->Kp; x++)
    {
      for (q = 0;  q < nq;            q++) om->sbv[x][q] = _mm_xor_si128(_mm_subs_epu8(tmp, om->rbv[x][q]), tmp2);
      for (q = nq; q < nq + p7O_EXTRA_SB; q++) om->sbv[x][q] = om->sbv[x][q % nq];
    }

  return eslOK;
#endif // HAVE_SSE2
#ifndef HAVE_SSE2
  return eslENORESULT;
#endif  
}

/* mf_conversion(): 
 * 
 * This builds the MSVFilter() parts of the profile <om>, scores
 * in lspace uchars (16-way parallel), by rescaling, rounding, and
 * casting the scores in <gm>.
 * 
 * Returns <eslOK> on success;
 * throws <eslEINVAL> if <om> hasn't been allocated properly.
 */
int
mf_conversion_sse(const P7_PROFILE *gm, P7_OPROFILE *om)
{
#ifdef HAVE_SSE2  
  int     M   = gm->M;		/* length of the query                                          */
  int     nq  = P7_NVB(M);     /* segment length; total # of striped vectors needed            */
  union { __m128i v; uint8_t i[16]; } tmp; /* used to align and load simd minivectors           */

  float   max = 0.0;		/* maximum residue score: used for unsigned emission score bias */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     z;			/* counter within elements of one SIMD minivector               */

  if (nq > om->allocQ16) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  /* First we determine the basis for the limited-precision MSVFilter scoring system. 
   * Default: 1/3 bit units, base offset 190:  range 0..255 => -190..65 => -63.3..21.7 bits
   * See J2/66, J4/138 for analysis.
   */
  for (x = 0; x < gm->abc->K; x++)  max = ESL_MAX(max, esl_vec_FMax(gm->rsc[x], (M+1)*2));
  om->scale_b = 3.0 / eslCONST_LOG2;                    /* scores in units of third-bits */
  om->base_b  = 190;
  om->bias_b  = unbiased_byteify(om, -1.0 * max);

  /* striped match costs: start at k=1.  */
  for (x = 0; x < gm->abc->Kp; x++)
    for (q = 0, k = 1; q < nq; q++, k++)
      {
  for (z = 0; z < 16; z++) tmp.i[z] = ((k+ z*nq <= M) ? biased_byteify(om, P7P_MSC(gm, k+z*nq, x)) : 255);
  om->rbv[x][q]   = tmp.v;  
      }

  /* transition costs */
  om->tbm_b = unbiased_byteify(om, logf(2.0f / ((float) gm->M * (float) (gm->M+1)))); /* constant B->Mk penalty        */
  om->tec_b = unbiased_byteify(om, logf(0.5f));                                       /* constant multihit E->C = E->J */
  om->tjb_b = unbiased_byteify(om, logf(3.0f / (float) (gm->L+3))); /* this adopts the L setting of the parent profile */

  sf_conversion_sse(om);

  return eslOK;
#endif //HAVE_SSE2
#ifndef HAVE_SSE2
  return eslENORESULT;
#endif  
}


/* vf_conversion(): 
 * 
 * This builds the ViterbiFilter() parts of the profile <om>, scores
 * in lspace swords (8-way parallel), by rescaling, rounding, and
 * casting the scores in <gm>.
 * 
 * Returns <eslOK> on success;
 * throws <eslEINVAL> if <om> hasn't been allocated properly.
 */
int
vf_conversion_sse(const P7_PROFILE *gm, P7_OPROFILE *om)
{
#ifdef HAVE_SSE2  
  int     M   = gm->M;		/* length of the query                                          */
  
  int     nq  = P7_NVW(M);     /* segment length; total # of striped vectors needed            */
  union { __m128i v; int16_t i[8]; } tmp; /* used to align and load simd minivectors            */
  if (nq > om->allocQ8) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     kb;			/* possibly offset base k for loading om's TSC vectors          */
  int     z;			/* counter within elements of one SIMD minivector               */
  int     t;			/* counter over transitions 0..7 = p7O_{BM,MM,IM,DM,MD,MI,II,DD}*/
  int     tg;			/* transition index in gm                                       */
  int     j;			/* counter in interleaved vector arrays in the profile          */
  int     ddtmp;		/* used in finding worst DD transition bound                    */
  int16_t  maxval;		/* used to prevent zero cost II                                 */
  int16_t  val;


  /* First set the basis for the limited-precision scoring system. 
   * Default: 1/500 bit units, base offset 12000:  range -32768..32767 => -44768..20767 => -89.54..41.53 bits
   * See J4/138 for analysis.
   */
  om->scale_w = 500.0 / eslCONST_LOG2;
  om->base_w  = 12000;

  /* striped match scores */
  for (x = 0; x < gm->abc->Kp; x++) {
    for (k = 1, q = 0; q < nq; q++, k++)
      {
	       for (z = 0; z < 8; z++) tmp.i[z] = ((k+ z*nq <= M) ? wordify(om, P7P_MSC(gm, k+z*nq, x)) : -32768);
	       om->rwv[x][q]   = tmp.v;
      }
  }

  /* Transition costs, all but the DD's. */ 
  for (j = 0, k = 1, q = 0; q < nq; q++, k++)
    {
      for (t = p7O_BM; t <= p7O_II; t++) /* this loop of 7 transitions depends on the order in p7o_tsc_e */
	{
	  switch (t) {
	  case p7O_BM: tg = p7P_LM;  kb = k-1; maxval =  0; break; /* gm has tLMk stored off by one! start from k=0 not 1   */
	  case p7O_MM: tg = p7P_MM;  kb = k-1; maxval =  0; break; /* MM, DM, IM vectors are rotated by -1, start from k=0  */
	  case p7O_IM: tg = p7P_IM;  kb = k-1; maxval =  0; break;
	  case p7O_DM: tg = p7P_DM;  kb = k-1; maxval =  0; break;
	  case p7O_MD: tg = p7P_MD;  kb = k;   maxval =  0; break; /* the remaining ones are straight up  */
	  case p7O_MI: tg = p7P_MI;  kb = k;   maxval =  0; break; 
	  case p7O_II: tg = p7P_II;  kb = k;   maxval = -1; break; 
	  }

	  for (z = 0; z < 8; z++) {
	    val      = ((kb+ z*nq < M) ? wordify(om, P7P_TSC(gm, kb+ z*nq, tg)) : -32768);
	    tmp.i[z] = (val <= maxval) ? val : maxval; /* do not allow an II transition cost of 0, or hell may occur. */
	  }
	  om->twv[j++] = tmp.v;
	}
    }

  /* Finally the DD's, which are at the end of the optimized tsc vector; (j is already sitting there) */
  for (k = 1, q = 0; q < nq; q++, k++)
    {
      for (z = 0; z < 8; z++) tmp.i[z] = ((k+ z*nq < M) ? wordify(om, P7P_TSC(gm, k+ z*nq, p7P_DD)) : -32768);
      om->twv[j++] = tmp.v;
    }

  /* Specials. (Actually in same order in om and gm, but we copy in general form anyway.)  */
  /* VF CC,NN,JJ transitions hardcoded zero; -3.0 nat approximation used instead; this papers
   * over a length independence problem, where the approximation weirdly outperforms the
   * exact solution, probably indicating that the model's Pascal distribution is problematic,
   * and the "approximation" is in fact closer to the One True Model, the mythic H4 supermodel.
   * [xref J5/36] 
   */
  om->xw[p7O_E][p7O_LOOP] = wordify(om, gm->xsc[p7P_E][p7P_LOOP]);  
  om->xw[p7O_E][p7O_MOVE] = wordify(om, gm->xsc[p7P_E][p7P_MOVE]);
  om->xw[p7O_N][p7O_MOVE] = wordify(om, gm->xsc[p7P_N][p7P_MOVE]);
  om->xw[p7O_N][p7O_LOOP] = 0;                                        /* was wordify(om, gm->xsc[p7P_N][p7P_LOOP]); */
  om->xw[p7O_C][p7O_MOVE] = wordify(om, gm->xsc[p7P_C][p7P_MOVE]);
  om->xw[p7O_C][p7O_LOOP] = 0;                                        /* was wordify(om, gm->xsc[p7P_C][p7P_LOOP]); */
  om->xw[p7O_J][p7O_MOVE] = wordify(om, gm->xsc[p7P_J][p7P_MOVE]);
  om->xw[p7O_J][p7O_LOOP] = 0;                                        /* was wordify(om, gm->xsc[p7P_J][p7P_LOOP]); */

  om->ncj_roundoff = 0.0; /* goes along with NN=CC=JJ=0, -3.0 nat approximation */
                          /* otherwise, would be = om->scale_w * gm->xsc[p7P_N][p7P_LOOP] -  om->xw[p7O_N][p7O_LOOP];   */
			  /* see J4/150 for discussion of VF error suppression, superceded by the -3.0 nat approximation */

  /* Transition score bound for "lazy F" DD path evaluation (xref J2/52) */
  om->ddbound_w = -32768;	
  for (k = 2; k < M-1; k++) 
    {
      ddtmp         = (int) wordify(om, P7P_TSC(gm, k,   p7P_DD));
      ddtmp        += (int) wordify(om, P7P_TSC(gm, k+1, p7P_DM));
      ddtmp        -= (int) wordify(om, P7P_TSC(gm, k+1, p7P_LM));
      om->ddbound_w = ESL_MAX(om->ddbound_w, ddtmp);
    }

  return eslOK;
#endif //HAVE_SSE2
#ifndef HAVE_SSE2
  return eslENORESULT;
#endif
}


/* fb_conversion()
 * This builds the Forward/Backward part of the optimized profile <om>,
 * where we use odds ratios (not log-odds scores).
 */
int
fb_conversion_sse(const P7_PROFILE *gm, P7_OPROFILE *om)
{
 #ifdef HAVE_SSE2 
  int     M   = gm->M;		/* length of the query                                          */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     kb;			/* possibly offset base k for loading om's TSC vectors          */
  int     z;			/* counter within elements of one SIMD minivector               */
  int     t;			/* counter over transitions 0..7 = p7O_{BM,MM,IM,DM,MD,MI,II,DD}*/
  int     tg;			/* transition index in gm                                       */
  int     j;			/* counter in interleaved vector arrays in the profile          */
  union { __m128 v; float x[4]; } tmp; /* used to align and load simd minivectors               */

  int     nq  = P7_NVF(M);     /* segment length; total # of striped vectors needed            */
  if (nq > om->allocQ4) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  /* striped match scores: start at k=1 */
  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 1, q = 0; q < nq; q++, k++)
      {
	for (z = 0; z < 4; z++) tmp.x[z] = (k+ z*nq <= M) ? P7P_MSC(gm, k+z*nq, x) : -eslINFINITY;
	om->rfv[x][q] = esl_sse_expf(tmp.v);
      }

  /* Transition scores, all but the DD's. */
  for (j = 0, k = 1, q = 0; q < nq; q++, k++)
    {
      for (t = p7O_BM; t <= p7O_II; t++) /* this loop of 7 transitions depends on the order in the definition of p7o_tsc_e */
	{
	  switch (t) {
	  case p7O_BM: tg = p7P_LM;  kb = k-1; break; /* gm has tBMk stored off by one! start from k=0 not 1 */
	  case p7O_MM: tg = p7P_MM;  kb = k-1; break; /* MM, DM, IM quads are rotated by -1, start from k=0  */
	  case p7O_IM: tg = p7P_IM;  kb = k-1; break;
	  case p7O_DM: tg = p7P_DM;  kb = k-1; break;
	  case p7O_MD: tg = p7P_MD;  kb = k;   break; /* the remaining ones are straight up  */
	  case p7O_MI: tg = p7P_MI;  kb = k;   break; 
	  case p7O_II: tg = p7P_II;  kb = k;   break; 
	  }

	  for (z = 0; z < 4; z++) tmp.x[z] = (kb+z*nq < M) ? P7P_TSC(gm, kb+z*nq, tg) : -eslINFINITY;
	  om->tfv[j++] = esl_sse_expf(tmp.v);
	}
    }

  /* And finally the DD's, which are at the end of the optimized tfv vector; (j is already there) */
  for (k = 1, q = 0; q < nq; q++, k++)
    {
      for (z = 0; z < 4; z++) tmp.x[z] = (k+z*nq < M) ? P7P_TSC(gm, k+z*nq, p7P_DD) : -eslINFINITY;
      om->tfv[j++] = esl_sse_expf(tmp.v);
    }

  /* Specials. (These are actually in exactly the same order in om and
   *  gm, but we copy in general form anyway.)
   */
  om->xf[p7O_E][p7O_LOOP] = expf(gm->xsc[p7P_E][p7P_LOOP]);  
  om->xf[p7O_E][p7O_MOVE] = expf(gm->xsc[p7P_E][p7P_MOVE]);
  om->xf[p7O_N][p7O_LOOP] = expf(gm->xsc[p7P_N][p7P_LOOP]);
  om->xf[p7O_N][p7O_MOVE] = expf(gm->xsc[p7P_N][p7P_MOVE]);
  om->xf[p7O_C][p7O_LOOP] = expf(gm->xsc[p7P_C][p7P_LOOP]);
  om->xf[p7O_C][p7O_MOVE] = expf(gm->xsc[p7P_C][p7P_MOVE]);
  om->xf[p7O_J][p7O_LOOP] = expf(gm->xsc[p7P_J][p7P_LOOP]);
  om->xf[p7O_J][p7O_MOVE] = expf(gm->xsc[p7P_J][p7P_MOVE]);

  return eslOK;
 #endif // HAVE_SSE2
 #ifndef HAVE_SSE2
  return eslENORESULT;
#endif 
}

/*------------ end, conversions to P7_OPROFILE ------------------*/

/*******************************************************************
 * 3. Conversion from optimized P7_OPROFILE to compact score arrays
 *******************************************************************/

/* Function:  p7_oprofile_GetFwdTransitionArray()
 * Synopsis:  Retrieve full 32-bit float transition probabilities from an
 *            optimized profile into a flat array
 *
 * Purpose:   Extract an array of <type> (e.g. p7O_II) transition probabilities
 *            from the underlying <om> profile. In SIMD implementations,
 *            these are striped and interleaved, making them difficult to
 *            directly access.
 *
 * Args:      <om>   - optimized profile, containing transition information
 *            <type> - transition type (e.g. p7O_II)
 *            <arr>  - preallocated array into which floats will be placed
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_oprofile_GetFwdTransitionArray_sse(const P7_OPROFILE *om, int type, float *arr )
{
#ifdef HAVE_SSE2
  int     nq  = P7_NVF(om->M);     /* # of striped vectors needed            */
  int i, j;
  union { __m128 v; float x[4]; } tmp; /* used to align and read simd minivectors               */


  for (i=0; i<nq; i++) {
    // because DD transitions are held at the end of the tfv array
    tmp.v = om->tfv[ (type==p7O_DD ?  nq*7+i :  type+7*i) ];
    for (j=0; j<4; j++)
      if ( i+1+ j*nq < om->M+1)
        arr[i+1+ j*nq]      = tmp.x[j];
  }

  return eslOK;
#endif
#ifndef HAVE_SSE2
  return eslENORESULT;
#endif
}

/* Function:  p7_oprofile_GetMSVEmissionScoreArray()
 * Synopsis:  Retrieve MSV residue emission scores from an optimized
 *            profile into an array
 *
 * Purpose:   Extract an implicitly 2D array of 8-bit int MSV residue
 *            emission scores from an optimized profile <om>. <arr> must
 *            be allocated by the calling function to be of size
 *            ( om->abc->Kp * ( om->M  + 1 )), and indexing into the array
 *            is done as  [om->abc->Kp * i +  c ] for character c at
 *            position i.
 *
 *            In SIMD implementations, the residue scores are striped
 *            and interleaved, making them somewhat difficult to
 *            directly access. Faster access is desired, for example,
 *            in SSV back-tracking of a high-scoring diagonal
 *
 * Args:      <om>   - optimized profile, containing transition information
 *            <arr>  - preallocated array into which scores will be placed
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_oprofile_GetMSVEmissionScoreArray_sse(const P7_OPROFILE *om, uint8_t *arr )
{
#ifdef HAVE_SSE2  
  int x, q, z, k;
 int      M   = om->M;    /* length of the query                                          */
  int      K   = om->abc->Kp;
  union { __m128i v; uint8_t i[16]; } tmp; /* used to align and read simd minivectors           */
  int      nq  = P7_NVB(M);     /* segment length; total # of striped vectors needed            */

  int cell_cnt = (om->M + 1) * K;
  for (x = 0; x < K ; x++) {
    for (q = 0, k = 1; q < nq; q++, k++) {
      tmp.v = om->rbv[x][q];
      for (z=0; z<16; z++)
        if (  (K * (k+z*nq) + x) < cell_cnt) 
          arr[ K * (k+z*nq) + x ] = tmp.i[z];
    }
  }

  return eslOK;
 #endif //HAVE_SSE2
 #ifndef HAVE_SSE2
  return eslENORESULT;
#endif
}


/* Function:  p7_oprofile_GetFwdEmissionScoreArray()
 * Synopsis:  Retrieve Fwd (float) residue emission scores from an optimized
 *            profile into an array
 *
 * Purpose:   Extract an implicitly 2D array of 32-bit float Fwd residue
 *            emission scores from an optimized profile <om>. <arr> must
 *            be allocated by the calling function to be of size
 *            ( om->abc->Kp * ( om->M  + 1 )), and indexing into the array
 *            is done as  [om->abc->Kp * i +  c ] for character c at
 *            position i.
 *
 *            In SIMD implementations, the residue scores are striped
 *            and interleaved, making them somewhat difficult to
 *            directly access.
 *
 * Args:      <om>   - optimized profile, containing transition information
 *            <arr>  - preallocated array into which scores will be placed
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_oprofile_GetFwdEmissionScoreArray_sse(const P7_OPROFILE *om, float *arr )
{
#ifdef HAVE_SSE2  
  int x, q, z, k;
  union { __m128 v; float f[4]; } tmp; /* used to align and read simd minivectors               */
  int      M   = om->M;    /* length of the query                                          */
  int      K   = om->abc->Kp;
  int      nq  = P7_NVF(M);     /* segment length; total # of striped vectors needed            */
  int cell_cnt = (om->M + 1) * K;

  for (x = 0; x < K; x++) {
      for (q = 0, k = 1; q < nq; q++, k++) {
        tmp.v = esl_sse_logf(om->rfv[x][q]);
        for (z = 0; z < 4; z++)
          if (  (K * (k+z*nq) + x) < cell_cnt)
            arr[ K * (k+z*nq) + x ] = tmp.f[z];
      }
  }
  return eslOK;
#endif //HAVE_SSE2
 #ifndef HAVE_SSE2
  return eslENORESULT;
#endif
}

/* Function:  p7_oprofile_GetFwdEmissionArray()
 * Synopsis:  Retrieve Fwd (float) residue emission values from an optimized
 *            profile into an array
 *
 * Purpose:   Extract an implicitly 2D array of 32-bit float Fwd residue
 *            emission values from an optimized profile <om>, converting
 *            back to emission values based on the background. <arr> must
 *            be allocated by the calling function to be of size
 *            ( om->abc->Kp * ( om->M  + 1 )), and indexing into the array
 *            is done as  [om->abc->Kp * i +  c ] for character c at
 *            position i.
 *
 *            In SIMD implementations, the residue scores are striped
 *            and interleaved, making them somewhat difficult to
 *            directly access.
 *
 * Args:      <om>   - optimized profile, containing transition information
 *            <bg>   - background frequencies
 *            <arr>  - preallocated array into which scores will be placed
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_oprofile_GetFwdEmissionArray_sse(const P7_OPROFILE *om, P7_BG *bg, float *arr )
{
#ifdef HAVE_SSE2
  int x, q, z, k;
  union { __m128 v; float f[4]; } tmp; /* used to align and read simd minivectors               */
  int      M   = om->M;    /* length of the query                                          */
  int      Kp  = om->abc->Kp;
  int      K   = om->abc->K;
  int      nq  = P7_NVF(M);     /* segment length; total # of striped vectors needed            */
  int cell_cnt = (om->M + 1) * Kp;

  for (x = 0; x < K; x++) {
      for (q = 0, k = 1; q < nq; q++, k++) {
        tmp.v = om->rfv[x][q];
        for (z = 0; z < 4; z++)
          if (  (Kp * (k+z*nq) + x) < cell_cnt)
            arr[ Kp * (k+z*nq) + x ] = tmp.f[z] * bg->f[x];
      }
  }

  //degeneracy emissions for each position
  for (x = 0; x <= M; x++)
    esl_abc_FExpectScVec(om->abc, arr+Kp*x, bg->f);
  return eslOK;
#endif //HAVE_SSE2
 #ifndef HAVE_SSE2
  return eslENORESULT;
#endif
}
/*------------ end, conversions from P7_OPROFILE ------------------*/


/*****************************************************************
 * 4. Debugging and development utilities.
 *****************************************************************/


/* oprofile_dump_mf()
 * 
 * Dump the MSVFilter part of a profile <om> to <stdout>.
 */
int
oprofile_dump_mf_sse(FILE *fp, const P7_OPROFILE *om)
{
#ifdef HAVE_SSE2  
  int     M   = om->M;		/* length of the query                                          */

  int     nq  = P7_NVB(M);     /* segment length; total # of striped vectors needed            */
  union { __m128i v; uint8_t i[16]; } tmp; /* used to align and read simd minivectors           */

  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* counter over nodes 1..M                                      */
  int     z;			/* counter within elements of one SIMD minivector               */

/* This will generate gibberish if more than one of p7_build_SSE, p7_build_AVX2, and p7_build_AVX512
  are set */
  /* Header (rearranged column numbers, in the vectors)  */
  fprintf(fp, "     ");
  for (k =1, q = 0; q < nq; q++, k++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < 16; z++) 
  if (k+z*nq <= M) fprintf(fp, "%4d ", k+z*nq);
  else             fprintf(fp, "%4s ", "xx");
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");

  /* Table of residue emissions */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 

      for (q = 0; q < nq; q++)
  {
    fprintf(fp, "[ ");
    _mm_store_si128(&tmp.v, om->rbv[x][q]);
    for (z = 0; z < 16; z++) fprintf(fp, "%4d ", tmp.i[z]);
    fprintf(fp, "]");
  }
      fprintf(fp, "\n");
    }
  fprintf(fp, "\n");

  fprintf(fp, "t_EC,EJ:    %4d\n",  om->tec_b);
  fprintf(fp, "t_NB,JB,CT: %4d\n",  om->tjb_b);
  fprintf(fp, "t_BMk:      %4d\n",  om->tbm_b);
  fprintf(fp, "scale:      %.2f\n", om->scale_b);
  fprintf(fp, "base:       %4d\n",  om->base_b);
  fprintf(fp, "bias:       %4d\n",  om->bias_b);
  fprintf(fp, "Q:          %4d\n",  nq);  
  fprintf(fp, "M:          %4d\n",  M);  
  return eslOK;
  #endif //HAVE_SSE2
 #ifndef HAVE_SSE2
  return eslENORESULT;
#endif
}


/* oprofile_dump_vf()
 * 
 * Dump the ViterbiFilter part of a profile <om> to <stdout>.
 */
int oprofile_dump_vf_sse(FILE *fp, const P7_OPROFILE *om)
{
 #ifdef HAVE_SSE2 
  int     M   = om->M;		/* length of the query                                          */
  int     nq  = P7_NVW(M);     /* segment length; total # of striped vectors needed            */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     kb;			/* possibly offset base k for loading om's TSC vectors          */
  int     z;			/* counter within elements of one SIMD minivector               */
  int     t;			/* counter over transitions 0..7 = p7O_{BM,MM,IM,DM,MD,MI,II,DD}*/
  int     j;			/* counter in interleaved vector arrays in the profile          */

  union { __m128i v; int16_t i[8]; } tmp; /* used to align and read simd minivectors           */

  /* Emission score header (rearranged column numbers, in the vectors)  */
  fprintf(fp, "     ");
  for (k =1, q = 0; q < nq; q++, k++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < 8; z++) 
	if (k+z*nq <= M) fprintf(fp, "%6d ", k+z*nq);
	else             fprintf(fp, "%6s ", "xx");
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");

  /* Table of residue emissions */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 

      /* Match emission scores (insert emissions are assumed zero by design) */
      for (q = 0; q < nq; q++)
	{
	  fprintf(fp, "[ ");
	  _mm_store_si128(&tmp.v, om->rwv[x][q]);
	  for (z = 0; z < 8; z++) fprintf(fp, "%6d ", tmp.i[z]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");
    }
  fprintf(fp, "\n");

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
	  case p7O_MM: kb = (1 + (nq+k-2)) % nq; break; /* MM, DM, IM quads rotated by +1  */
	  case p7O_IM: kb = (1 + (nq+k-2)) % nq; break;  
	  case p7O_DM: kb = (1 + (nq+k-2)) % nq; break;  
	  case p7O_MD: kb = k;                 break; /* the remaining ones are straight up  */
	  case p7O_MI: kb = k;                 break; 
	  case p7O_II: kb = k;                 break; 
	  }
	  fprintf(fp, "[ ");
	  for (z = 0; z < 8; z++) 
	    if (kb+z*nq <= M) fprintf(fp, "%6d ", kb+z*nq);
	    else              fprintf(fp, "%6s ", "xx");
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n     ");	  
      for (q = 0; q < nq; q++)
	{
	  fprintf(fp, "[ ");
	  _mm_store_si128(&tmp.v, om->twv[q*7 + t]);
	  for (z = 0; z < 8; z++) fprintf(fp, "%6d ", tmp.i[z]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");	  
    }

  /* DD transitions */
  fprintf(fp, "\ntDD: ");
  for (k =1, q = 0; q < nq; q++, k++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < 8; z++) 
	if (k+z*nq <= M) fprintf(fp, "%6d ", k+z*nq);
	else             fprintf(fp, "%6s ", "xx");
      fprintf(fp, "]");
    }
  fprintf(fp, "\n     ");	  
  for (j = nq*7, q = 0; q < nq; q++, j++)
    {
      fprintf(fp, "[ ");
      _mm_store_si128(&tmp.v, om->twv[j]);
      for (z = 0; z < 8; z++) fprintf(fp, "%6d ", tmp.i[z]);
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");	  

  fprintf(fp, "E->C: %6d    E->J: %6d\n", om->xw[p7O_E][p7O_MOVE], om->xw[p7O_E][p7O_LOOP]);
  fprintf(fp, "N->B: %6d    N->N: %6d\n", om->xw[p7O_N][p7O_MOVE], om->xw[p7O_N][p7O_LOOP]);
  fprintf(fp, "J->B: %6d    J->J: %6d\n", om->xw[p7O_J][p7O_MOVE], om->xw[p7O_J][p7O_LOOP]);
  fprintf(fp, "C->T: %6d    C->C: %6d\n", om->xw[p7O_C][p7O_MOVE], om->xw[p7O_C][p7O_LOOP]);

  fprintf(fp, "scale: %6.2f\n", om->scale_w);
  fprintf(fp, "base:  %6d\n",   om->base_w);
  fprintf(fp, "bound: %6d\n",   om->ddbound_w);
  fprintf(fp, "Q:     %6d\n",   nq);  
  fprintf(fp, "M:     %6d\n",   M);  


  return eslOK;
  #endif //HAVE_SSE2
 #ifndef HAVE_SSE2
  return eslENORESULT;
#endif
}


/* oprofile_dump_fb()
 * 
 * Dump the Forward/Backward part of a profile <om> to <stdout>.
 * <width>, <precision> control the floating point output:
 *  8,5 is a reasonable choice for prob space,
 *  5,2 is reasonable for log space.
 */
int
oprofile_dump_fb_sse(FILE *fp, const P7_OPROFILE *om, int width, int precision)
{
#ifdef HAVE_SSE2  
  int     M   = om->M;		/* length of the query                                          */
  int     nq  = P7_NVF(M);     /* segment length; total # of striped vectors needed            */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     kb;			/* possibly offset base k for loading om's TSC vectors          */
  int     z;			/* counter within elements of one SIMD minivector               */
  int     t;			/* counter over transitions 0..7 = p7O_{BM,MM,IM,DM,MD,MI,II,DD}*/
  int     j;			/* counter in interleaved vector arrays in the profile          */
  union { __m128 v; float x[4]; } tmp; /* used to align and read simd minivectors               */
  /* Residue emissions */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 
      for (k =1, q = 0; q < nq; q++, k++)
	{
	  fprintf(fp, "[ ");
	  for (z = 0; z < 4; z++) 
	    if (k+z*nq <= M) fprintf(fp, "%*d ", width, k+z*nq);
	    else             fprintf(fp, "%*s ", width, "xx");
	  fprintf(fp, "]");
	}
      fprintf(fp, "\nmat: ");
      for (q = 0; q < nq; q++)
	{
	  fprintf(fp, "[ ");
	  tmp.v = om->rfv[x][q];
	  for (z = 0; z < 4; z++) fprintf(fp, "%*.*f ", width, precision, tmp.x[z]);
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
	  case p7O_MM:/* MM, DM, IM quads rotated by +1  */
	  case p7O_IM:
	  case p7O_DM:
		  kb = (1 + (nq+k-2)) % nq;
		  break;
	  case p7O_BM:/* the remaining ones are straight up  */
	  case p7O_MD:
	  case p7O_MI:
	  case p7O_II:
		  kb = k;
		  break;
	  }
	  fprintf(fp, "[ ");
	  for (z = 0; z < 4; z++) 
	    if (kb+z*nq <= M) fprintf(fp, "%*d ", width, kb+z*nq);
	    else              fprintf(fp, "%*s ", width, "xx");
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n     ");	  
      for (q = 0; q < nq; q++)
	{
	  fprintf(fp, "[ ");
	  tmp.v = om->tfv[q*7 + t];
	  for (z = 0; z < 4; z++) fprintf(fp, "%*.*f ", width, precision, tmp.x[z]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");	  
    }

  /* DD transitions */
  fprintf(fp, "\ntDD: ");
  for (k =1, q = 0; q < nq; q++, k++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < 4; z++) 
	if (k+z*nq <= M) fprintf(fp, "%*d ", width, k+z*nq);
	else             fprintf(fp, "%*s ", width, "xx");
      fprintf(fp, "]");
    }
  fprintf(fp, "\n     ");	  
  for (j = nq*7, q = 0; q < nq; q++, j++)
    {
      fprintf(fp, "[ ");
      tmp.v = om->tfv[j];
      for (z = 0; z < 4; z++) fprintf(fp, "%*.*f ", width, precision, tmp.x[z]);
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");	  
  
  /* Specials */
  fprintf(fp, "E->C: %*.*f    E->J: %*.*f\n", width, precision, om->xf[p7O_E][p7O_MOVE], width, precision, om->xf[p7O_E][p7O_LOOP]);
  fprintf(fp, "N->B: %*.*f    N->N: %*.*f\n", width, precision, om->xf[p7O_N][p7O_MOVE], width, precision, om->xf[p7O_N][p7O_LOOP]);
  fprintf(fp, "J->B: %*.*f    J->J: %*.*f\n", width, precision, om->xf[p7O_J][p7O_MOVE], width, precision, om->xf[p7O_J][p7O_LOOP]);
  fprintf(fp, "C->T: %*.*f    C->C: %*.*f\n", width, precision, om->xf[p7O_C][p7O_MOVE], width, precision, om->xf[p7O_C][p7O_LOOP]);
  fprintf(fp, "Q:     %d\n",   nq);  
  fprintf(fp, "M:     %d\n",   M);  

  return eslOK;
#endif //HAVE_SSE2
 #ifndef HAVE_SSE2
  return eslENORESULT;
#endif
}

/* Function:  p7_oprofile_Compare()
 * Synopsis:  Compare two optimized profiles for equality.
 * Incept:    SRE, Wed Jan 21 13:29:10 2009 [Janelia]
 *
 * Purpose:   Compare the contents of <om1> and <om2>; return 
 *            <eslOK> if they are effectively identical profiles,
 *            or <eslFAIL> if not.
 * 
 *            Floating point comparisons are done to a tolerance
 *            of <tol> using <esl_FCompare()>.
 *            
 *            If a comparison fails, an informative error message is
 *            left in <errmsg> to indicate why.
 *            
 *            Internal allocation sizes are not compared, only the
 *            data.b
 *            
 * Args:      om1    - one optimized profile to compare
 *            om2    - the other
 *            tol    - floating point comparison tolerance; see <esl_FCompare()>
 *            errmsg - ptr to array of at least <eslERRBUFSIZE> characters.
 *            
 * Returns:   <eslOK> on effective equality;  <eslFAIL> on difference.
 */
int
p7_oprofile_Compare_sse(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg)
{
  if(om1->simd != om2->simd){
    ESL_FAIL(eslFAIL, errmsg, "profiles were built for different SIMD ISAs");
  }

#ifdef HAVE_SSE2 
  int q, r, x, y;
 int Q4  = P7_NVF(om1->M);
  int Q8  = P7_NVW(om1->M);
  int Q16 = P7_NVB(om1->M); 
  union { __m128i v; uint8_t c[16]; } a16, b16;
  union { __m128i v; int16_t w[8];  } a8,  b8;
  union { __m128  v; float   x[4];  } a4,  b4;

  if (om1->mode      != om2->mode)      ESL_FAIL(eslFAIL, errmsg, "comparison failed: mode");
  if (om1->L         != om2->L)         ESL_FAIL(eslFAIL, errmsg, "comparison failed: L");
  if (om1->M         != om2->M)         ESL_FAIL(eslFAIL, errmsg, "comparison failed: M");
  if (om1->nj        != om2->nj)        ESL_FAIL(eslFAIL, errmsg, "comparison failed: nj");
  if (om1->abc->type != om2->abc->type) ESL_FAIL(eslFAIL, errmsg, "comparison failed: alphabet type");

  /* MSVFilter part */
  for (x = 0; x < om1->abc->Kp; x++)
    for (q = 0; q < Q16; q++)
      {
	a16.v = om1->rbv[x][q]; b16.v = om2->rbv[x][q];
	for (r = 0; r < 16; r++) if (a16.c[r] != b16.c[r]) ESL_FAIL(eslFAIL, errmsg, "comparison failed: rb[%d] elem %d", q, r);
      }

  if (om1->tbm_b     != om2->tbm_b)     ESL_FAIL(eslFAIL, errmsg, "comparison failed: tbm_b");
  if (om1->tec_b     != om2->tec_b)     ESL_FAIL(eslFAIL, errmsg, "comparison failed: tec_b");
  if (om1->tjb_b     != om2->tjb_b)     ESL_FAIL(eslFAIL, errmsg, "comparison failed: tjb_b");
  if (om1->scale_b   != om2->scale_b)   ESL_FAIL(eslFAIL, errmsg, "comparison failed: scale_b");
  if (om1->base_b    != om2->base_b)    ESL_FAIL(eslFAIL, errmsg, "comparison failed: base_b");
  if (om1->bias_b    != om2->bias_b)    ESL_FAIL(eslFAIL, errmsg, "comparison failed: bias_b");

  /* ViterbiFilter() part */
  for (x = 0; x < om1->abc->Kp; x++)
    for (q = 0; q < Q8; q++)
      {
	a8.v = om1->rwv[x][q]; b8.v = om2->rwv[x][q];
	for (r = 0; r < 8; r++) if (a8.w[r] != b8.w[r]) ESL_FAIL(eslFAIL, errmsg, "comparison failed: rw[%d] elem %d", q, r);
      }
  for (q = 0; q < 8*Q16; q++)
    {
      a8.v = om1->twv[q]; b8.v = om2->twv[q];
      for (r = 0; r < 8; r++) if (a8.w[r] != b8.w[r]) ESL_FAIL(eslFAIL, errmsg, "comparison failed: tw[%d] elem %d", q, r);
    }

  for (x = 0; x < p7O_NXSTATES; x++)
    for (y = 0; y < p7O_NXTRANS; y++)
      if (om1->xw[x][y] != om2->xw[x][y]) ESL_FAIL(eslFAIL, errmsg, "comparison failed: xw[%d][%d]", x, y);
  if (om1->scale_w   != om2->scale_w)   ESL_FAIL(eslFAIL, errmsg, "comparison failed: scale");
  if (om1->base_w    != om2->base_w)    ESL_FAIL(eslFAIL, errmsg, "comparison failed: base");
  if (om1->ddbound_w != om2->ddbound_w) ESL_FAIL(eslFAIL, errmsg, "comparison failed: ddbound_w");
  
  /* Forward/Backward part */
  for (x = 0; x < om1->abc->Kp; x++)
    for (q = 0; q < Q4; q++)
      {
	a4.v = om1->rfv[x][q]; b4.v = om2->rfv[x][q];
	for (r = 0; r < 4; r++) if (esl_FCompare(a4.x[r], b4.x[r], tol) != eslOK)  ESL_FAIL(eslFAIL, errmsg, "comparison failed: rf[%d] elem %d", q, r);
      }
  for (q = 0; q < 8*Q4; q++)
    {
      a4.v = om1->tfv[q]; b4.v = om2->tfv[q];
      for (r = 0; r < 4; r++) if (a4.x[r] != b4.x[r]) ESL_FAIL(eslFAIL, errmsg, "comparison failed: tf[%d] elem %d", q, r);
    }

  for (x = 0; x < p7O_NXSTATES; x++)
    if (esl_vec_FCompare(om1->xf[x], om2->xf[x], p7O_NXTRANS, tol) != eslOK) ESL_FAIL(eslFAIL, errmsg, "comparison failed: xf[%d] vector", x);

   for (x = 0; x < p7_NOFFSETS; x++)
     if (om1->offs[x] != om2->offs[x]) ESL_FAIL(eslFAIL, errmsg, "comparison failed: offs[%d]", x);

   if (esl_strcmp(om1->name,      om2->name)      != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: name");
   if (esl_strcmp(om1->acc,       om2->acc)       != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: acc");
   if (esl_strcmp(om1->desc,      om2->desc)      != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: desc");
   if (esl_strcmp(om1->rf,        om2->rf)        != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: ref");
   if (esl_strcmp(om1->mm,        om2->mm)        != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: mm");
   if (esl_strcmp(om1->cs,        om2->cs)        != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: cs");
   if (esl_strcmp(om1->consensus, om2->consensus) != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: consensus");
   
   if (esl_vec_FCompare(om1->evparam, om2->evparam, p7_NEVPARAM, tol) != eslOK) ESL_FAIL(eslFAIL, errmsg, "comparison failed: evparam vector");
   if (esl_vec_FCompare(om1->cutoff,  om2->cutoff,  p7_NCUTOFFS, tol) != eslOK) ESL_FAIL(eslFAIL, errmsg, "comparison failed: cutoff vector");
   if (esl_vec_FCompare(om1->compo,   om2->compo,   p7_MAXABET,  tol) != eslOK) ESL_FAIL(eslFAIL, errmsg, "comparison failed: compo vector");

   return eslOK;
   #endif //HAVE_SSE2
 #ifndef HAVE_SSE2
  return eslENORESULT;
#endif
}

/*------------ end, P7_OPROFILE debugging tools  ----------------*/

/*****************************************************************
 * @LICENSE@
 *   
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
