/* Routines for the P7_OPROFILE structure:  
 * a search profile in an optimized implementation.
 * 
 * Contents:
 *   1. The P7_OPROFILE object: allocation, initialization, destruction.
 *   2. Conversion from generic P7_PROFILE to optimized P7_OPROFILE
 *   3. Debugging and development utilities.
 *   4. Benchmark driver.
 *   5. Unit tests.
 *   6. Test driver.
 *   7. Example.
 *   8. Copyright and license information.
 *   
 * SRE, Wed Jul 30 11:00:04 2008 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <string.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_random.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_sse.h"



/*****************************************************************
 * 1. The P7_OPROFILE structure: a score profile.
 *****************************************************************/

/* Function:  p7_oprofile_Create()
 * Synopsis:  Allocate an optimized profile structure.
 * Incept:    SRE, Sun Nov 25 12:03:19 2007 [Casa de Gatos]
 *
 * Purpose:   Allocate for profiles of up to <allocM> nodes for digital alphabet <abc>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_OPROFILE *
p7_oprofile_Create(int allocM, const ESL_ALPHABET *abc)
{
  int          status;
  P7_OPROFILE *om  = NULL;
  int          nqu = p7O_NQU(allocM); /* # of uchar vectors needed for query */
  int          nqf = p7O_NQF(allocM); /* # of float vectors needed for query */
  int          x;

  /* level 0 */
  ESL_ALLOC(om, sizeof(P7_OPROFILE));
  om->tu_mem = NULL;
  om->tf_mem = NULL;
  om->ru_mem = NULL;
  om->rm_mem = NULL;
  om->rf_mem = NULL;
  om->ru     = NULL;
  om->tf     = NULL;
  om->rf     = NULL;

  /* level 1 */
  ESL_ALLOC(om->tu_mem, sizeof(__m128i) * nqu  * p7O_NTRANS       +15);   
  ESL_ALLOC(om->tf_mem, sizeof(__m128)  * nqf  * p7O_NTRANS       +15);    
  ESL_ALLOC(om->ru_mem, sizeof(__m128i) * nqu  * p7O_NR * abc->Kp +15);                     
  ESL_ALLOC(om->rm_mem, sizeof(__m128i) * nqu           * abc->Kp +15);                     
  ESL_ALLOC(om->rf_mem, sizeof(__m128i) * nqf  * p7O_NR * abc->Kp +15);                     
  ESL_ALLOC(om->ru, sizeof(__m128i *) * abc->Kp); 
  ESL_ALLOC(om->rm, sizeof(__m128i *) * abc->Kp); 
  ESL_ALLOC(om->rf, sizeof(__m128 *)  * abc->Kp); 

  /* memory alignment */
  om->tu    = (__m128i *) (((unsigned long int) om->tu_mem + 15) & (~0xf));
  om->tf    = (__m128  *) (((unsigned long int) om->tf_mem + 15) & (~0xf));
  om->ru[0] = (__m128i *) (((unsigned long int) om->ru_mem + 15) & (~0xf));
  om->rm[0] = (__m128i *) (((unsigned long int) om->rm_mem + 15) & (~0xf));
  om->rf[0] = (__m128  *) (((unsigned long int) om->rf_mem + 15) & (~0xf));

  /* the rest of the row pointers */
  for (x = 1; x < abc->Kp; x++) {
    om->ru[x] = om->ru[0] + (x * nqu * p7O_NR);
    om->rm[x] = om->rm[0] + (x * nqu);
    om->rf[x] = om->rf[0] + (x * nqf * p7O_NR);
  }
  om->allocQ16  = nqu;
  om->allocQ4   = nqf;

  /* Remaining initializations */
  om->tbm       = 0;
  om->tec       = 0;
  om->tjb       = 0;

  om->ddbound_u = 0;
  om->scale     = 0.0f;
  om->base      = 0;
  om->bias      = 0;

  om->ddbound_f = 0.0f;
  om->lspace_f  = FALSE;

  for (x = 0; x < p7_NOFFSETS; x++) om->offs[x]    = -1;
  for (x = 0; x < p7_NEVPARAM; x++) om->evparam[x] = 0.0f;
  for (x = 0; x < p7_NCUTOFFS; x++) om->cutoff[x]  = 0.0f;

  om->name      = NULL;
  om->acc       = NULL;
  om->desc      = NULL;

  ESL_ALLOC(om->ref,       sizeof(char) * (allocM+2)); 
  ESL_ALLOC(om->cs,        sizeof(char) * (allocM+2));
  ESL_ALLOC(om->consensus, sizeof(char) * (allocM+2));
  om->ref[0]       = '\0';
  om->cs[0]        = '\0';
  om->consensus[0] = '\0';

  om->mode      = p7_NO_MODE;
  om->L         = 0;
  om->allocM    = allocM;
  om->M         = 0;
  om->nj        = 0.0f;
  om->abc       = abc;
  return om;

 ERROR:
  p7_oprofile_Destroy(om);
  return NULL;
}

/* Function:  p7_oprofile_IsLocal()
 * Synopsis:  Returns TRUE if profile is in local alignment mode.
 * Incept:    SRE, Sat Aug 16 08:46:00 2008 [Janelia]
 */
int
p7_oprofile_IsLocal(const P7_OPROFILE *om)
{
  if (om->mode == p7_LOCAL || om->mode == p7_UNILOCAL) return TRUE;
  return FALSE;
}



/* Function:  p7_oprofile_Destroy()
 * Synopsis:  Frees an optimized profile structure.
 * Incept:    SRE, Sun Nov 25 12:22:21 2007 [Casa de Gatos]
 */
void
p7_oprofile_Destroy(P7_OPROFILE *om)
{
  if (om == NULL) return;

  if (om->tu_mem    != NULL) free(om->tu_mem);
  if (om->ru_mem    != NULL) free(om->ru_mem);
  if (om->rm_mem    != NULL) free(om->rm_mem);
  if (om->tf_mem    != NULL) free(om->tf_mem);
  if (om->rf_mem    != NULL) free(om->rf_mem);
  if (om->ru        != NULL) free(om->ru);
  if (om->rm        != NULL) free(om->rm);
  if (om->rf        != NULL) free(om->rf);
  if (om->name      != NULL) free(om->name);
  if (om->acc       != NULL) free(om->acc);
  if (om->desc      != NULL) free(om->desc);
  if (om->ref       != NULL) free(om->ref);
  if (om->cs        != NULL) free(om->cs);
  if (om->consensus != NULL) free(om->consensus);
  free(om);
}
/*----------------- end, P7_OPROFILE structure ------------------*/



/*****************************************************************
 * 2. Conversion from generic P7_PROFILE to optimized P7_OPROFILE
 *****************************************************************/
static uint8_t biased_charify  (P7_OPROFILE *om, float sc);
static uint8_t unbiased_charify(P7_OPROFILE *om, float sc);
static int     lspace_uchar_conversion(const P7_PROFILE *gm, P7_OPROFILE *om);
static int     lspace_float_conversion(const P7_PROFILE *gm, P7_OPROFILE *om);
static int     pspace_float_conversion(const P7_PROFILE *gm, P7_OPROFILE *om);


/* Function:  p7_oprofile_Convert()
 * Synopsis:  Converts standard profile to an optimized one.
 * Incept:    SRE, Mon Nov 26 07:38:57 2007 [Janelia]
 *
 * Purpose:   Convert a standard profile <gm> to an optimized profile <om>,
 *            where <om> has already been allocated for a profile of at 
 *            least <gm->M> nodes and the same emission alphabet <gm->abc>.
 *
 * Args:      gm - profile to optimize
 *            om - allocated optimized profile for holding the result.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <gm>, <om> aren't compatible. 
 *            <eslEMEM> on allocation failure.
 */
int
p7_oprofile_Convert(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int status, z;

  if (gm->abc->type != om->abc->type)  ESL_EXCEPTION(eslEINVAL, "alphabets of the two profiles don't match");  
  if (gm->M         >  om->allocM)     ESL_EXCEPTION(eslEINVAL, "oprofile is too small");  

  if ((status =  lspace_uchar_conversion(gm, om)) != eslOK) return status;   /* ViterbiFilter()'s information */
  if ((status =  pspace_float_conversion(gm, om)) != eslOK) return status;   /* ForwardFilter()'s information */

  if ((status = esl_strdup(gm->name, -1, &(om->name))) != eslOK) goto ERROR;
  if ((status = esl_strdup(gm->acc,  -1, &(om->acc)))  != eslOK) goto ERROR;
  if ((status = esl_strdup(gm->desc, -1, &(om->desc))) != eslOK) goto ERROR;
  strcpy(om->ref,       gm->rf);
  strcpy(om->cs,        gm->cs);
  strcpy(om->consensus, gm->consensus);
  for (z = 0; z < p7_NEVPARAM; z++) om->evparam[z] = gm->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) om->cutoff[z]  = gm->cutoff[z];

  /* MSVFilter's constants */
  om->tbm  = unbiased_charify(om, logf(2.0f / ((float) gm->M * (float) (gm->M+1)))); /* constant B->Mk penalty        */
  om->tec  = unbiased_charify(om, logf(0.5f));                                       /* constant multihit E->C = E->J */
  /* tjb is length dependent; see ReconfigLength() for setting, below */

  om->mode = gm->mode;
  om->L    = gm->L;
  om->M    = gm->M;
  om->nj   = gm->nj;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_oprofile_ReconfigLength()
 * Synopsis:  Set the target sequence length of a model.
 * Incept:    SRE, Thu Dec 20 09:56:40 2007 [Janelia]
 *
 * Purpose:   Given an already configured model <om>, quickly reset its
 *            expected length distribution for a new mean target sequence
 *            length of <L>. 
 *            
 *            This doesn't affect the length distribution of the null
 *            model. That must also be reset, using <p7_bg_SetLength()>.
 *            
 *            We want this routine to run as fast as possible, because
 *            this call is in the critical path: it must be called at
 *            each new target sequence in a database search.
 *
 * Returns:   <eslOK> on success. Costs/scores for N,C,J transitions are set
 *            here.
 */
int
p7_oprofile_ReconfigLength(P7_OPROFILE *om, int L)
{
  float pmove, ploop;
  
  pmove = (2.0f + om->nj) / ((float) L + 2.0f + om->nj); /* 2/(L+2) for sw; 3/(L+3) for fs */
  ploop = 1.0f - pmove;

  if (om->lspace_f) {  /* ViterbiScore(): lspace floats */
    om->xf[p7O_N][p7O_LOOP] =  om->xf[p7O_C][p7O_LOOP] = om->xf[p7O_J][p7O_LOOP] = logf(ploop);
    om->xf[p7O_N][p7O_MOVE] =  om->xf[p7O_C][p7O_MOVE] = om->xf[p7O_J][p7O_MOVE] = logf(pmove);
  } else {   /* ForwardFilter() parameters: pspace floats */
    om->xf[p7O_N][p7O_LOOP] =  om->xf[p7O_C][p7O_LOOP] = om->xf[p7O_J][p7O_LOOP] = ploop;
    om->xf[p7O_N][p7O_MOVE] =  om->xf[p7O_C][p7O_MOVE] = om->xf[p7O_J][p7O_MOVE] = pmove;
  }
  
  /* MSVFilter() */
  om->tjb = unbiased_charify(om, logf(3.0f / (float) (L+3)));

  /* ViterbiFilter() parameters: lspace uchars; the LOOP costs are zero  */
  om->xu[p7O_N][p7O_MOVE] =  om->xu[p7O_C][p7O_MOVE] = om->xu[p7O_J][p7O_MOVE] = unbiased_charify(om, logf(pmove));
  om->L = L;
  return eslOK;
}

/* Function:  p7_oprofile_ReconfigMultihit()
 * Synopsis:  Quickly reconfig model into multihit mode for target length <L>.
 * Incept:    SRE, Thu Aug 21 10:04:07 2008 [Janelia]
 *
 * Purpose:   Given a profile <om> that's already been configured once,
 *            quickly reconfigure it into a multihit mode for target 
 *            length <L>. 
 *            
 *            This gets called in domain definition, when we need to
 *            flip the model in and out of unihit mode to
 *            process individual domains.
 *            
 * Note:      You can't just flip uni/multi mode alone, because that
 *            parameterization also affects target length
 *            modeling. You need to make sure uni vs. multi choice is
 *            made before the length model is set, and you need to
 *            make sure the length model is recalculated if you change
 *            the uni/multi mode. Hence, these functions call
 *            <p7_oprofile_ReconfigLength()>.
 */
int
p7_oprofile_ReconfigMultihit(P7_OPROFILE *om, int L)
{
  if (om->lspace_f) {
    om->xf[p7O_E][p7O_MOVE] = -eslCONST_LOG2;
    om->xf[p7O_E][p7O_LOOP] = -eslCONST_LOG2;
    om->nj = 1.0f;
  } else {
    om->xf[p7O_E][p7O_MOVE] = 0.5;
    om->xf[p7O_E][p7O_LOOP] = 0.5;
    om->nj = 1.0f;
  }

  om->xu[p7O_E][p7O_MOVE] = unbiased_charify(om, -eslCONST_LOG2);
  om->xu[p7O_E][p7O_LOOP] = unbiased_charify(om, -eslCONST_LOG2);

  return p7_oprofile_ReconfigLength(om, L);
}

/* Function:  p7_oprofile_ReconfigUnihit()
 * Synopsis:  Quickly reconfig model into unihit mode for target length <L>.
 * Incept:    SRE, Thu Aug 21 10:10:32 2008 [Janelia]
 *
 * Purpose:   Given a profile <om> that's already been configured once,
 *            quickly reconfigure it into a unihit mode for target 
 *            length <L>. 
 *            
 *            This gets called in domain definition, when we need to
 *            flip the model in and out of unihit <L=0> mode to
 *            process individual domains.
 */
int
p7_oprofile_ReconfigUnihit(P7_OPROFILE *om, int L)
{
  
  if (om->lspace_f) {
    om->xf[p7O_E][p7O_MOVE] = 0.0f;
    om->xf[p7O_E][p7O_LOOP] = -eslINFINITY;
    om->nj = 0.0f;
  } else {
    om->xf[p7O_E][p7O_MOVE] = 1.0f;
    om->xf[p7O_E][p7O_LOOP] = 0.0f;
    om->nj = 0.0f;
  }

  om->xu[p7O_E][p7O_MOVE] = 255;
  om->xu[p7O_E][p7O_LOOP] = 0;

  return p7_oprofile_ReconfigLength(om, L);
}



/* Function:  p7_oprofile_Logify()
 * Synopsis:  Convert existing model's float scores to lspace.
 * Incept:    SRE, Sun Aug  3 13:24:38 2008 [St. Louis]
 *
 * Purpose:   Convert a model <om> that currently has its float scores
 *            (<om->tf>, <om->rf>, <om->xf>) in probability space (the
 *            default) to one in which these scores are in log 
 *            probability space, suitable for a <p7_ViterbiScore()>
 *            call. This means just taking the log of all the
 *            floating-point scores.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_oprofile_Logify(P7_OPROFILE *om)
{
  int nq  = p7O_NQF(om->M);   /* segment length; total # of striped vectors needed            */
  int x;
  int j;

  if (om->lspace_f == TRUE) return eslOK;

  for (x = 0; x < om->abc->Kp; x++)
    for (j = 0; j < nq*2; j++)
      om->rf[x][j] = esl_sse_logf(om->rf[x][j]);

  for (j = 0; j < nq*p7O_NTRANS; j++)
    om->tf[j] = esl_sse_logf(om->tf[j]);

  om->xf[p7O_E][p7O_LOOP] = logf(om->xf[p7O_E][p7O_LOOP]);  
  om->xf[p7O_E][p7O_MOVE] = logf(om->xf[p7O_E][p7O_MOVE]);
  om->xf[p7O_N][p7O_LOOP] = logf(om->xf[p7O_N][p7O_LOOP]);
  om->xf[p7O_N][p7O_MOVE] = logf(om->xf[p7O_N][p7O_MOVE]);
  om->xf[p7O_C][p7O_LOOP] = logf(om->xf[p7O_C][p7O_LOOP]);
  om->xf[p7O_C][p7O_MOVE] = logf(om->xf[p7O_C][p7O_MOVE]);
  om->xf[p7O_J][p7O_LOOP] = logf(om->xf[p7O_J][p7O_LOOP]);
  om->xf[p7O_J][p7O_MOVE] = logf(om->xf[p7O_J][p7O_MOVE]);
  om->lspace_f = TRUE;
  return eslOK;
}

/* Function:  p7_oprofile_Probify()
 * Synopsis:  Convert existing model's float scores to pspace.
 * Incept:    SRE, Sun Aug  3 13:24:38 2008 [St. Louis]
 *
 * Purpose:   Convert a model <om> that currently has its float scores
 *            (<om->tf>, <om->rf>, <om->xf>) in log probability space 
 *            back to the default probability space, by exponentiating
 *            each individual score.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_oprofile_Probify(P7_OPROFILE *om)
{
  int nq  = p7O_NQF(om->M);   /* segment length; total # of striped vectors needed            */
  int x;
  int j;

  if (om->lspace_f == FALSE) return eslOK;

  for (x = 0; x < om->abc->Kp; x++)
    for (j = 0; j < nq*2; j++)
      om->rf[x][j] = esl_sse_expf(om->rf[x][j]);

  for (j = 0; j < nq*p7O_NTRANS; j++)
    om->tf[j] = esl_sse_expf(om->tf[j]);

  om->xf[p7O_E][p7O_LOOP] = expf(om->xf[p7O_E][p7O_LOOP]);  
  om->xf[p7O_E][p7O_MOVE] = expf(om->xf[p7O_E][p7O_MOVE]);
  om->xf[p7O_N][p7O_LOOP] = expf(om->xf[p7O_N][p7O_LOOP]);
  om->xf[p7O_N][p7O_MOVE] = expf(om->xf[p7O_N][p7O_MOVE]);
  om->xf[p7O_C][p7O_LOOP] = expf(om->xf[p7O_C][p7O_LOOP]);
  om->xf[p7O_C][p7O_MOVE] = expf(om->xf[p7O_C][p7O_MOVE]);
  om->xf[p7O_J][p7O_LOOP] = expf(om->xf[p7O_J][p7O_LOOP]);
  om->xf[p7O_J][p7O_MOVE] = expf(om->xf[p7O_J][p7O_MOVE]);
  om->lspace_f = FALSE;
  return eslOK;
}



/* biased_charify()
 * converts a log-odds residue score to a rounded biased uchar cost.
 * e.g. a score of +3.2, with scale 3.0 and bias 12, becomes 2.
 *    3.2*3 = 9.6; rounded = 10; bias-10 = 2.
 * When used, we add the bias, then subtract this cost.
 *  
 */
static uint8_t
biased_charify(P7_OPROFILE *om, float sc)
{
  uint8_t b;

  sc  = -1.0f * roundf(om->scale * sc);        	        /* ugh. sc is now an integer cost represented in a float...           */
  b   = (sc > 255.) ? 255 : (uint8_t) sc + om->bias;	/* and now we cast, saturate, and bias it to an unsigned char cost... */
  return b;
}
 
/* unbiased_charify()
 * converts log probability score to a rounded uchar cost.
 * e.g. a score of -2.1, with scale 3.0, becomes a cost of 6.
 */
static uint8_t 
unbiased_charify(P7_OPROFILE *om, float sc)
{
  sc  = -1. * roundf(om->scale * sc);          
  return ((sc > 255.) ? 255 : (uint8_t) sc);	
}

/* lspace_uchar_conversion(): 
 * 
 * This builds the ViterbiFilter() parts of the profile <om>, scores
 * in lspace uchars (16-way parallel), by rescaling, rounding, and
 * casting the scores in <gm>.
 * 
 * Returns <eslOK> on success, or exceptions.
 */
static int
lspace_uchar_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int     M   = gm->M;		/* length of the query                                          */
  int     nq  = p7O_NQU(M);     /* segment length; total # of striped vectors needed            */
  float   max = 0.0;		/* maximum residue score: used for unsigned emission score bias */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     kb;			/* possibly offset base k for loading om's TSC vectors          */
  int     z;			/* counter within elements of one SIMD minivector               */
  int     t;			/* counter over transitions 0..7 = p7O_{BM,MM,IM,DM,MD,MI,II,DD}*/
  int     tg;			/* transition index in gm                                       */
  int     j;			/* counter in interleaved vector arrays in the profile          */
  int     ddtmp;		/* used in finding worst DD transition bound                    */
  union { __m128i v; uint8_t i[16]; } tmp; /* used to align and load simd minivectors           */

  if (nq > om->allocQ16) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  /* First we determine the basis for the limited-precision scoring system. 
   * We could implement a more smooth func for determining om->scale below,  
   * and squeaking more precision out; [xref J2/66] for considerations.
   */
  for (x = 0; x < gm->abc->K; x++)
    max = ESL_MAX(max, esl_vec_FMax(gm->rsc[x], (M+1)*2));
  if (gm->M <=800)  om->scale = 3.0 / eslCONST_LOG2;    /* scores in units of third-bits for most models [xref J2/66]     */
  else              om->scale = 2.0 / eslCONST_LOG2;    /* scores in units of half-bits for large models                  */
  om->base  = 195;		                        /* bias in 0..255 scores in DP matrix: scores range -65..20 bits  */
  om->bias  = -1 * unbiased_charify(om, max);           /* match, insert emission costs all negative, by subtracting bias */

  /* striped and interleaved match, insert emission costs: start at k=1.  */
  for (x = 0; x < gm->abc->Kp; x++)
    for (j = 0, k = 1, q = 0; q < nq; q++, k++)
      {
	for (z = 0; z < 16; z++) tmp.i[z] = ((k+ z*nq <= M) ? biased_charify(om, p7P_MSC(gm, k+z*nq, x)) : 255);
	om->ru[x][j++] = tmp.v;

	om->rm[x][q]   = tmp.v;	/* MSVFilter's Mk scores are stored directly, not interleaved. */

	for (z = 0; z < 16; z++) tmp.i[z] = ((k+ z*nq <= M) ? biased_charify(om, p7P_ISC(gm, k+z*nq, x)) : 255);
	om->ru[x][j++] = tmp.v;
      }

  /* Transition costs, all but the DD's. */
  for (j = 0, k = 1, q = 0; q < nq; q++, k++)
    {
      for (t = p7O_BM; t <= p7O_II; t++) /* this loop of 7 transitions depends on the order in p7o_tsc_e */
	{
	  switch (t) {
	  case p7O_BM: tg = p7P_BM;  kb = k-1; break; /* gm has tBMk stored off by one! start from k=0 not 1   */
	  case p7O_MM: tg = p7P_MM;  kb = k-1; break; /* MM, DM, IM vectors are rotated by -1, start from k=0  */
	  case p7O_IM: tg = p7P_IM;  kb = k-1; break;
	  case p7O_DM: tg = p7P_DM;  kb = k-1; break;
	  case p7O_MD: tg = p7P_MD;  kb = k;   break; /* the remaining ones are straight up  */
	  case p7O_MI: tg = p7P_MI;  kb = k;   break; 
	  case p7O_II: tg = p7P_II;  kb = k;   break; 
	  }

	  for (z = 0; z < 16; z++) tmp.i[z] = ((kb+ z*nq < M) ? unbiased_charify(om, p7P_TSC(gm, kb+ z*nq, tg)) : 255);
	  om->tu[j++] = tmp.v;
	}
    }

  /* Finally the DD's, which are at the end of the optimized tsc vector; (j is already sitting there) */
  for (k = 1, q = 0; q < nq; q++, k++)
    {
      for (z = 0; z < 16; z++) tmp.i[z] = ((k+ z*nq < M) ? unbiased_charify(om, p7P_TSC(gm, k+ z*nq, p7P_DD)) : 255);
      om->tu[j++] = tmp.v;
    }

  /* Specials. (Actually in same order in om and gm, but we copy in general form anyway.)
   * NN, CC, JJ are all hardcoded 0; part of precision-maximizing strategy [xref J2/66].
   */
  om->xu[p7O_E][p7O_LOOP] = unbiased_charify(om, gm->xsc[p7P_E][p7P_LOOP]);  
  om->xu[p7O_E][p7O_MOVE] = unbiased_charify(om, gm->xsc[p7P_E][p7P_MOVE]);
  om->xu[p7O_N][p7O_LOOP] = 0;
  om->xu[p7O_N][p7O_MOVE] = unbiased_charify(om, gm->xsc[p7P_N][p7P_MOVE]);
  om->xu[p7O_C][p7O_LOOP] = 0;
  om->xu[p7O_C][p7O_MOVE] = unbiased_charify(om, gm->xsc[p7P_C][p7P_MOVE]);
  om->xu[p7O_J][p7O_LOOP] = 0;
  om->xu[p7O_J][p7O_MOVE] = unbiased_charify(om, gm->xsc[p7P_J][p7P_MOVE]);

  /* Transition score bound for "lazy F" DD path evaluation (xref J2/52) */
  om->ddbound_u = -9999;	
  for (k = 2; k < M-1; k++) 
    {
      ddtmp         = (int) unbiased_charify(om, p7P_TSC(gm, k+1, p7P_BM));
      ddtmp        -= (int) unbiased_charify(om, p7P_TSC(gm, k,   p7P_DD));
      ddtmp        -= (int) unbiased_charify(om, p7P_TSC(gm, k+1, p7P_DM));
      om->ddbound_u = ESL_MAX(om->ddbound_u, ddtmp);
    }

  return eslOK;
}

/* lspace_float_conversion()
 * 
 * This rearranges log-odds and log-prob scores from <gm> into the
 * striped, interleaved optimized profile <om>.
 * 
 * This is the form used by ViterbiScore(), and it's also currently an
 * intermediate for pspace_float_conversion() for ForwardFilter().
 */
static int
lspace_float_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int     M   = gm->M;		/* length of the query                                          */
  int     nq  = p7O_NQF(M);     /* segment length; total # of striped vectors needed            */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     kb;			/* possibly offset base k for loading om's TSC vectors          */
  int     z;			/* counter within elements of one SIMD minivector               */
  int     t;			/* counter over transitions 0..7 = p7O_{BM,MM,IM,DM,MD,MI,II,DD}*/
  int     tg;			/* transition index in gm                                       */
  int     j;			/* counter in interleaved vector arrays in the profile          */
  union { __m128 v; float x[4]; } tmp; /* used to align and load simd minivectors               */

  if (nq > om->allocQ4) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  /* striped and interleaved match, insert emission scores: start at k=1 */
  for (x = 0; x < gm->abc->Kp; x++)
    for (j = 0, k = 1, q = 0; q < nq; q++, k++)
      {
	for (z = 0; z < 4; z++) tmp.x[z] = (k+ z*nq <= M) ? p7P_MSC(gm, k+z*nq, x) : -eslINFINITY;
	om->rf[x][j++] = tmp.v;

	for (z = 0; z < 4; z++) tmp.x[z] = (k+ z*nq <= M) ? p7P_ISC(gm, k+z*nq, x) : -eslINFINITY;
	om->rf[x][j++] = tmp.v;
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

	  for (z = 0; z < 4; z++) tmp.x[z] = (kb+z*nq < M) ? p7P_TSC(gm, kb+z*nq, tg) : -eslINFINITY;
	  om->tf[j++] = tmp.v;
	}
    }

  /* And finally the DD's, which are at the end of the optimized tsc vector; (j is already there) */
  for (k = 1, q = 0; q < nq; q++, k++)
    {
      for (z = 0; z < 4; z++) tmp.x[z] = (k+z*nq < M) ? p7P_TSC(gm, k+z*nq, p7P_DD) : -eslINFINITY;
      om->tf[j++] = tmp.v;
    }

  /* Specials. (These are actually in exactly the same order in om and
   *  gm, but we copy in general form anyway.)
   */
  om->xf[p7O_E][p7O_LOOP] = gm->xsc[p7P_E][p7P_LOOP];  
  om->xf[p7O_E][p7O_MOVE] = gm->xsc[p7P_E][p7P_MOVE];
  om->xf[p7O_N][p7O_LOOP] = gm->xsc[p7P_N][p7P_LOOP];
  om->xf[p7O_N][p7O_MOVE] = gm->xsc[p7P_N][p7P_MOVE];
  om->xf[p7O_C][p7O_LOOP] = gm->xsc[p7P_C][p7P_LOOP];
  om->xf[p7O_C][p7O_MOVE] = gm->xsc[p7P_C][p7P_MOVE];
  om->xf[p7O_J][p7O_LOOP] = gm->xsc[p7P_J][p7P_LOOP];
  om->xf[p7O_J][p7O_MOVE] = gm->xsc[p7P_J][p7P_MOVE];

  /* Transition score bound for "lazy F" DD path evaluation (xref J2/52) */
  om->ddbound_f = -eslINFINITY;	
  for (k = 2; k < M-1; k++)
    om->ddbound_f = ESL_MAX(p7P_TSC(gm, k, p7P_DD) + p7P_TSC(gm, k+1, p7P_DM) - p7P_TSC(gm, k+1, p7P_BM), om->ddbound_f);

  om->lspace_f = TRUE;
  om->mode     = gm->mode;
  om->M        = M;
  return eslOK;
}


/* pspace_float_conversion()
 * 
 * This builds the ForwardFilter() parts of the profile <om> --
 * scores in pspace floats, 4-way parallel - by exponentiating
 * the log-odds and log-prob scores in <gm>.
 * 
 * NOTE: For now, we do this by first loading the log-odds and
 * log-prob scores with lspace_float_conversion(), then exponentiating
 * everything (which we can do in SSE).  When we get to hmmpfam, where
 * profile conversion time will be in critical path, this routine
 * becomes an optimization target.  We can optimize by replacing exp()
 * of <gm> with calculation of odds ratios using the original
 * <hmm>. To do this, we might want <gm> to store the B->Mk
 * distribution as a probability distribution, not just as log probs.
 */
static int
pspace_float_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int     nq  = p7O_NQF(gm->M);  /* segment length; total # of striped vectors needed            */
  int     status;

  if (nq > om->allocQ4) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  /* First stripe and interleave the scores from <gm> */
  if ((status = lspace_float_conversion(gm, om)) != eslOK) return status;

  /* Then exponentiate them all quickly and in stride (using SSE) */
  if ((status = p7_oprofile_Probify(om))         != eslOK) return status;

  om->mode     = gm->mode;
  om->M        = gm->M;
  return eslOK;
}
/*------------ end, conversions to P7_OPROFILE ------------------*/



/*****************************************************************
 * 3. Debugging and development utilities.
 *****************************************************************/
static int oprofile_dump_uchar(FILE *fp, const P7_OPROFILE *om);
static int oprofile_dump_float(FILE *fp, const P7_OPROFILE *om, int width, int precision);

/* Function:  p7_oprofile_Sample()
 * Synopsis:  Sample a random profile.
 * Incept:    SRE, Wed Jul 30 13:11:52 2008 [Janelia]
 *
 * Purpose:   Sample a random profile of <M> nodes for alphabet <abc>,
 *            using <r> as the source of random numbers. Parameterize
 *            it for generation of target sequences of mean length
 *            <L>. Calculate its log-odds scores using background
 *            model <bg>.
 *            
 * Args:      r       - random number generator
 *            abc     - emission alphabet 
 *            bg      - background frequency model
 *            M       - size of sampled profile, in nodes
 *            L       - configured target seq mean length
 *            opt_hmm - optRETURN: sampled HMM
 *            opt_gm  - optRETURN: sampled normal profile
 *            opt_om  - RETURN: optimized profile
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_oprofile_Sample(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, const P7_BG *bg, int M, int L,
		   P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_OPROFILE **ret_om)
{
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_OPROFILE    *om   = NULL;
  int             status;

  if ((gm = p7_profile_Create (M, abc)) == NULL)  { status = eslEMEM; goto ERROR; }
  if ((om = p7_oprofile_Create(M, abc)) == NULL)  { status = eslEMEM; goto ERROR; }

  if ((status = p7_hmm_Sample(r, M, abc, &hmm))             != eslOK) goto ERROR;
  if ((status = p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL)) != eslOK) goto ERROR;
  if ((status = p7_oprofile_Convert(gm, om))                != eslOK) goto ERROR;
  if ((status = p7_oprofile_ReconfigLength(om, L))          != eslOK) goto ERROR;

  if (opt_hmm != NULL) *opt_hmm = hmm; else p7_hmm_Destroy(hmm);
  if (opt_gm  != NULL) *opt_gm  = gm;  else p7_profile_Destroy(gm);
  *ret_om = om;
  return eslOK;

 ERROR:
  if (opt_hmm != NULL) *opt_hmm = NULL;
  if (opt_gm  != NULL) *opt_gm  = NULL;
  *ret_om = NULL;
  return status;
}

/* Function:  p7_oprofile_SameRounding()
 * Synopsis:  Round a generic lspace profile to match lspace uchar oprofile.
 * Incept:    SRE, Wed Jul 30 13:37:48 2008 [Janelia]
 *
 * Purpose:   Round all the scores in a generic (lspace) <P7_PROFILE> in
 *            exactly the same way that the scores in a lspace uchar
 *            <P7_OPROFILE> were rounded.  Then the two profiles
 *            should give identical internal scores in testing, say,
 *            <p7_ViterbiFilter()> against <p7_GViterbi()>.
 *
 *            <gm> must be the same profile that <om> was constructed from.
 * 
 *            <gm> is irrevocably altered by this call. 
 *            
 *            Do not call this more than once on any given <gm>! 
 *
 * Args:      <om>  - optimized profile, containing scale information.
 *            <gm>  - generic profile that <om> was built from.          
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_oprofile_SameRounding(const P7_OPROFILE *om, P7_PROFILE *gm)
{
  int k;
  int x;

  /* Transitions */
  /* <= -eslINFINITY test is used solely to silence compiler. really testing == -eslINFINITY */
  for (x = 0; x < gm->M*p7P_NTRANS; x++)
      gm->tsc[x] = (gm->tsc[x] <= -eslINFINITY) ? -255 : roundf(om->scale * gm->tsc[x]);
  
  /* Emissions */
  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 0; k <= p7P_NR*gm->M; k++)
      gm->rsc[x][k] = (gm->rsc[x][k] <= -eslINFINITY) ? -255 : roundf(om->scale * gm->rsc[x][k]);

  /* Specials */
  for (k = 0; k < p7P_NXSTATES; k++)
    for (x = 0; x < p7P_NXTRANS; x++)
      gm->xsc[k][x] = (gm->xsc[k][x] <= -eslINFINITY) ? -255 : roundf(om->scale * gm->xsc[k][x]);

  /* NN, CC, JJ are hardcoded 0 in uchar limited precision */
  gm->xsc[p7P_N][p7P_LOOP] = 0;
  gm->xsc[p7P_C][p7P_LOOP] = 0;
  gm->xsc[p7P_J][p7P_LOOP] = 0;

  return eslOK;
}

/* Function:  p7_oprofile_SameMSV()
 * Synopsis:  Set a generic profile's scores to give MSV scores.
 * Incept:    SRE, Wed Jul 30 13:42:49 2008 [Janelia]
 *
 * Purpose:   Set a generic profile's scores so that the normal <dp_generic> DP 
 *            algorithms will give the same score as <p7_MSVFilter()>:
 *            all t_MM scores = 0; all other core transitions = -inf;
 *            multihit local mode; all <t_BMk> entries uniformly <log 2/(M(M+1))>;
 *            <tCC, tNN, tJJ> scores 0; total approximated later as -3;
 *            rounded in the same way as the 8-bit limited precision.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_oprofile_SameMSV(const P7_OPROFILE *om, P7_PROFILE *gm)
{
  int    k;

  esl_vec_FSet(gm->tsc, p7P_NTRANS * gm->M, -eslINFINITY);
  for (k = 1; k <  gm->M; k++) p7P_TSC(gm, k, p7P_MM) = 0.0f;
  for (k = 0; k <  gm->M; k++) p7P_TSC(gm, k, p7P_BM) = log(2.0f / ((float) gm->M * (float) (gm->M+1)));
  
  gm->xsc[p7P_N][p7P_LOOP] =  gm->xsc[p7P_J][p7P_LOOP] =  gm->xsc[p7P_C][p7P_LOOP] = 0;

  return p7_oprofile_SameRounding(om, gm);
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
 * 
 * Note:      For now, it dumps both uchar and float parts of the
 *            model, and the float part is 8.5 formatted for
 *            ForwardFilter's pspace scores. We might eventually
 *            want to specialize a bit more: allow dumping only
 *            one part of the model, or dumping in a 5.2 format
 *            better suited for ViterbiScore()'s lspace scores.
 */
int
p7_oprofile_Dump(FILE *fp, const P7_OPROFILE *om)
{
  int status;

  fprintf(fp, "Dump of a <P7_OPROFILE> ::\n");
  /* Dump float part, ForwardFilter pspace scores, in %8.5 float format */
  fprintf(fp, "  -- float part, odds for ForwardFilter():\n");
  if ((status = oprofile_dump_float(fp, om, 8, 5)) != eslOK) return status;

  /* Dump uchar part, ViterbiFilter lspace scores */
  fprintf(fp, "  -- uchar part, log odds for ViterbiFilter(): \n");
  if ((status = oprofile_dump_uchar(fp, om))       != eslOK) return status;

  return eslOK;
}


/* oprofile_dump_uchar()
 * 
 * Dump the uchar part of a profile <om> to <stdout>.
 */
static int
oprofile_dump_uchar(FILE *fp, const P7_OPROFILE *om)
{
  int     M   = om->M;		/* length of the query                                          */
  int     nq  = p7O_NQU(M);     /* segment length; total # of striped vectors needed            */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     kb;			/* possibly offset base k for loading om's TSC vectors          */
  int     z;			/* counter within elements of one SIMD minivector               */
  int     t;			/* counter over transitions 0..7 = p7O_{BM,MM,IM,DM,MD,MI,II,DD}*/
  int     j;			/* counter in interleaved vector arrays in the profile          */
  union { __m128i v; uint8_t i[16]; } tmp; /* used to align and read simd minivectors           */

  /* Table of residue emissions */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 

      /* Header (rearranged column numbers, in the vectors)  */
      for (k =1, q = 0; q < nq; q++, k++)
	{
	  fprintf(fp, "[ ");
	  for (z = 0; z < 16; z++) 
	    if (k+z*nq <= M) fprintf(fp, "%4d ", k+z*nq);
	    else             fprintf(fp, "%4s ", "xx");
	  fprintf(fp, "]");
	}

      /* Match emission scores */
      fprintf(fp, "\nmat: ");
      for (j = 0, q = 0; q < nq; q++, j+=2)
	{
	  fprintf(fp, "[ ");
	  _mm_store_si128(&tmp.v, om->ru[x][j]);
	  for (z = 0; z < 16; z++) fprintf(fp, "%4d ", tmp.i[z]);
	  fprintf(fp, "]");
	}

      fprintf(fp, "\nins: ");
      for (j = 1, q = 0; q < nq; q++, j+=2)
	{
	  fprintf(fp, "[ ");
	  _mm_store_si128(&tmp.v, om->ru[x][j]);
	  for (z = 0; z < 16; z++) fprintf(fp, "%4d ", tmp.i[z]);
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
	  for (z = 0; z < 16; z++) 
	    if (kb+z*nq <= M) fprintf(fp, "%4d ", kb+z*nq);
	    else              fprintf(fp, "%4s ", "xx");
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n     ");	  
      for (q = 0; q < nq; q++)
	{
	  fprintf(fp, "[ ");
	  _mm_store_si128(&tmp.v, om->tu[q*7 + t]);
	  for (z = 0; z < 16; z++) fprintf(fp, "%4d ", tmp.i[z]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");	  
    }

  /* DD transitions */
  fprintf(fp, "\ntDD: ");
  for (k =1, q = 0; q < nq; q++, k++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < 16; z++) 
	if (k+z*nq <= M) fprintf(fp, "%4d ", k+z*nq);
	else             fprintf(fp, "%4s ", "xx");
      fprintf(fp, "]");
    }
  fprintf(fp, "\n     ");	  
  for (j = nq*7, q = 0; q < nq; q++, j++)
    {
      fprintf(fp, "[ ");
      _mm_store_si128(&tmp.v, om->tu[j]);
      for (z = 0; z < 16; z++) fprintf(fp, "%4d ", tmp.i[z]);
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");	  

  fprintf(fp, "E->C: %4d    E->J: %4d\n", om->xu[p7O_E][p7O_MOVE], om->xu[p7O_E][p7O_LOOP]);
  fprintf(fp, "N->B: %4d    N->N: %4d\n", om->xu[p7O_N][p7O_MOVE], om->xu[p7O_N][p7O_LOOP]);
  fprintf(fp, "J->B: %4d    J->J: %4d\n", om->xu[p7O_J][p7O_MOVE], om->xu[p7O_J][p7O_LOOP]);
  fprintf(fp, "C->T: %4d    C->C: %4d\n", om->xu[p7O_C][p7O_MOVE], om->xu[p7O_C][p7O_LOOP]);

  fprintf(fp, "bound: %4d\n",  om->ddbound_u);
  fprintf(fp, "scale: %.2f\n", om->scale);
  fprintf(fp, "base:  %4d\n",  om->base);
  fprintf(fp, "bias:  %4d\n",  om->bias);
  fprintf(fp, "Q:     %d\n",   nq);  
  fprintf(fp, "M:     %d\n",   M);  
  return eslOK;
}


/* oprofile_dump_float()
 * 
 * Dump the float part of a profile <om> to <stdout>.
 * <width>, <precision> control the floating point output:
 *  8,5 is a reasonable choice for prob space,
 *  5,2 is reasonable for log space.
 */
static int
oprofile_dump_float(FILE *fp, const P7_OPROFILE *om, int width, int precision)
{
  int     M   = om->M;		/* length of the query                                          */
  int     nq  = p7O_NQF(M);     /* segment length; total # of striped vectors needed            */
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
      for (j = 0, q = 0; q < nq; q++, j+=2)
	{
	  fprintf(fp, "[ ");
	  tmp.v = om->rf[x][j];
	  for (z = 0; z < 4; z++) fprintf(fp, "%*.*f ", width, precision, tmp.x[z]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\nins: ");
      for (j = 1, q = 0; q < nq; q++, j+=2)
	{
	  fprintf(fp, "[ ");
	  tmp.v = om->rf[x][j];
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
	  case p7O_BM: kb = k;                 break; 
	  case p7O_MM: kb = 1 + (nq+k-2) % nq; break; /* MM, DM, IM quads rotated by +1  */
	  case p7O_IM: kb = 1 + (nq+k-2) % nq; break;  
	  case p7O_DM: kb = 1 + (nq+k-2) % nq; break;  
	  case p7O_MD: kb = k;                 break; /* the remaining ones are straight up  */
	  case p7O_MI: kb = k;                 break; 
	  case p7O_II: kb = k;                 break; 
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
	  tmp.v = om->tf[q*7 + t];
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
      tmp.v = om->tf[j];
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
}
/*------------ end, debugging dumps of P7_OPROFILE --------------*/


/*****************************************************************
 * 4. Benchmark driver.
 *****************************************************************/

#ifdef p7OPROFILE_BENCHMARK
/* Timing profile conversion.
   gcc -o benchmark-oprofile -std=gnu99 -g -Wall -msse2 -I.. -L.. -I../../easel -L../../easel -Dp7OPROFILE_BENCHMARK\
      p7_oprofile.c -lhmmer -leasel -lm 
   icc -o benchmark-oprofile -O3 -static -I.. -L.. -I../../easel -L../../easel -Dp7OPROFILE_BENCHMARK p7_oprofile.c -lhmmer -leasel -lm 
   ./benchmark-sse <hmmfile>         runs benchmark
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-L",        eslARG_INT,    "400", NULL, NULL,  NULL,  NULL, NULL, "length of target sequence",                        0 },
  { "-N",        eslARG_INT, "100000", NULL, NULL,  NULL,  NULL, NULL, "number of conversions to time",                    0 },
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
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  int             i;

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    p7_oprofile_Convert(gm, om);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M = %d\n", gm->M);

  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7OPROFILE_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/




  
/*****************************************************************
 * 5. Unit tests
 *****************************************************************/
#ifdef p7OPROFILE_TESTDRIVE


#endif /*p7OPROFILE_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/




/*****************************************************************
 * 6. Test driver
 *****************************************************************/
#ifdef p7OPROFILE_TESTDRIVE


#endif /*p7OPROFILE_TESTDRIVE*/
/*------------------- end, test driver --------------------------*/


/*****************************************************************
 * 7. Example
 *****************************************************************/
#ifdef p7OPROFILE_EXAMPLE


#endif /*p7OPROFILE_EXAMPLE*/
/*------------------- end, test driver --------------------------*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
