/* SSE implementation of fast filters, both Viterbi and Forward.
 * 
 * A filter module provides a standard API to the p7_ViterbiFilter()
 * and p7_ForwardFilter() calls. The API requires that the module
 * implement support for two objects, a P7_OPROFILE optimized score
 * profile and a P7_OMX optimized dynamic programming matrix.
 * 
 * Table of contents:
 *   1. The P7_OPROFILE structure: a score profile
 *   2. The P7_OMX structure: a dynamic programming matrix
 *   3. Debugging dumps of P7_OPROFILE structures
 *   4. Debugging dumps of P7_OMX structures
 *   5. Conversions into P7_OPROFILE format
 *   6. MSV filter implementation
 *   7. Viterbi filter DP implementation
 *   8. Forward filter DP implementation
 *   9. Viterbi score DP implementation
 *  10. Benchmark drivers
 *  11. Unit tests
 *  12. Test driver
 *  13. Example
 *  14. Copyright and license information
 * 
 * SRE, Sun Nov 25 11:26:48 2007 [Casa de Gatos]
 * SVN $Id$
 */
#include "p7_config.h"
#if defined (p7_IMPL_VMX)

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include <altivec.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"
#include "esl_vmx.h" 

#include "hmmer.h"
#include "impl_vmx.h"

/* ?MX(q) access macros work for either uchar or float, so long as you
 * init your "dp" to point to the appropriate array.
 */
#define MMX(q) (dp[(q) * p7X_NSCELLS + p7X_M])
#define DMX(q) (dp[(q) * p7X_NSCELLS + p7X_D])
#define IMX(q) (dp[(q) * p7X_NSCELLS + p7X_I])

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
  om->tu = NULL;
  om->ru = NULL;
  om->tf = NULL;
  om->rf = NULL;

  /* level 1 */
  ESL_ALLOC(om->tu_mem, sizeof(vector unsigned char)   * nqu  * p7O_NTRANS + 15);    
  ESL_ALLOC(om->tf_mem, sizeof(vector float)           * nqf  * p7O_NTRANS + 15);    
  ESL_ALLOC(om->ru, sizeof(vector unsigned char *) * abc->Kp); 
  ESL_ALLOC(om->rm, sizeof(vector unsigned char *) * abc->Kp); 
  ESL_ALLOC(om->rf, sizeof(vector float *)         * abc->Kp); 
  om->tu = (vector unsigned char *) ((((size_t) om->tu_mem) + 15) & (~0xf));
  om->tf = (vector float         *) ((((size_t) om->tf_mem) + 15) & (~0xf));
  om->ru[0] = NULL;
  om->rm[0] = NULL;
  om->rf[0] = NULL;

  /* level 2 */
  ESL_ALLOC(om->ru_mem, sizeof(vector unsigned char) * nqu  * p7O_NR * abc->Kp + 15);                     
  ESL_ALLOC(om->rm_mem, sizeof(vector unsigned char) * nqu           * abc->Kp + 15);                     
  ESL_ALLOC(om->rf_mem, sizeof(vector float)         * nqf  * p7O_NR * abc->Kp + 15);                     
  om->ru[0] = (vector unsigned char *) ((((size_t) om->ru_mem) + 15) & (~0xf));
  om->rm[0] = (vector unsigned char *) ((((size_t) om->rm_mem) + 15) & (~0xf));
  om->rf[0] = (vector float         *) ((((size_t) om->rf_mem) + 15) & (~0xf));
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

  om->name      = NULL;
  om->nj        = 0.0f;
  om->mode      = p7_NO_MODE;
  om->allocM    = allocM;
  om->M         = 0;
  om->abc       = abc;
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

  if (om->name   != NULL) free(om->name);
  if (om->tu_mem != NULL) free(om->tu_mem);
  if (om->tf_mem != NULL) free(om->tf_mem);

  if (om->ru != NULL) { if (om->ru_mem != NULL) free(om->ru_mem);  free(om->ru); }
  if (om->rm != NULL) { if (om->rm_mem != NULL) free(om->rm_mem);  free(om->rm); }
  if (om->rf != NULL) { if (om->rf_mem != NULL) free(om->rf_mem);  free(om->rf); }
  free(om);
}
/*----------------- end, P7_OPROFILE structure ------------------*/




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
  P7_OMX  *ox  = NULL;
  int      nqf = p7O_NQF(allocM);	       /* segment length; total # of striped vectors for uchar */
  int      nqu = p7O_NQU(allocM);	       /* segment length; total # of striped vectors for float */
  int      status;

  ESL_ALLOC(ox, sizeof(P7_OMX));
  ox->dpu = NULL;
  ox->dpf = NULL;

  ESL_ALLOC(ox->dpu_mem,  sizeof(vector signed int) * p7X_NSCELLS * nqu + 15); 
  ox->dpu = (vector unsigned char *) ((((size_t) ox->dpu_mem) + 15) & (~0xf));
  ox->allocQ16  = nqu;
  ox->Q16       = 0;

  ESL_ALLOC(ox->dpf_mem,  sizeof(vector float)      * p7X_NSCELLS * nqf + 15);
  ox->dpf = (vector float *) ((((size_t) ox->dpf_mem) + 15) & (~0xf));
  ox->allocQ4   = nqf;
  ox->Q4        = 0;

  ox->allocM    = allocM;
  ox->M         = 0;
#ifdef p7_DEBUGGING
  ox->debugging = FALSE;
  ox->dfp       = NULL;
#endif
  return ox;

 ERROR:
  p7_omx_Destroy(ox);
  return NULL;
}


/* Function:  p7_omx_GrowTo()
 * Synopsis:  Assure that a DP matrix is big enough.
 * Incept:    SRE, Thu Dec 20 09:27:07 2007 [Janelia]
 *
 * Purpose:   Assures that an optimized DP matrix <ox> is allocated for
 *            a model up to <allocM> in length; if not, reallocate to
 *            make it so.
 *            
 *            Because the optimized matrix is one-row, only the model
 *            length matters; the target sequence length isn't
 *            relevant.
 *
 * Returns:   <eslOK> on success, and <gx> may be reallocated upon
 *            return; any data that may have been in <gx> must be 
 *            assumed to be invalidated.
 *
 * Throws:    <eslEMEM> on allocation failure, and any data that may
 *            have been in <gx> must be assumed to be invalidated.
 */
int
p7_omx_GrowTo(P7_OMX *ox, int allocM)
{
  if (allocM <= ox->allocM) return eslOK;

  void *p;
  int   nqf = p7O_NQF(allocM);	       /* segment length; total # of striped vectors for uchar */
  int   nqu = p7O_NQU(allocM);	       /* segment length; total # of striped vectors for float */
  int   status;
 
  ESL_RALLOC(ox->dpu_mem, p, sizeof(vector signed int) * p7X_NSCELLS * nqu + 15);
  ox->dpu = (vector unsigned char *) ((((size_t) ox->dpu_mem) + 15) & (~0xf));
  ox->allocQ16 = nqu;

  ESL_RALLOC(ox->dpf_mem, p, sizeof(vector float)      * p7X_NSCELLS * nqf + 15);
  ox->dpf = (vector float *) ((((size_t) ox->dpf_mem) + 15) & (~0xf));
  ox->allocQ4  = nqf;

  ox->allocM = allocM;
  ox->M      = 0; 
  return     eslOK;

 ERROR:
  return status;
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
  if (ox->dpu_mem != NULL) free(ox->dpu_mem);
  if (ox->dpf_mem != NULL) free(ox->dpf_mem);
  free(ox);
  return;
}
/*------------------- end, P7_OMX structure ---------------------*/




/*****************************************************************
 * 3. Debugging dumps of P7_OPROFILE structures
 *****************************************************************/

/* oprofile_dump_uchar()
 * 
 * Dump the uchar part of a profile <om> to <stdout>.
 */
static int
oprofile_dump_uchar(FILE *fp, P7_OPROFILE *om)
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
  union { vector unsigned char v; uint8_t i[16]; } tmp; /* used to align and read simd minivectors           */

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
	  //_mm_store_si128(&tmp.v, om->ru[x][j]);
          vec_st(om->ru[x][j],0,&tmp.v);
	  for (z = 0; z < 16; z++) fprintf(fp, "%4d ", tmp.i[z]);
	  fprintf(fp, "]");
	}

      fprintf(fp, "\nins: ");
      for (j = 1, q = 0; q < nq; q++, j+=2)
	{
	  fprintf(fp, "[ ");
	  //_mm_store_si128(&tmp.v, om->ru[x][j]);
          vec_st(om->ru[x][j],0,&tmp.v);
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
	  //_mm_store_si128(&tmp.v, om->tu[q*7 + t]);
          vec_st(om->tu[q*7 + t],0,&tmp.v);
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
      //_mm_store_si128(&tmp.v, om->tu[j]);
      vec_st(om->tu[j],0,&tmp.v);
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
oprofile_dump_float(FILE *fp, P7_OPROFILE *om, int width, int precision)
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
  union { vector float v; float x[4]; } tmp; /* used to align and read simd minivectors               */

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
p7_oprofile_Dump(FILE *fp, P7_OPROFILE *om)
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
/*------------ end, debugging dumps of P7_OPROFILE --------------*/




/*****************************************************************
 * 4. Debugging dumps of P7_OMX structures
 *****************************************************************/

/* Because the P7_OMX is a one-row DP matrix, we can't just run a DP
 * calculation and then dump a whole matrix; we have to dump each row
 * one at a time, as the DP calculation is progressing. Thus we need
 * to call the dump from *within* the DP routine. We'd rather not have
 * anything like this in production code - not even a flag check.  So,
 * we use a compile-time debugging idiom, with conditionally compiled
 * debugging code that's added to the DP routines to check a
 * debugging flag in the P7_OMX structure; if it's up, we dump a row.
 *
 * Therefore, the externally exposed API call is p7_omx_SetDumpMode(),
 * rather than the dumping routine itself; and all p7_omx_SetDumpMode()
 * does is sets the debugging flag in <ox>.
 */


/* Function:  p7_omx_SetDumpMode()
 * Synopsis:  Set an optimized DP matrix to be dumped for debugging.
 * Incept:    SRE, Thu Dec 13 10:24:38 2007 [Janelia]
 *
 * Purpose:   Sets debugging mode for DP matrix <ox>.  If <truefalse>
 *            flag is <TRUE>, then whenever a dynamic programming
 *            calculation is run, dump DP matrix <ox> to stream <fp>
 *            for diagnostics.
 *            
 *            When the dump mode is on, the DP routine itself actually
 *            does the dumping, because it has to dump after every row
 *            is calculated. (We're doing an optimized one-row
 *            calculation.)
 *            
 *            If the code has not been compiled with the
 *            <p7_DEBUGGING> flag up, this function is a no-op.
 *
 * Args:      fp        - output stream for diagnostics (stdout, perhaps)
 *            ox        - DP matrix to set debugging mode
 *            truefalse - TRUE to set dumping, FALSE to unset
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      J2/62.
 */
int
p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox, int truefalse)
{
#if p7_DEBUGGING
  ox->debugging = truefalse;
  ox->dfp       = fp;
#endif
  return eslOK;
}


#ifdef p7_DEBUGGING
/* omx_dump_uchar_row()
 *
 * Dump current row of uchar part of DP matrix <ox> for diagnostics,
 * and include the values of specials <xE>, etc. The index <rowi> for
 * the current row is used as a row label.
 * 
 * If <rowi> is 0, print a header first too.
 * 
 * The output format is coordinated with <p7_gmx_Dump()> to
 * facilitate comparison to a known answer.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
static int
omx_dump_uchar_row(P7_OMX *ox, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC)
{
  vector unsigned char *dp = ox->dpu;	/* must set <dp> before using {MDI}MX macros */
  int      Q  = ox->Q16;
  int      M  = ox->M;
  uint8_t *v  = NULL;		/* array of unstriped, uninterleaved scores  */
  int      q,z,k;
  union { vector unsigned char v; uint8_t i[16]; } tmp;
  int      status;

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

  /* Unpack and unstripe, then print M's. */
  for (q = 0; q < Q; q++) {
    tmp.v = MMX(q);
    for (z = 0; z < 16; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d M ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", v[k]);

  /* The specials */
  fprintf(ox->dfp, "%3d %3d %3d %3d %3d\n", xE, xN, xJ, xB, xC);

  /* Unpack and unstripe, then print I's. */
  for (q = 0; q < Q; q++) {
    tmp.v = IMX(q);
    for (z = 0; z < 16; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d I ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", v[k]);
  fprintf(ox->dfp, "\n");

  /* Unpack, unstripe, then print D's. */
  for (q = 0; q < Q; q++) {
    tmp.v = DMX(q);
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

/* omx_dump_float_row()
 *
 * Dump current row of float part of DP matrix <ox> for diagnostics,
 * and include the values of specials <xE>, etc. The index <rowi> for
 * the current row is used as a row label. The output format of the
 * floats is controlled by <width>, <precision>; 8,5 is good for
 * pspace, 5,2 is fine for lspace.
 * 
 * If <rowi> is 0, print a header first too.
 * 
 * If <logify> is TRUE, then scores are printed as log(score); this is
 * useful for comparing DP with pspace scores with other DP matrices
 * (like generic P7_GMX ones) that use log-odds scores.
 * 
 * The output format is coordinated with <p7_gmx_Dump()> to
 * facilitate comparison to a known answer.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
static int
omx_dump_float_row(P7_OMX *ox, int logify, int rowi, int width, int precision, float xE, float xN, float xJ, float xB, float xC)
{
  vector float *dp  = ox->dpf;	/* must set <dp> before using {MDI}MX macros */
  int      Q  = ox->Q4;
  int      M  = ox->M;
  float   *v  = NULL;		/* array of uninterleaved, unstriped scores  */
  int      q,z,k;
  union { vector float v; float x[16]; } tmp;
  int      status;

  ESL_ALLOC(v, sizeof(float) * ((Q*4)+1));
  v[0] = 0.;

  if (rowi == 0)
    {
      fprintf(ox->dfp, "      ");
      for (k = 0; k <= M;  k++) fprintf(ox->dfp, "%*d ", width, k);
      fprintf(ox->dfp, "%*s %*s %*s %*s %*s\n", width, "E", width, "N", width, "J", width, "B", width, "C");
      fprintf(ox->dfp, "      ");
      for (k = 0; k <= M+5;  k++) fprintf(ox->dfp, "%*s ", width, "--------");
      fprintf(ox->dfp, "\n");
    }

  /* Unpack, unstripe, then print M's. */
  for (q = 0; q < Q; q++) {
    tmp.v = MMX(q);
    for (z = 0; z < 4; z++) v[q+Q*z+1] = tmp.x[z];
  }
  fprintf(ox->dfp, "%3d M ", rowi);
  if (logify) for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k] == 0. ? -eslINFINITY : log(v[k]));
  else        for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k]);

 /* The specials */
  if (logify) fprintf(ox->dfp, "%*.*f %*.*f %*.*f %*.*f %*.*f\n",
		      width, precision, xE == 0. ? -eslINFINITY : log(xE),
		      width, precision, xN == 0. ? -eslINFINITY : log(xN),
		      width, precision, xJ == 0. ? -eslINFINITY : log(xJ),
		      width, precision, xB == 0. ? -eslINFINITY : log(xB), 
		      width, precision, xC == 0. ? -eslINFINITY : log(xC));
  else        fprintf(ox->dfp, "%*.*f %*.*f %*.*f %*.*f %*.*f\n",
		      width, precision, xE,   width, precision, xN, width, precision, xJ, 
		      width, precision, xB,   width, precision, xC);

  /* Unpack, unstripe, then print I's. */
  for (q = 0; q < Q; q++) {
    tmp.v = IMX(q);
    for (z = 0; z < 4; z++) v[q+Q*z+1] = tmp.x[z];
  }
  fprintf(ox->dfp, "%3d I ", rowi);
  if (logify) for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k] == 0. ? -eslINFINITY : log(v[k]));
  else        for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k]);
  fprintf(ox->dfp, "\n");

  /* Unpack, unstripe, then print D's. */
  for (q = 0; q < Q; q++) {
    tmp.v = DMX(q);
    for (z = 0; z < 4; z++) v[q+Q*z+1] = tmp.x[z];
  }
  fprintf(ox->dfp, "%3d D ", rowi);
  if (logify) for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k] == 0. ? -eslINFINITY : log(v[k]));
  else        for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k]);
  fprintf(ox->dfp, "\n\n");

  free(v);
  return eslOK;

ERROR:
  free(v);
  return status;
}

/* omx_dump_msv_row()
 *
 * Dump current row of uchar part of DP matrix <ox> for diagnostics,
 * and include the values of specials <xE>, etc. The index <rowi> for
 * the current row is used as a row label. This routine has to be
 * specialized for the layout of the MSVFilter() row, because it's
 * all match scores dp[0..q..Q-1], rather than triplets of M,D,I.
 *
* If <rowi> is 0, print a header first too.
 *
 * The output format is coordinated with <p7_gmx_Dump()> to
 * facilitate comparison to a known answer.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
static int
omx_dump_msv_row(P7_OMX *ox, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC)
{
   vector unsigned char *dp = ox->dpu;
   int      Q  = ox->Q16;
   int      M  = ox->M;
   uint8_t *v  = NULL;          /* array of unstriped scores  */
   int      q,z,k;
   union { vector unsigned char v; uint8_t i[16]; } tmp;
   int      status;
 
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
 
   /* Unpack and unstripe, then print M's. */
   for (q = 0; q < Q; q++) {
     tmp.v = dp[q];
     for (z = 0; z < 16; z++) v[q+Q*z+1] = tmp.i[z];
   }
   fprintf(ox->dfp, "%4d M ", rowi);
   for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", v[k]);
 
   /* The specials */
   fprintf(ox->dfp, "%3d %3d %3d %3d %3d\n", xE, xN, xJ, xB, xC);
 
   /* I's are all 0's; print just to facilitate comparison. */
   fprintf(ox->dfp, "%4d I ", rowi);
   for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", 0);
   fprintf(ox->dfp, "\n");
 
   /* D's are all 0's too */
   fprintf(ox->dfp, "%4d D ", rowi);
   for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", 0);
   fprintf(ox->dfp, "\n\n");
 
   free(v);
   return eslOK;
 
ERROR:
   free(v);
   return status;
}
#endif /*p7_DEBUGGING*/
/*------------- end, debugging dumps of P7_OMX ------------------*/




/*****************************************************************
 * 5. Conversions into P7_OPROFILE format
 *****************************************************************/

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

  sc  = -1.0f * roundf(om->scale * sc);        	        /* ugh. sc is now an integer cost represented in a float...    */
  b   = (sc > 255.) ? 255 : (uint8_t) sc + om->bias;	/* and now we cast and saturate it to an unsigned char cost... */
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
lspace_uchar_conversion(P7_PROFILE *gm, P7_OPROFILE *om)
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
  union { vector unsigned char v; uint8_t i[16]; } tmp; /* used to align and load simd minivectors           */

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
lspace_float_conversion(P7_PROFILE *gm, P7_OPROFILE *om)
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
  union { vector float v; float x[4]; } tmp; /* used to align and load simd minivectors               */

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

/* pspace_to_lspace_float() 
 * 
 * Take the log of all the scores. 
 * This is needed to convert a ForwardFilter() model into
 * one that can call ViterbiScore(), for example.
 */
static int
pspace_to_lspace_float(P7_OPROFILE *om)
{
  int nq  = p7O_NQF(om->M);   /* segment length; total # of striped vectors needed            */
  int x;
  int j;

  for (x = 0; x < om->abc->Kp; x++)
    for (j = 0; j < nq*2; j++)
      om->rf[x][j] = esl_vmx_logf(om->rf[x][j]);

  for (j = 0; j < nq*p7O_NTRANS; j++)
    om->tf[j] = esl_vmx_logf(om->tf[j]);

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

/* lspace_to_pspace_float() 
 * 
 * Exponentiate all the scores. 
 * This is needed to convert a ViterbiScore() model
 * back to one that can be used in ForwardFilter(),
 * for example.
 */
static int
lspace_to_pspace_float(P7_OPROFILE *om)
{
  int nq  = p7O_NQF(om->M);   /* segment length; total # of striped vectors needed            */
  int x;
  int j;

  for (x = 0; x < om->abc->Kp; x++)
    for (j = 0; j < nq*2; j++)
      om->rf[x][j] = esl_vmx_expf(om->rf[x][j]);

  for (j = 0; j < nq*p7O_NTRANS; j++)
    om->tf[j] = esl_vmx_expf(om->tf[j]);

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
pspace_float_conversion(P7_PROFILE *gm, P7_OPROFILE *om)
{
  int     nq  = p7O_NQF(gm->M);  /* segment length; total # of striped vectors needed            */
  int     status;

  if (nq > om->allocQ4) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  /* First stripe and interleave the scores from <gm> */
  if ((status = lspace_float_conversion(gm, om)) != eslOK) return status;

  /* Then exponentiate them all quickly and in stride (using SSE) */
  if ((status = lspace_to_pspace_float(om))      != eslOK) return status;

  om->lspace_f = FALSE;
  om->mode     = gm->mode;
  om->M        = gm->M;
  return eslOK;
}





/* Function:  p7_oprofile_Convert()
 * Synopsis:  Converts standard profile to an optimized one.
 * Incept:    SRE, Mon Nov 26 07:38:57 2007 [Janelia]
 */
int
p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om)
{
  int status, z;

  if (gm->abc->type != om->abc->type)  ESL_EXCEPTION(eslEINVAL, "alphabets of the two profiles don't match");  

  if ((status =  lspace_uchar_conversion(gm, om)) != eslOK) return status;   /* ViterbiFilter()'s information */
  if ((status =  pspace_float_conversion(gm, om)) != eslOK) return status;   /* ForwardFilter()'s information */

  if ((status = esl_strdup(gm->name, -1, &(om->name))) != eslOK) goto ERROR;
  for (z = 0; z < p7_NEVPARAM; z++) om->evparam[z] = gm->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) om->cutoff[z]  = gm->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) om->compo[z]   = gm->compo[z];

  /* MSVFilter's constants */
  om->tbm  = unbiased_charify(om, logf(2.0f / ((float) gm->M * (float) (gm->M+1)))); /* constant B->Mk penalty        */
  om->tec  = unbiased_charify(om, logf(0.5f));                                       /* constant multihit E->C = E->J */
  /* tjb is length dependent; see ReconfigLength() for setting, below */

  om->nj   = gm->nj;
  om->mode = gm->mode;
  om->M    = gm->M;
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
  return eslOK;
}
/*------------ end, conversions to P7_OPROFILE ------------------*/





/*****************************************************************
 * 6. The p7_MSVFilter() DP implementation.
 *****************************************************************/

/* Returns TRUE if any a[z] > b[z] .
 * This is an incantation. SSE provides no cmpgt_epu8 instruction!
 * Note that cmpeq_epi8 works fine for unsigned ints (there is no 
 * cmpeq_epu8 instruction either). 
 */
//static int 
//sse_any_gt_epu8(__m128i a, __m128i b)
//{
//  __m128i mask    = _mm_cmpeq_epi8(_mm_max_epu8(a,b), b); /* anywhere a>b, mask[z] = 0x0; elsewhere 0xff */
//  int   maskbits  = _mm_movemask_epi8(_mm_xor_si128(mask,  _mm_cmpeq_epi8(mask, mask)));
//  return maskbits != 0;
//}
// should be replaced by Altivec predicate vec_any_gt()

/* Returns maximum element \max_z a[z] in epu8 vector */
static unsigned char
vmx_hmax_vecuchar(vector unsigned char a)
{
  union { vector unsigned char v; uint8_t i[16]; } tmp;

  vector unsigned char   onevec = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  vector unsigned char shiftvec = {1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3,
                                   1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3};
  
  tmp.v = vec_max(a    , vec_slo(a    , shiftvec));
  shiftvec = vec_sl(shiftvec, onevec);
  tmp.v = vec_max(tmp.v, vec_slo(tmp.v, shiftvec));
  shiftvec = vec_sl(shiftvec, onevec);
  tmp.v = vec_max(tmp.v, vec_slo(tmp.v, shiftvec));
  shiftvec = vec_sl(shiftvec, onevec);
  tmp.v = vec_max(tmp.v, vec_slo(tmp.v, shiftvec));
  return tmp.i[0];
}



/* Function:  p7_MSVFilter()
 * Synopsis:  Calculates MSV score, vewy vewy fast, in limited precision.
 * Incept:    SRE, Wed Dec 26 15:12:25 2007 [Janelia]
 *
 * Purpose:   Calculates an approximation of the MSV score for sequence
 *            <dsq> of length <L> residues, using optimized profile <om>,
 *            and a preallocated one-row DP matrix <ox>. Return the 
 *            estimated MSV score (in nats) in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow. 
 *            
 *            The model may be in any mode, because only its match
 *            emission scores will be used. The MSV filter inherently
 *            assumes a multihit local mode, and uses its own special
 *            state transition scores, not the scores in the profile.
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - DP matrix
 *            ret_sc  - RETURN: Viterbi score (in nats)          
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_MSVFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register vector unsigned char mpv;       /* previous row values                                       */
  register vector unsigned char xEv;	   /* E state: keeps max for Mk->E as we go                     */
  register vector unsigned char xBv;	   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register vector unsigned char sv;	   /* temp storage of 1 curr row value in progress              */
  register vector unsigned char biasv;	   /* emission bias in a vector                                 */
  uint8_t  xE, xB, xC;             	   /* special states' scores                                    */
  int i;			   	   /* counter over sequence positions 1..L                      */
  int q;			   	   /* counter over vectors 0..nq-1                              */
  int Q        = p7O_NQU(om->M);   	   /* segment length: # of vectors                              */
  vector unsigned char *dp  = ox->dpu;	   /* we're going to use dp[0..q..Q-1], not {MDI}MX(q) macros   */
  vector unsigned char *rsc;		   /* will point at om->ru[x] for residue x[i]                  */

  vector unsigned char shiftvec = {1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3,
                                   1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3};

  /* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ16)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M   = om->M;
  ox->Q16 = Q;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base.
   */
  biasv = (vector unsigned char) {(int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias,
                                  (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias,
                                  (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias,
                                  (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias};
  for (q = 0; q < Q; q++)
    dp[q] = (vector unsigned char) {0};
  xB   = om->base - om->tjb;                /* remember, all values are costs to be subtracted. */
  xC   = 0;

#if p7_DEBUGGING
  if (ox->debugging) omx_dump_msv_row(ox, 0, 0, 0, xC, xB, xC);
#endif

  for (i = 1; i <= L; i++)
    {
      rsc = om->rm[dsq[i]];
      xEv = (vector unsigned char) {0};      
      xBv = (vector unsigned char) {(int8_t) (xB - om->tbm), (int8_t) (xB - om->tbm), (int8_t) (xB - om->tbm), (int8_t) (xB - om->tbm),
                                    (int8_t) (xB - om->tbm), (int8_t) (xB - om->tbm), (int8_t) (xB - om->tbm), (int8_t) (xB - om->tbm),
                                    (int8_t) (xB - om->tbm), (int8_t) (xB - om->tbm), (int8_t) (xB - om->tbm), (int8_t) (xB - om->tbm),
                                    (int8_t) (xB - om->tbm), (int8_t) (xB - om->tbm), (int8_t) (xB - om->tbm), (int8_t) (xB - om->tbm)};

      /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically, which is our -infinity.
       */
      //mpv = _mm_slli_si128(dp[Q-1], 1);   
      mpv = vec_sro(dp[Q-1], shiftvec);
      for (q = 0; q < Q; q++)
	{
	  /* Calculate new MMX(i,q); don't store it yet, hold it in sv. */
          sv   = vec_max(mpv, xBv);
          sv   = vec_adds(sv, biasv);
          sv   = vec_subs(sv, *rsc); rsc++;
          xEv  = vec_max(xEv, sv);

	  mpv   = dp[q];   	  /* Load {MDI}(i-1,q) into mpv */
	  dp[q] = sv;       	  /* Do delayed store of M(i,q) now that memory is usable */
	}	  

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      xE = vmx_hmax_vecuchar(xEv);
      if (xE >= 255 - om->bias) { *ret_sc = eslINFINITY; return eslOK; }	/* immediately detect overflow */

      xC = ESL_MAX(xC,        xE  - om->tec);
      xB = ESL_MAX(om->base,  xC) - om->tjb;
	  
#if p7_DEBUGGING
      if (ox->debugging) omx_dump_msv_row(ox, i, xE, 0, xC, xB, xC);   
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T, and add our missing precision on the NN,CC,JJ back */
  *ret_sc = ((float) (xC - om->tjb) - (float) om->base);
  *ret_sc /= om->scale;
  *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
  return eslOK;
}
/*------------------ end, p7_MSVFilter() ------------------------*/






/*****************************************************************
 * 7. The p7_ViterbiFilter() DP implementation.
 *****************************************************************/

/* Function:  p7_ViterbiFilter()
 * Synopsis:  Calculates Viterbi score, vewy vewy fast, in limited precision.
 * Incept:    SRE, Tue Nov 27 09:15:24 2007 [Janelia]
 *
 * Purpose:   Calculates an approximation of the Viterbi score for sequence
 *            <dsq> of length <L> residues, using optimized profile <om>,
 *            and a preallocated one-row DP matrix <ox>. Return the 
 *            estimated Viterbi score (in nats) in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow. 
 *            
 *            The model must be in a local alignment mode; other modes
 *            cannot provide the necessary guarantee of no underflow.
 *            
 *            This is a striped SIMD Viterbi implementation using Intel
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
  register vector unsigned char mpv, dpv, ipv;  /* previous row values                                       */
  register vector unsigned char sv;		   /* temp storage of 1 curr row value in progress              */
  register vector unsigned char dcv;		   /* delayed storage of D(i,q+1)                               */
  register vector unsigned char xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register vector unsigned char xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register vector unsigned char Dmaxv;          /* keeps track of maximum D cell on row                      */
  vector unsigned char  biasv;		   /* emission bias in a vector                                 */
  uint8_t  xE, xB, xC, xJ;	   /* special states' scores                                    */
  uint8_t  Dmax;		   /* maximum D cell score on row                               */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = p7O_NQU(om->M);   /* segment length: # of vectors                              */
  vector unsigned char *dp  = ox->dpu;	   /* using {MDI}MX(q) macro requires initialization of <dp>    */
  vector unsigned char *rsc;			   /* will point at om->ru[x] for residue x[i]                  */
  vector unsigned char *tsc;			   /* will point into (and step thru) om->tu                    */

  vector unsigned char shiftvec = { 1 << 3,  1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3,
                                    1 << 3,  1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3, 1 << 3};

  /* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ16)                                ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  if (om->mode != p7_LOCAL && om->mode != p7_UNILOCAL) ESL_EXCEPTION(eslEINVAL, "Fast filter only works for local alignment");
  ox->M   = om->M;
  ox->Q16 = Q;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base.
   */
  biasv = (vector unsigned char) {(int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias,
                                  (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias,
                                  (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias,
                                  (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias, (int8_t) om->bias};
  for (q = 0; q < Q; q++)
    MMX(q) = IMX(q) = DMX(q) = (vector unsigned char) {0};
  xB   = om->base - om->xu[p7O_N][p7O_MOVE]; /* remember, all values are costs to be subtracted. */
  xJ   = 0;
  xC   = 0;
  xE   = 0;

#if p7_DEBUGGING
  if (ox->debugging) omx_dump_uchar_row(ox, 0, xE, 0, xJ, xB, xC); /* first 0 is <rowi>: do header. second 0 is xN: always 0 here. */
#endif

  for (i = 1; i <= L; i++)
    {
      rsc   = om->ru[dsq[i]];
      tsc   = om->tu;
      dcv   = (vector unsigned char) {0};      /* "-infinity" */
      xEv   = (vector unsigned char) {0};     
      Dmaxv = (vector unsigned char) {0};     
      xBv   = (vector unsigned char) {(char) xB, (char) xB, (char) xB, (char) xB, (char) xB, (char) xB, (char) xB, (char) xB,
                                      (char) xB, (char) xB, (char) xB, (char) xB, (char) xB, (char) xB, (char) xB, (char) xB};

      /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically, which is our -infinity.
       */
      //mpv = MMX(Q-1);  mpv = _mm_slli_si128(mpv, 1);  
      //dpv = DMX(Q-1);  dpv = _mm_slli_si128(dpv, 1);  
      //ipv = IMX(Q-1);  ipv = _mm_slli_si128(ipv, 1);  
      mpv = MMX(Q-1);  mpv = vec_sro(mpv, shiftvec);  
      dpv = DMX(Q-1);  dpv = vec_sro(dpv, shiftvec);  
      ipv = IMX(Q-1);  ipv = vec_sro(ipv, shiftvec);  

      for (q = 0; q < Q; q++)
	{
	  /* Calculate new MMX(i,q); don't store it yet, hold it in sv. */
	  sv   =              vec_subs(xBv, *tsc);  tsc++;
	  sv   = vec_max (sv, vec_subs(mpv, *tsc)); tsc++;
	  sv   = vec_max (sv, vec_subs(ipv, *tsc)); tsc++;
	  sv   = vec_max (sv, vec_subs(dpv, *tsc)); tsc++;
	  sv   = vec_adds(sv, biasv);     
	  sv   = vec_subs(sv, *rsc);                     rsc++;
	  xEv  = vec_max(xEv, sv);
	  
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
	  dcv   = vec_subs(sv, *tsc);  tsc++;
	  Dmaxv = vec_max(dcv, Dmaxv);

	  /* Calculate and store I(i,q) */
	  sv     =                   vec_subs(mpv, *tsc);  tsc++;
	  sv     = vec_max (sv, vec_subs(ipv, *tsc)); tsc++;
	  sv     = vec_adds(sv, biasv);
	  IMX(q) = vec_subs(sv, *rsc);                     rsc++;
	}	  

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      xE = vmx_hmax_vecuchar(xEv);
      if (xE >= 255 - om->bias) { *ret_sc = eslINFINITY; return eslOK; }	/* immediately detect overflow */
      xC = ESL_MAX(xC, xE - om->xu[p7O_E][p7O_MOVE]);
      xJ = ESL_MAX(xJ, xE - om->xu[p7O_E][p7O_LOOP]);
      xB = ESL_MAX(xJ - om->xu[p7O_J][p7O_MOVE],  om->base - om->xu[p7O_N][p7O_MOVE]);
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
      Dmax = vmx_hmax_vecuchar(Dmaxv);
      if ((int) Dmax + om->ddbound_u > (int) xB) 
	{
	  /* Now we're obligated to do at least one complete DD path to be sure. */
	  /* dcv has carried through from end of q loop above */
	  //dcv = _mm_slli_si128(dcv, 1);
          dcv = vec_sro(dcv, shiftvec);
	  tsc = om->tu + 7*Q;	/* set tsc to start of the DD's */
	  for (q = 0; q < Q; q++) 
	    {
	      DMX(q) = vec_max(dcv, DMX(q));	
	      dcv    = vec_subs(DMX(q), *tsc); tsc++;
	    }

	  /* We may have to do up to three more passes; the check
	   * is for whether crossing a segment boundary can improve
	   * our score. 
	   */
	  do {
	    //dcv = _mm_slli_si128(dcv, 1);
            dcv = vec_sro(dcv, shiftvec);
	    tsc = om->tu + 7*Q;	/* set tsc to start of the DD's */
	    for (q = 0; q < Q; q++) 
	      {
		if (! vec_any_gt(dcv, DMX(q))) break;
		DMX(q) = vec_max(dcv, DMX(q));	
		dcv    = vec_subs(DMX(q), *tsc);   tsc++;
	      }	    
	  } while (q == Q);
	}
      else  /* not calculating DD? then just store the last M->D vector calc'ed.*/
	//DMX(0) = _mm_slli_si128(dcv, 1);
        DMX(0) = vec_sro(dcv, shiftvec);
	  
#if p7_DEBUGGING
      if (ox->debugging) omx_dump_uchar_row(ox, i, xE, 0, xJ, xB, xC);   
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T, and add our missing precision on the NN,CC,JJ back */
  *ret_sc = ((float) (xC - om->xu[p7O_C][p7O_MOVE]) - (float) om->base);
  *ret_sc /= om->scale;
  if      (om->mode == p7_UNILOCAL) *ret_sc -= 2.0; /* that's ~ L \log \frac{L}{L+2}, for our NN,CC,JJ */
  else if (om->mode == p7_LOCAL)    *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
  return eslOK;
}
/*---------------- end, p7_ViterbiFilter() ----------------------*/



/*****************************************************************
 * 8. The p7_ForwardFilter() DP implementation.
 *****************************************************************/


/* Function:  p7_ForwardFilter()
 * Synopsis:  Calculates Forward score, vewy vewy fast, with limited upper range.
 * Incept:    SRE, Thu Dec 13 08:54:07 2007 [Janelia]
 *
 * Purpose:   Calculates the Forward score for sequence <dsq> of length <L> 
 *            residues, using optimized profile <om>, and a preallocated
 *            one-row DP matrix <ox>. Return the Forward score (in nats)
 *            in <ret_sc>.
 *            
 *            The Forward score may overflow, and will, on
 *            high-scoring sequences. Range is limited to -88 to +88
 *            nats (-127 to 127 bits). Scores will not underflow, for
 *            models configured in local mode, within HMMER's design
 *            limits ($L \leq 10^5$; $M \leq 10^4$).
 *            
 *            The model must be in a local mode; other modes cannot
 *            guarantee that we will not underflow.
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - DP matrix
 *            ret_sc  - RETURN: Forward score (in nats)          
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small, or if the profile
 *            isn't in local alignment mode.
 */
int
p7_ForwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register vector float mpv, dpv, ipv;   /* previous row values                                       */
  register vector float sv;		   /* temp storage of 1 curr row value in progress              */
  register vector float dcv;		   /* delayed storage of D(i,q+1)                               */
  register vector float xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register vector float xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  vector float   zerov;		   /* splatted 0.0's in a vector                                */
  float    xN, xE, xB, xC, xJ;	   /* special states' scores                                    */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over quads 0..nq-1                                */
  int j;			   /* counter over DD iterations (4 is full serialization)      */
  int Q       = p7O_NQF(om->M);	   /* segment length: # of vectors                              */
  vector float *dp  = ox->dpf;           /* using {MDI}MX(q) macro requires initialization of <dp>    */
  vector float *rp;			   /* will point at om->rf[x] for residue x[i]                  */
  vector float *tp;			   /* will point into (and step thru) om->tf                    */

  vector unsigned char shiftvec = { 4 << 3,  4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3,
                                    4 << 3,  4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3};


/* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ4)                                 ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  if (om->mode != p7_LOCAL && om->mode != p7_UNILOCAL) ESL_EXCEPTION(eslEINVAL, "Fast filter only works for local alignment");
  ox->M  = om->M;
  ox->Q4 = Q;

  /* Initialization.
   */
  zerov = (vector float) {0.0};
  for (q = 0; q < Q; q++)
    MMX(q) = IMX(q) = DMX(q) = zerov;
  xE    = 0.;
  xN    = 1.;
  xJ    = 0.;
  xB    = om->xf[p7O_N][p7O_MOVE];
  xC    = 0.;

#if p7_DEBUGGING
  if (ox->debugging) omx_dump_float_row(ox, TRUE, 0, 9, 5, xE, xN, xJ, xB, xC);	/* logify=TRUE, <rowi>=0, width=8, precision=5*/
#endif

  for (i = 1; i <= L; i++)
    {
      rp    = om->rf[dsq[i]];
      tp    = om->tf;
      dcv   = (vector float) {0.0};
      xEv   = (vector float) {0.0};
      xBv   = (vector float) {xB, xB, xB, xB};

      /* Right shifts by 4 bytes. 4,8,12,x becomes x,4,8,12.  Shift zeros on.
       */
      /* No, your OTHER right!  Literally, the SSE code shifts the vector left
         so the VMX code should _actually_ go right ... I think.  Maybe. */
      //mpv = MMX(Q-1);  mpv = _mm_shuffle_ps(mpv, mpv, _MM_SHUFFLE(2, 1, 0, 0));   mpv = _mm_move_ss(mpv, zerov);
      //dpv = DMX(Q-1);  dpv = _mm_shuffle_ps(dpv, dpv, _MM_SHUFFLE(2, 1, 0, 0));   dpv = _mm_move_ss(dpv, zerov);
      //ipv = IMX(Q-1);  ipv = _mm_shuffle_ps(ipv, ipv, _MM_SHUFFLE(2, 1, 0, 0));   ipv = _mm_move_ss(ipv, zerov);
      mpv = MMX(Q-1);  mpv = vec_sro(mpv, shiftvec);
      dpv = DMX(Q-1);  dpv = vec_sro(dpv, shiftvec);
      ipv = IMX(Q-1);  ipv = vec_sro(ipv, shiftvec); 
      
      for (q = 0; q < Q; q++)
	{
	  /* Calculate new MMX(i,q); don't store it yet, hold it in sv. */
	  sv   =             vec_madd(xBv, *tp, zerov);  tp++;
	  sv   = vec_add(sv, vec_madd(mpv, *tp, zerov)); tp++;
	  sv   = vec_add(sv, vec_madd(ipv, *tp, zerov)); tp++;
	  sv   = vec_add(sv, vec_madd(dpv, *tp, zerov)); tp++;
	  sv   = vec_madd(sv, *rp, zerov);               rp++;
	  xEv  = vec_add(xEv, sv);
	  
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
	  dcv   = vec_madd(sv, *tp, zerov); tp++;

	  /* Calculate and store I(i,q) */
	  sv     =             vec_madd(mpv, *tp, zerov);  tp++;
	  sv     = vec_add(sv, vec_madd(ipv, *tp, zerov)); tp++;
	  IMX(q) = vec_madd(sv, *rp, zerov);               rp++;
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
      //dcv    = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
      //dcv    = _mm_move_ss(dcv, zerov);
      dcv = vec_sro(dcv, shiftvec);
      DMX(0) = zerov;
      tp     = om->tf + 7*Q;	/* set tp to start of the DD's */
      for (q = 0; q < Q; q++) 
	{
	  DMX(q) = vec_add(dcv, DMX(q));	
	  dcv    = vec_madd(DMX(q), *tp, zerov); tp++; /* extend DMX(q), so we include M->D and D->D paths */
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
	  for (j = 1; j < 4; j++)
	    {
	      //dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
	      //dcv = _mm_move_ss(dcv, zerov);
              dcv = vec_sro(dcv, shiftvec);
	      tp = om->tf + 7*Q;	/* set tp to start of the DD's */
	      for (q = 0; q < Q; q++) 
		{
		  DMX(q) = vec_add(dcv, DMX(q));	
		  dcv    = vec_madd(dcv, *tp, zerov);   tp++; /* note, extend dcv, not DMX(q); only adding DD paths now */
		}	    
	    }
	} 
      else
	{			/* Slightly parallelized version, but which incurs some overhead */
	  for (j = 1; j < 4; j++)
	    {
	      //register __m128 cv;	/* keeps track of whether any DD's change DMX(q) */
              register vector float cv;

	      //dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
	      //dcv = _mm_move_ss(dcv, zerov);
              dcv = vec_sro(dcv, shiftvec);
	      tp  = om->tf + 7*Q;	/* set tp to start of the DD's */
	      cv  = zerov;
	      for (q = 0; q < Q; q++) 
		{
		  sv     = vec_add(dcv, DMX(q));	
		  cv     = vec_or(cv, vec_cmpgt(sv, DMX(q))); /* remember if DD paths changed any DMX(q): *without* conditional branch */
		  DMX(q) = sv;	                                    /* store new DMX(q) */
		  dcv    = vec_madd(dcv, *tp, zerov);   tp++;            /* note, extend dcv, not DMX(q); only adding DD paths now */
		}	    
              // FIXME?  not sure about this one...
	      //if (! _mm_movemask_ps(cv)) break; /* DD's didn't change any DMX(q)? Then we're done, break out. */
              if (! vec_any_lt(cv, (vector float) {0.0})) break;
	    }
	}

      /* Add D's to xEv */
      for (q = 0; q < Q; q++) xEv = vec_add(DMX(q), xEv);

      /* Finally the "special" states, which start from Mk->E (->C, ->J->B) */
      /* The following incantation is a horizontal sum of xEv's elements  */
      /* These must follow DD calculations, because D's contribute to E in Forward
       * (as opposed to Viterbi)
       */
      xEv = vec_add(xEv, vec_perm(xEv, xEv, ((vector unsigned char) { 4, 5, 6, 7,  8, 9,10,11, 12,13,14,15,  0, 1, 2, 3})));
      xEv = vec_add(xEv, vec_perm(xEv, xEv, ((vector unsigned char) { 8, 9,10,11, 12,13,14,15,  0, 1, 2, 3,  4, 5, 6, 7})));
      //_mm_store_ss(&xE, xEv);
      vec_ste(xEv, 0, &xE);

      xN =  xN * om->xf[p7O_N][p7O_LOOP];
      xC = (xC * om->xf[p7O_C][p7O_LOOP]) +  (xE * om->xf[p7O_E][p7O_MOVE]);
      xJ = (xJ * om->xf[p7O_J][p7O_LOOP]) +  (xE * om->xf[p7O_E][p7O_LOOP]);
      xB = (xJ * om->xf[p7O_J][p7O_MOVE]) +  (xN * om->xf[p7O_N][p7O_MOVE]);
      /* and now xB will carry over into next i, and xC carries over after i=L */

#if p7_DEBUGGING
      if (ox->debugging) omx_dump_float_row(ox, TRUE, i, 9, 5, xE, xN, xJ, xB, xC);	/* logify=TRUE, <rowi>=i, width=8, precision=5*/
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T, and flip back to log space (nats) */
  /* On overflow, xC is inf or nan (nan arises because inf*0 = nan). */
  /* On an underflow (which shouldn't happen), we counterintuitively return infinity:
   * the effect of this is to force the caller to rescore us with full range.
   */
  if       (isnan(xC))      *ret_sc = eslINFINITY;
  else if  (xC == 0.0)      *ret_sc = eslINFINITY; /* on underflow, force caller to rescore us! */
  else if  (isinf(xC) == 1) *ret_sc = eslINFINITY;
  else                      *ret_sc = log(xC * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}
/*------------------ end, p7_ForwardFilter() --------------------*/



/*****************************************************************
 * 9. Viterbi score DP implementation
 *****************************************************************/

/* Return TRUE if any a[i] > b[i]; from Apple's Altivec/SSE migration guide */
//static int 
//sse_any_gt_ps(__m128 a, __m128 b)
//{
//  __m128 mask    = _mm_cmpgt_ps(a,b);
//  int   maskbits = _mm_movemask_ps( mask );
//  return maskbits != 0;
//}
// Replace with VMX internal vec_any_gt()

/* Function:  p7_ViterbiScore()
 * Synopsis:  Calculates Viterbi score, correctly, and vewy vewy fast.
 * Incept:    SRE, Tue Nov 27 09:15:24 2007 [Janelia]
 *
 * Purpose:   Calculates the Viterbi score for sequence <dsq> of length <L> 
 *            residues, using optimized profile <om>, and a preallocated
 *            one-row DP matrix <ox>. Return the Viterbi score (in nats)
 *            in <ret_sc>.
 *            
 *            The model <om> must be configured specially to have
 *            lspace float scores, not its usual pspace float scores for
 *            <p7_ForwardFilter()>.
 *            
 *            As with all <*Score()> implementations, the score is
 *            accurate (full range and precision) and can be
 *            calculated on models in any mode, not only local modes.
 *            
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - DP matrix
 *            ret_sc  - RETURN: Viterbi score (in nats)          
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_ViterbiScore(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register vector float mpv, dpv, ipv;   /* previous row values                                       */
  register vector float sv;		 /* temp storage of 1 curr row value in progress              */
  register vector float dcv;		 /* delayed storage of D(i,q+1)                               */
  register vector float xEv;		 /* E state: keeps max for Mk->E as we go                     */
  register vector float xBv;		 /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register vector float Dmaxv;           /* keeps track of maximum D cell on row                      */
  vector float  infv;			 /* -eslINFINITY in a vector                                  */
  float    xN, xE, xB, xC, xJ;	   	 /* special states' scores                                    */
  float    Dmax;		   	 /* maximum D cell on row                                     */
  int i;			   	 /* counter over sequence positions 1..L                      */
  int q;			   	 /* counter over vectors 0..nq-1                              */
  int Q       = p7O_NQF(om->M);	   	 /* segment length: # of vectors                              */
  vector float *dp  = ox->dpf;	   	 /* using {MDI}MX(q) macro requires initialization of <dp>    */
  vector float *rsc;			 /* will point at om->rf[x] for residue x[i]                  */
  vector float *tsc;			 /* will point into (and step thru) om->tf                    */

  //vector unsigned char shiftvec = { 4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3,
  //                                  4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3, 4 << 3};
  vector unsigned char selectvec= { 0x10, 0x11, 0x12, 0x13,  0x00, 0x01, 0x02, 0x03,
                                    0x04, 0x05, 0x06, 0x07,  0x08, 0x09, 0x0A, 0x0B};

  /* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ4) ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M  = om->M;
  ox->Q4 = Q;

  /* Initialization. */
  infv = (vector float) {-eslINFINITY, -eslINFINITY, -eslINFINITY, -eslINFINITY};
  for (q = 0; q < Q; q++)
    MMX(q) = IMX(q) = DMX(q) = infv;
  xN   = 0.;
  xB   = om->xf[p7O_N][p7O_MOVE];
  xE   = -eslINFINITY;
  xJ   = -eslINFINITY;
  xC   = -eslINFINITY;

#if p7_DEBUGGING
  if (ox->debugging) omx_dump_float_row(ox, FALSE, 0, 5, 2, xE, xN, xJ, xB, xC); /* logify=FALSE, <rowi>=0, width=5, precision=2*/
#endif

  for (i = 1; i <= L; i++)
    {
      rsc   = om->rf[dsq[i]];
      tsc   = om->tf;
      dcv   = infv;
      xEv   = infv;
      Dmaxv = infv;
      xBv   = (vector float) {xB, xB, xB, xB};

      /* Right shifts by 4 bytes. 4,8,12,x becomes x,4,8,12. 
       */
      //mpv = MMX(Q-1);  mpv = _mm_shuffle_ps(mpv, mpv, _MM_SHUFFLE(2, 1, 0, 0));   mpv = _mm_move_ss(mpv, infv);
      //dpv = DMX(Q-1);  dpv = _mm_shuffle_ps(dpv, dpv, _MM_SHUFFLE(2, 1, 0, 0));   dpv = _mm_move_ss(dpv, infv);
      //ipv = IMX(Q-1);  ipv = _mm_shuffle_ps(ipv, ipv, _MM_SHUFFLE(2, 1, 0, 0));   ipv = _mm_move_ss(ipv, infv);
      mpv = MMX(Q-1);  mpv = vec_perm(mpv, infv, selectvec);
      dpv = DMX(Q-1);  dpv = vec_perm(dpv, infv, selectvec);
      ipv = IMX(Q-1);  ipv = vec_perm(ipv, infv, selectvec);

      
      for (q = 0; q < Q; q++)
	{
	  /* Calculate new MMX(i,q); don't store it yet, hold it in sv. */
	  sv   =             vec_add(xBv, *tsc);  tsc++;
	  sv   = vec_max(sv, vec_add(mpv, *tsc)); tsc++;
	  sv   = vec_max(sv, vec_add(ipv, *tsc)); tsc++;
	  sv   = vec_max(sv, vec_add(dpv, *tsc)); tsc++;
	  sv   = vec_add(sv, *rsc);               rsc++;
	  xEv  = vec_max(xEv, sv);
	  
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
	  dcv   = vec_add(sv, *tsc); tsc++;
	  Dmaxv = vec_max(dcv, Dmaxv);

	  /* Calculate and store I(i,q) */
	  sv     =             vec_add(mpv, *tsc);  tsc++;
	  sv     = vec_max(sv, vec_add(ipv, *tsc)); tsc++;
	  IMX(q) = vec_add(sv, *rsc);               rsc++;
	}	  

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      /* The following incantation takes the max of xEv's elements  */
      //xEv = vec_max(xEv, vec_perm(xEv, xEv, ((vector unsigned char) {0, 3, 2, 1, 0,0,0,0, 0,0,0,0, 0,0,0,0})));
      //xEv = vec_max(xEv, vec_perm(xEv, xEv, ((vector unsigned char) {1, 0, 3, 2, 0,0,0,0, 0,0,0,0, 0,0,0,0})));
      xEv = vec_max(xEv, vec_perm(xEv, xEv, ((vector unsigned char) { 4, 5, 6, 7,  8, 9,10,11, 12,13,14,15,  0, 1, 2, 3})));
      xEv = vec_max(xEv, vec_perm(xEv, xEv, ((vector unsigned char) { 8, 9,10,11, 12,13,14,15,  0, 1, 2, 3,  4, 5, 6, 7})));
      //_mm_store_ss(&xE, xEv);
      vec_ste(xEv, 0, &xE);

      xN = xN +  om->xf[p7O_N][p7O_LOOP];
      xC = ESL_MAX(xC + om->xf[p7O_C][p7O_LOOP],  xE + om->xf[p7O_E][p7O_MOVE]);
      xJ = ESL_MAX(xJ + om->xf[p7O_J][p7O_LOOP],  xE + om->xf[p7O_E][p7O_LOOP]);
      xB = ESL_MAX(xJ + om->xf[p7O_J][p7O_MOVE],  xN + om->xf[p7O_N][p7O_MOVE]);
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
      //Dmaxv = vec_max(Dmaxv, vec_perm(Dmaxv, Dmaxv, ((vector unsigned char) {0, 3, 2, 1, 0,0,0,0, 0,0,0,0, 0,0,0,0})));
      //Dmaxv = vec_max(Dmaxv, vec_perm(Dmaxv, Dmaxv, ((vector unsigned char) {1, 0, 3, 2, 0,0,0,0, 0,0,0,0, 0,0,0,0})));
      Dmaxv = vec_max(Dmaxv, vec_perm(Dmaxv, Dmaxv, ((vector unsigned char) { 4, 5, 6, 7,  8, 9,10,11, 12,13,14,15,  0, 1, 2, 3})));
      Dmaxv = vec_max(Dmaxv, vec_perm(Dmaxv, Dmaxv, ((vector unsigned char) { 8, 9,10,11, 12,13,14,15,  0, 1, 2, 3,  4, 5, 6, 7})));
      //_mm_store_ss(&Dmax, Dmaxv);
      vec_ste(Dmaxv, 0, &Dmax);
      if (Dmax + om->ddbound_f > xB) 
	{
	  /* Now we're obligated to do at least one complete DD path to be sure. */
	  /* dcv has carried through from end of q loop above */
	  //dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
	  //dcv = _mm_move_ss(dcv, infv);
          dcv = vec_perm(dcv, infv, selectvec);
	  tsc = om->tf + 7*Q;	/* set tsc to start of the DD's */
	  for (q = 0; q < Q; q++) 
	    {
	      DMX(q) = vec_max(dcv, DMX(q));	
	      dcv    = vec_add(DMX(q), *tsc); tsc++;
	    }

	  /* We may have to do up to three more passes; the check
	   * is for whether crossing a segment boundary can improve
	   * our score. 
	   */
	  do {
	    //dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
	    //dcv = _mm_move_ss(dcv, infv);
            dcv = vec_perm(dcv, infv, selectvec);
	    tsc = om->tf + 7*Q;	/* set tsc to start of the DD's */
	    for (q = 0; q < Q; q++) 
	      {
		if (! vec_any_gt(dcv, DMX(q))) break;
		DMX(q) = vec_max(dcv, DMX(q));	
		dcv    = vec_add(DMX(q), *tsc);   tsc++;
	      }	    
	  } while (q == Q);
	}
      else
	{ /* not calculating DD? then just store that last MD vector we calc'ed. */
	  //dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
	  //dcv = _mm_move_ss(dcv, infv);
          dcv = vec_perm(dcv, infv, selectvec);
	  DMX(0) = dcv;
	}

#if p7_DEBUGGING
      if (ox->debugging) omx_dump_float_row(ox, FALSE, i, 5, 2, xE, xN, xJ, xB, xC); /* logify=FALSE, <rowi>=i, width=5, precision=2*/
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T */
  *ret_sc = xC + om->xf[p7O_C][p7O_MOVE];
  return eslOK;
}
/*------------------ end, p7_ViterbiScore() ---------------------*/



/*****************************************************************
 * 10. Benchmark drivers.
 *****************************************************************/

#if defined(p7IMPL_VMX_BENCHMARK) || defined(p7IMPL_VMX_TESTDRIVE) || defined(p7IMPL_VMX_EXAMPLE)
/* Here's a couple of useful debugging functions, used in both the
 * benchmark and testdriver. (They're used in the benchmark for
 * additional manual testing purposes.)
 */

/* round_profile()
 * Round all the scores in a generic (lspace) P7_PROFILE in exactly the same
 * way that the scores in a lspace uchar P7_OPROFILE were rounded.  Then the
 * two profiles should give identical internal scores.
 *
 * Do not call this  more than once on any given <gm>!
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
}



/* simulate_msv_in_generic_profile()
 * Set a generic profile's scores so that the normal dp_generic DP 
 * algorithms will give the same score as MSVFilter():
 *   1. All t_MM scores = 0
 *   2. All other core transitions = -inf
 *   3. multihit local mode
 *   4. All t_BMk entries uniformly log 2/(M(M+1))
 *   5. tCC, tNN, tJJ scores 0; total approximated later as -3
 *   6. rounded in the same way as the 8-bit limited precision.
 * 
 * Scores in model <gm> are modified accordingly.
 *
 * Do not call this more than once on any given <gm>!
 * xref J2/79
 */
static void
simulate_msv_in_generic_profile(P7_PROFILE *gm, P7_OPROFILE *om, int L)
{
  int    k;

  esl_vec_FSet(gm->tsc, p7P_NTRANS * gm->M, -eslINFINITY);
  for (k = 1; k <  gm->M; k++) p7P_TSC(gm, k, p7P_MM) = 0.0f;
  for (k = 0; k <  gm->M; k++) p7P_TSC(gm, k, p7P_BM) = log(2.0f / ((float) gm->M * (float) (gm->M+1)));
  
  gm->xsc[p7P_N][p7P_LOOP] =  gm->xsc[p7P_J][p7P_LOOP] =  gm->xsc[p7P_C][p7P_LOOP] = 0;

  round_profile(om, gm);
}
#endif



/* There are two benchmark drivers.  The first benches DP algorithms,
 * the main optimization target; the second benches profile
 * conversion, which becomes part of the critical path in hmmpfam.
 * 
 * The -c option is useful in debugging, for comparing to known
 * (extensively tested) answers from the GViterbi() and GForward()
 * algorithms. However, watch out:
 *    - for testing ViterbiFilter(), you want to use -cx; the -x option
 *      rounds the scores in a generic profile the same way used
 *      in the optimized profile. Otherwise, you'll see the
 *      differences expected from the lack of precision in uchars.
 *      
 *    - for testing ForwardFilter(), you need to go over to 
 *      logsum.c::p7_FLogsum() and have it calculate log(exp(s1) + exp(s2),
 *      rather than using its lookup table, otherwise you'll see
 *      differences caused by lack of precision in p7_FLogsum().
 */
#ifdef p7IMPL_VMX_BENCHMARK
/* 
   gcc -o benchmark-sse -std=gnu99 -g -Wall -msse2 -I. -L. -I../easel -L../easel -Dp7IMPL_SSE_BENCHMARK impl_sse.c -lhmmer -leasel -lm 
   icc -o benchmark-sse -O3 -static -I. -L. -I../easel -L../easel -Dp7IMPL_SSE_BENCHMARK impl_sse.c -lhmmer -leasel -lm 

   ./benchmark-sse <hmmfile>       runs benchmark on ViterbiFilter() (-M for MSVFilter; -F for ForwardFilter, -S for ViterbiScore)
   ./benchmark-sse -b <hmmfile>    gets baseline time to subtract: just random seq generation
   ./benchmark-sse -c <hmmfile>    compare scores of SSE to generic impl

   ./benchmark-sse -Mx -N100 <hmmfile>     test that MSVFilter scores match Viterbi
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_vmx.h"

#define ALGOPTS "-V,-F,-S"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-b",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "baseline timing: don't run DP at all",             0 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "compare scores of generic, VMX DP        (debug)", 0 }, 
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                  0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose: show individual scores",               0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "round generic profile, make scores match (debug)", 0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  { "-M",        eslARG_NONE,   FALSE, NULL, NULL,ALGOPTS, NULL, NULL, "benchmark p7_MSVFilter()",                         0 },
  { "-V",        eslARG_NONE,"default",NULL, NULL,ALGOPTS, NULL, NULL, "benchmark p7_ViterbiFilter()",                     0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,ALGOPTS, NULL, NULL, "benchmark p7_ForwardFilter()",                     0 },
  { "-S",        eslARG_NONE,   FALSE, NULL, NULL,ALGOPTS, NULL, NULL, "benchmark p7_ViterbiScore()",                      0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for the VMX DP implementations";


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
  p7_oprofile_ReconfigLength(om, L);
  if (esl_opt_GetBoolean(go, "-x")) {
    if (esl_opt_GetBoolean(go, "-M")) simulate_msv_in_generic_profile(gm, om, L);
    else                              round_profile(om, gm);
  }
  if (esl_opt_GetBoolean(go, "-S")) pspace_to_lspace_float(om);

  ox = p7_omx_Create(gm->M);
  gx = p7_gmx_Create(gm->M, L);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      if (! esl_opt_GetBoolean(go, "-b")) {
	if      (esl_opt_GetBoolean(go, "-M")) p7_MSVFilter    (dsq, L, om, ox, &sc1);   
	else if (esl_opt_GetBoolean(go, "-F")) p7_ForwardFilter(dsq, L, om, ox, &sc1);   
	else if (esl_opt_GetBoolean(go, "-S")) p7_ViterbiScore (dsq, L, om, ox, &sc1);   
	else                                   p7_ViterbiFilter(dsq, L, om, ox, &sc1);   

	/* -c option: compare generic to fast score */
	if (esl_opt_GetBoolean(go, "-c")) 
	  {
	    if       (esl_opt_GetBoolean(go, "-F"))    p7_GForward(dsq, L, gm, gx, &sc2); 
	    else if  (esl_opt_GetBoolean(go, "-M"))    p7_GMSV    (dsq, L, gm, gx, &sc2); 
	    else                                       p7_GViterbi(dsq, L, gm, gx, &sc2); 

	    printf("%.4f %.4f\n", sc1, sc2);  
	  }

	/* -x option: compare generic to fast score in a way that should give exactly the same result */
	if (esl_opt_GetBoolean(go, "-x")) {
	  if       (esl_opt_GetBoolean(go, "-F")) p7_GForward(dsq, L, gm, gx, &sc2); 
	  else {
	    p7_GViterbi(dsq, L, gm, gx, &sc2); 
	    sc2 /= om->scale;
	    if (om->mode == p7_UNILOCAL)   sc2 -= 2.0; /* that's ~ L \log \frac{L}{L+2}, for our NN,CC,JJ */
	    else if (om->mode == p7_LOCAL) sc2 -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
	  }
	  printf("%.4f %.4f\n", sc1, sc2);  
	}
      }
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);

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
#endif /*p7IMPL_VNX_BENCHMARK*/




#ifdef p7IMPL_VMX_BENCHMARK2
/* The second benchmark driver is for timing profile conversion.
   gcc -o benchmark-sse -std=gnu99 -g -Wall -msse2 -I. -L. -I../easel -L../easel -Dp7IMPL_VMX_BENCHMARK2 impl_sse.c -lhmmer -leasel -lm 
   icc -o benchmark-sse -O3 -static -I. -L. -I../easel -L../easel -Dp7IMPL_VMX_BENCHMARK2 impl_sse.c -lhmmer -leasel -lm 

   ./benchmark-sse <hmmfile>         runs benchmark
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_vmx.h"

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
#endif /*p7IMPL_VMX_BENCHMARK2*/



/*****************************************************************
 * 11. Unit tests
 *****************************************************************/

#ifdef p7IMPL_VMX_TESTDRIVE
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

static void
make_random_profile(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L,
		    P7_HMM **ret_hmm, P7_PROFILE **ret_gm, P7_OPROFILE **ret_om)
{
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_OPROFILE    *om   = NULL;

  if (p7_hmm_Sample(r, M, abc, &hmm)                != eslOK) esl_fatal("failed to sample an HMM");
  if ((gm = p7_profile_Create(hmm->M, abc))         == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL)    != eslOK) esl_fatal("failed to config profile");
  if (p7_hmm_Validate    (hmm, NULL, 0.0001)        != eslOK) esl_fatal("whoops, HMM is bad!");
  if (p7_profile_Validate(gm,  NULL, 0.0001)        != eslOK) esl_fatal("whoops, profile is bad!");
  if ((om = p7_oprofile_Create(M, abc))             == NULL)  esl_fatal("failed to create optimized profile");
  if (p7_oprofile_Convert(gm, om)                   != eslOK) esl_fatal("failed to convert profile to optimized form");
  if (p7_oprofile_ReconfigLength(om, L)             != eslOK) esl_fatal("failed to config length of oprofile");

  *ret_hmm = hmm;
  *ret_gm  = gm;
  *ret_om  = om;
}

/* MSVFilter() unit test
 * 
 * We can check that scores are identical (within machine error) to
 * scores of generic DP with scores rounded the same way.  Do this for
 * a random model of length <M>, for <N> test sequences of length <L>.
 * 
 * We assume that we don't accidentally generate a high-scoring random
 * sequence that overflows MSVFilter()'s limited range.
 * 
 */
static void
utest_msv_filter(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *ox  = p7_omx_Create(M);
  P7_GMX      *gx  = p7_gmx_Create(M, L);
  float sc1, sc2;

  make_random_profile(r, abc, bg, M, L, &hmm, &gm, &om);
  simulate_msv_in_generic_profile(gm, om, L);

#ifdef p7_DEBUGGING
  p7_oprofile_Dump(stdout, om);              //dumps the optimized profile
  p7_omx_SetDumpMode(stdout, ox, TRUE);      //makes the fast DP algorithms dump their matrices
#endif

  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_MSVFilter(dsq, L, om, ox, &sc1);
      p7_GViterbi (dsq, L, gm, gx, &sc2);

#ifdef p7_DEBUGGING
      p7_gmx_Dump(stdout, gx);           //dumps a generic DP matrix
#endif

      sc2 = sc2 / om->scale - 3.0f;
      if (fabs(sc1-sc2) > 0.001) esl_fatal("msv filter unit test failed: scores differ (%.2f, %.2f)", sc1, sc2);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}


/* ViterbiFilter() unit test
 * 
 * We can check that scores are identical (within machine error) to
 * scores of generic DP with scores rounded the same way.  Do this for
 * a random model of length <M>, for <N> test sequences of length <L>.
 * 
 * We assume that we don't accidentally generate a high-scoring random
 * sequence that overflows ViterbiFilter()'s limited range.
 * 
 */
static void
utest_viterbi_filter(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *ox  = p7_omx_Create(M);
  P7_GMX      *gx  = p7_gmx_Create(M, L);
  float sc1, sc2;

  make_random_profile(r, abc, bg, M, L, &hmm, &gm, &om);
  round_profile(om, gm);	/* round and scale the scores in <gm> the same as in <om> */

#ifdef p7_DEBUGGING
  p7_oprofile_Dump(stdout, om);              // dumps the optimized profile
  p7_omx_SetDumpMode(stdout, ox, TRUE);      // makes the fast DP algorithms dump their matrices
#endif

  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_ViterbiFilter(dsq, L, om, ox, &sc1);
      p7_GViterbi     (dsq, L, gm, gx, &sc2);

#ifdef p7_DEBUGGING
      p7_gmx_Dump(stdout, gx);              // dumps a generic DP matrix
#endif

      sc2 /= om->scale;
      if (om->mode == p7_UNILOCAL)   sc2 -= 2.0; /* that's ~ L \log \frac{L}{L+2}, for our NN,CC,JJ */
      else if (om->mode == p7_LOCAL) sc2 -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */

      if (fabs(sc1-sc2) > 0.001) esl_fatal("viterbi filter unit test failed: scores differ (%.2f, %.2f)", sc1, sc2);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}

/* ForwardFilter() unit test
 * 
 * The generic Forward() implementation uses FLogsum(), which incurs a
 * certain amount of discretization error from its lookup table, so we
 * can't compare scores too closely. (If we had a way of replacing
 * FLogsum() with a slow but accurate log(exp+exp) version, we could.)
 */
static void
utest_forward_filter(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *ox  = p7_omx_Create(M);
  P7_GMX      *gx  = p7_gmx_Create(M, L);
  float sc1, sc2;

  make_random_profile(r, abc, bg, M, L, &hmm, &gm, &om);
  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_ForwardFilter(dsq, L, om, ox, &sc1);
      p7_GForward     (dsq, L, gm, gx, &sc2);

      if (fabs(sc1-sc2) > 1.0) esl_fatal("forward filter unit test failed: scores differ (%.2f, %.2f)", sc1, sc2);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}


/* ViterbiScore() unit test
 * 
 * We can compare these scores to GViterbi() almost exactly; the only
 * differences should be negligible roundoff errors. Must convert
 * the optimized profile to lspace, though, rather than pspace.
 */
static void
utest_viterbi_score(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *ox  = p7_omx_Create(M);
  P7_GMX      *gx  = p7_gmx_Create(M, L);
  float sc1, sc2;

  make_random_profile(r, abc, bg, M, L, &hmm, &gm, &om);
  pspace_to_lspace_float(om);
  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_ViterbiScore(dsq, L, om, ox, &sc1);
      p7_GViterbi    (dsq, L, gm, gx, &sc2);

      if (fabs(sc1-sc2) > 0.001) esl_fatal("viterbi score unit test failed: scores differ (%.2f, %.2f)", sc1, sc2);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7IMPL_VMX_TESTDRIVE*/

/*****************************************************************
 * 12. Test driver
 *****************************************************************/
#ifdef p7IMPL_VMX_TESTDRIVE
/* 
   gcc -g -Wall -msse2 -std=gnu99 -I. -L. -I../easel -L../easel -o impl_sse_utest -Dp7IMPL_SSE_TESTDRIVE impl_sse.c -lhmmer -leasel -lm
   ./impl_sse_utest
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"
#include "impl_vmx.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the VMX implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = NULL;
  ESL_ALPHABET   *abc  = NULL;
  P7_BG          *bg   = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("MSVFilter() tests, DNA\n");
  utest_msv_filter(r, abc, bg, M, L, N);   /* normal sized models */
  utest_msv_filter(r, abc, bg, 1, L, 10);  /* size 1 models       */
  utest_msv_filter(r, abc, bg, M, 1, 10);  /* size 1 sequences    */

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter() tests, DNA\n");
  utest_viterbi_filter(r, abc, bg, M, L, N);   
  utest_viterbi_filter(r, abc, bg, 1, L, 10);  
  utest_viterbi_filter(r, abc, bg, M, 1, 10);  

  if (esl_opt_GetBoolean(go, "-v")) printf("ForwardFilter() tests, DNA\n");
  utest_forward_filter(r, abc, bg, M, L, N);   
  utest_forward_filter(r, abc, bg, 1, L, 10);  
  utest_forward_filter(r, abc, bg, M, 1, 10);  

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiScore() tests, DNA\n");
  utest_viterbi_score(r, abc, bg, M, L, N);   
  utest_viterbi_score(r, abc, bg, 1, L, 10);  
  utest_viterbi_score(r, abc, bg, M, 1, 10);  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("MSVFilter() tests, protein\n");
  utest_msv_filter(r, abc, bg, M, L, N);   
  utest_msv_filter(r, abc, bg, 1, L, 10);  
  utest_msv_filter(r, abc, bg, M, 1, 10);  

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter() tests, protein\n");
  utest_viterbi_filter(r, abc, bg, M, L, N); 
  utest_viterbi_filter(r, abc, bg, 1, L, 10);
  utest_viterbi_filter(r, abc, bg, M, 1, 10);

  if (esl_opt_GetBoolean(go, "-v")) printf("ForwardFilter() tests, protein\n");
  utest_forward_filter(r, abc, bg, M, L, N);   
  utest_forward_filter(r, abc, bg, 1, L, 10);  
  utest_forward_filter(r, abc, bg, M, 1, 10);  

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiScore() tests, protein\n");
  utest_viterbi_score(r, abc, bg, M, L, N);   
  utest_viterbi_score(r, abc, bg, 1, L, 10);  
  utest_viterbi_score(r, abc, bg, M, 1, 10);  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*IMPL_VMX_TESTDRIVE*/

/*****************************************************************
 * 13. Example
 *****************************************************************/

#ifdef p7IMPL_VMX_EXAMPLE
/* A minimal example.
   Also useful for debugging on small HMMs and sequences.

   gcc -g -Wall -msse2 -std=gnu99 -I. -L. -I../easel -L../easel -o example -Dp7IMPL_SSE_EXAMPLE impl_sse.c -lhmmer -leasel -lm
   ./example <hmmfile> <seqfile>
 */ 
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "impl_vmx.h"

int 
main(int argc, char **argv)
{
  char           *hmmfile = argv[1];
  char           *seqfile = argv[2];
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
  float           sc;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  /* Read in one sequence */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
  if  (esl_sqio_Read(sqfp, sq) != eslOK) p7_Fail("Failed to read sequence");

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  /* allocate DP matrices, both a generic and an optimized one */
  ox = p7_omx_Create(gm->M);
  gx = p7_gmx_Create(gm->M, sq->n);

  /* Useful to place and compile in for debugging:  */
#ifdef p7_DEBUGGING
     p7_oprofile_Dump(stdout, om);      //dumps the optimized profile
     p7_omx_SetDumpMode(stdout,ox, TRUE);      //makes the fast DP algorithms dump their matrices
     p7_gmx_Dump(stdout, gx);           //dumps a generic DP matrix
#endif
  /**/

  simulate_msv_in_generic_profile(gm, om, sq->n);

  /* take your pick: */
  p7_MSVFilter    (sq->dsq, sq->n, om, ox, &sc);  printf("msv filter score:     %.2f nats\n", sc);
  p7_ViterbiFilter(sq->dsq, sq->n, om, ox, &sc);  printf("viterbi filter score: %.2f nats\n", sc);
  p7_ForwardFilter(sq->dsq, sq->n, om, ox, &sc);  printf("forward filter score: %.2f nats\n", sc);
  p7_GViterbi     (sq->dsq, sq->n, gm, gx, &sc);  printf("viterbi (generic):    %.2f nats\n", sc);
  p7_GForward     (sq->dsq, sq->n, gm, gx, &sc);  printf("forward (generic):    %.2f nats\n", sc);

  /* Viterbi score requires a special config of the optimized profile.
   * This isn't the final design of our API: the pspace_ call is an internal function. */
  pspace_to_lspace_float(om);
  p7_ViterbiScore (sq->dsq, sq->n, om, ox, &sc);  printf("viterbi score (VMX):  %.2f nats\n", sc);

  /* now in a real app, you'd need to convert raw nat scores to final bit
   * scores, by subtracting the null model score and rescaling.
   */

  /* cleanup */
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
  return 0;
}
#endif /*p7IMPL_VMX_EXAMPLE*/

#else /*! p7_IMPL_VMX*/
/* The remainder of the file is just bookkeeping, for what to do when
 * we aren't compiling with VMX/Altivec instructions.
 */

/*
 * Provide a successful unit test on platforms where we don't have VMX instructions.
 */
#ifdef p7IMPL_VMX_TESTDRIVE
int main(void) { return 0; }
#endif

#endif /*p7_IMPL_VMX or not*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
