/* H4_COUNTS: profile HMM in double-precision counts-collection form.
 * 
 * Contents:
 *    1. H4_COUNTS structure
 *    2. Debugging and development tools
 */
#include <h4_config.h>

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_matrixops.h"

#include "h4_counts.h"

/*****************************************************************
 * 1. H4_COUNTS structure
 *****************************************************************/ 

/* Function:  h4_counts_Create()
 * Synopsis:  Allocate a new count-collection profile.
 * Incept:    SRE, Mon 20 May 2019
 *
 * Purpose:   Allocates a new count-collection profile of 
 *            length <M> consensus positions, for alphabet <abc>. 
 *            
 *            All transition counts <t> and emission counts <e> are
 *            initialized to zero; <M> and <abc> fields are set.
 */
H4_COUNTS *
h4_counts_Create(const ESL_ALPHABET *abc, int M)
{
  H4_COUNTS *ctm = NULL;
  int        status;

  ESL_ALLOC(ctm, sizeof(H4_COUNTS));
  ctm->t = ctm->e = NULL;

  if (( ctm->t = esl_mat_DCreate( M+1, h4_NT))  == NULL) goto ERROR;
  if (( ctm->e = esl_mat_DCreate( M+1, abc->K)) == NULL) goto ERROR;

  esl_mat_DSet(ctm->t, M+1, h4_NT,  0.);
  esl_mat_DSet(ctm->e, M+1, abc->K, 0.);

  ctm->M   = M;
  ctm->abc = abc;
  return ctm;

 ERROR:
  h4_counts_Destroy(ctm);
  return NULL;
}

/* Function:  h4_counts_Destroy()
 * Synopsis:  Free a H4_COUNTS structure
 * Incept:    SRE, Mon 20 May 2019
 */
void
h4_counts_Destroy(H4_COUNTS *ctm)
{
  if (ctm)
    {
      esl_mat_DDestroy(ctm->t);
      esl_mat_DDestroy(ctm->e);
      free(ctm);
    }
}
  


/*****************************************************************
 * 2. Debugging and development tools
 *****************************************************************/

/* Function:  h4_counts_Dump()
 * Synopsis:  Dump contents of H4_COUNTS for inspection
 * Incept:    SRE, Mon 20 May 2019
 */
int
h4_counts_Dump(FILE *fp, H4_COUNTS *ctm)
{
  int k,a,z;

  fprintf(fp, "Emission counts:\n");
  fprintf(fp, "     ");
  for (a = 0; a < ctm->abc->K; a++)
    fprintf(fp, "         %c%c", ctm->abc->sym[a], a == ctm->abc->K-1 ? '\n':' ');
  for (k = 1; k <= ctm->M; k++)
    {
      fprintf(fp, "%4d ", k);
      for (a = 0; a < ctm->abc->K; a++)
	fprintf(fp, "%10.4f%c", ctm->e[k][a], a == ctm->abc->K-1 ? '\n':' ');
    }

  fprintf(fp, "Transition counts:\n");
  fprintf(fp, "     ");
  fprintf(fp, "%10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
	  "TMM", "TMI", "TMD", "TIM", "TII", "TID", "TDM", "TDI", "TDD");
  for (k = 0; k < ctm->M; k++)  // include t[0] because GM, GD (in MM, MD) are data-dependent
    {                           // exclude t[M] which has no data dependent transitions.
      fprintf(fp, "%4d ", k);
      for (z = 0; z < h4_NT; z++)
	fprintf(fp, "%10.4f%c", ctm->t[k][z], z == h4_NT-1 ? '\n':' ');
    }
  return eslOK;
}
