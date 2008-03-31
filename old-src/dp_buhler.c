/* dp_buhler.c
 * 
 * Jeremy Buhler's optimized DP implementation for HMMER.
 * For information, please contact jbuhler@cse.wustl.edu.
 *
 * Part 1:  P7_OPROFILE profile structure API
 * Part 2:  P7_OMX DP matrix API
 * Part 3:  dynamic programming routine API 
 * 
 * Copyright (C) 2005 Washington University School of Medicine
 * Originally written by Christopher Swope and Jeremy Buhler;
 * minor modifications by SRE.
 * SVN $Id$
 */
#include "config.h"
#include "profile.h"		/* generic P7_PROFILE */
#include "dp_buhler.h"		/* optimized DP       */

/*****************************************************************
 * Part 1. The optimized profile structure (P7_OPROFILE)
 *****************************************************************/

/* Function: p7_oprofile_Create()
 * Date:     CRS, 10 June 2005 [J. Buhler's student, St. Louis]
 * 
 * Purpose:  Create a new optimized score profile. 
 * 
 * Args:     gm - the generic profile (usually <hmm->gm>)
 *
 * Returns:  a ptr to the new <P7_OPROFILE> space. 
 * 
 * Note:     Yes, you need separate Create() and Convert() functions;
 *           you may someday Create() and Read() directly into the structure.
 */
inline P7_OPROFILE *
p7_oprofile_Create(int M)
{ 
  P7_OPROFILE *om;
  unsigned int a;
  int         *miscBase;

  /* Allocation.
   */
  om       = MallocOrDie(sizeof(P7_OPROFILE));
  om->tsc  = MallocOrDie((gm->M + 1) * N_TMOVES * sizeof(int));
  om->misc = MallocOrDie(MAXCODE * sizeof(int *));
  miscBase = MallocOrDie(MAXCODE * (gm->M + 1) * N_MISC * sizeof(int));
  for (a = 0; a < MAXCODE; a++)
    om->misc[a] = miscBase + a * (gm->M + 1) * N_MISC;
  om->M    = M;
  return om;
}

/* Function:  p7_oprofile_Convert()
 * Incept:    CRS, 13 July 2005 [J. Buhler's student, St. Louis]
 *
 * Purpose:   Transfers data from a generic profile to the optimized
 *            profile.
 */
inline void
p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om)
{ 
  unsigned int k, a;
  
  /* Main model transitions
   */
  /* we only use values 0..M-1; 0 values are always -INFTY */
  for (k = 0; k < gm->M; k++)
    {
      int *tsck = om->tsc + k * N_TMOVES;
      tsck[_TMM] = gm->tsc[TMM][k];
      tsck[_TMI] = gm->tsc[TMI][k];
      tsck[_TMD] = gm->tsc[TMD][k];
      tsck[_TDM] = gm->tsc[TDM][k];
      tsck[_TDD] = gm->tsc[TDD][k];
      tsck[_TIM] = gm->tsc[TIM][k];
      tsck[_TII] = gm->tsc[TII][k];
    }
  
  /* Begin/end
   */
  /* we only use values 1..M */
  for (k = 1; k <= gm->M; k++)
    {
      int *tsck = om->tsc + k * N_TMOVES;
      tsck[_BSC] = gm->bsc[k];
      tsck[_ESC] = gm->esc[k];
    }
  
  /* Emissions (match/insert)
   */
  for (a = 0; a < MAXCODE; a++)
    {
      /* we only use values 1..M */
      for (k = 1; k <= gm->M; k++)
        {
          int *misck = om->misc[a] + k * N_MISC;
          misck[_MSC] = gmm->msc[a][k];
          misck[_ISC] = gmm->isc[a][k];
        }
    }
  
  om->xsc[_XTE][_LOOP] = gm->xsc[XTE][LOOP];  
  om->xsc[_XTE][_MOVE] = gm->xsc[XTE][MOVE];
  om->xsc[_XTN][_LOOP] = gm->xsc[XTN][LOOP];
  om->xsc[_XTN][_MOVE] = gm->xsc[XTN][MOVE];
  om->xsc[_XTC][_LOOP] = gm->xsc[XTC][LOOP];
  om->xsc[_XTC][_MOVE] = gm->xsc[XTC][MOVE];
  om->xsc[_XTJ][_LOOP] = gm->xsc[XTJ][LOOP];
  om->xsc[_XTJ][_MOVE] = gm->xsc[XTJ][MOVE];
  return;
}


/* Function: p7_oprofile_Destroy()
 * Date:     CRS, 10 June 2005 [J. Buhler's student, St. Louis]
 *
 * Purpose:  Frees a <P7_OPROFILE>.
 *
 * Returns:  (void)
 */
inline void 
p7_oprofile_Destroy(P7_OPROFILE *om)
{
  free(om->tsc);
  free(om->misc[0]);
  free(om->misc);
  free(om);
}



/*****************************************************************
 * Part 2. The P7_OMX DP matrix API.
 *****************************************************************/





/*****************************************************************
 * Part 3. The dynamic programming routine API.
 *****************************************************************/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
