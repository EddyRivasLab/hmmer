/* P7_GMX implementation: a generic dynamic programming matrix
 *
 * Contents:
 *   1. The <P7_GMX> object
 *   2. Debugging aids
 *   3. Unit tests
 *   4. Test driver
 *   5. Copyright and license information
 * 
 * SRE, Tue Jan 30 11:14:11 2007 [Einstein's, in St. Louis]
 * SVN $Id$
 */

#include "p7_config.h"
#include "hmmer.h"

/*****************************************************************
 *= 1. The <P7_GMX> object.
 *****************************************************************/

/* Function:  p7_gmx_Create()
 * Incept:    SRE, Tue Jan 30 11:20:33 2007 [Einstein's, in St. Louis]
 *
 * Purpose:   Allocate a reusable, resizeable <P7_GMX> for models up to
 *            size <allocM> and sequences up to length <allocL>.
 *            
 *            We've set this up so it should be easy to allocate
 *            aligned memory, though we're not doing this yet.
 *
 * Returns:   a pointer to the new <P7_GMX>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_GMX *
p7_gmx_Create(int allocM, int allocL)
{
  int     status;
  P7_GMX *gx = NULL;
  int     i;

  /* level 1: the structure itself */
  ESL_ALLOC(gx, sizeof(P7_GMX));
  gx->dp     = NULL;
  gx->xmx    = NULL;
  gx->dp_mem = NULL;

  /* level 2: row pointers, 0.1..L; and dp cell memory  */
  ESL_ALLOC(gx->dp,      sizeof(float *) * (allocL+1));
  ESL_ALLOC(gx->xmx,     sizeof(float)   * (allocL+1) * p7G_NXCELLS);
  ESL_ALLOC(gx->dp_mem,  sizeof(float)   * (allocL+1) * (allocM+1) * p7G_NSCELLS);

  /* Set the row pointers. */
  for (i = 0; i <= allocL; i++) 
    gx->dp[i] = gx->dp_mem + i * (allocM+1) * p7G_NSCELLS;

  /* Initialize memory that's allocated but unused, only to keep
   * valgrind and friends happy.
   */
  for (i = 0; i <= allocL; i++) 
    {
      gx->dp[i][0      * p7G_NSCELLS + p7G_M] = -eslINFINITY; /* M_0 */
      gx->dp[i][0      * p7G_NSCELLS + p7G_I] = -eslINFINITY; /* I_0 */      
      gx->dp[i][0      * p7G_NSCELLS + p7G_D] = -eslINFINITY; /* D_0 */
      gx->dp[i][1      * p7G_NSCELLS + p7G_D] = -eslINFINITY; /* D_1 */
      gx->dp[i][allocM * p7G_NSCELLS + p7G_I] = -eslINFINITY; /* I_M */
    }

  gx->M      = 0;
  gx->L      = 0;
  gx->allocW = allocM+1;
  gx->allocR = allocL+1;
  gx->validR = allocL+1;
  gx->ncells = (uint64_t) (allocM+1)* (uint64_t) (allocL+1);
  return gx;

 ERROR:
  if (gx != NULL) p7_gmx_Destroy(gx);
  return NULL;
}

/* Function:  p7_gmx_GrowTo()
 * Synopsis:  Assure that DP matrix is big enough.
 * Incept:    SRE, Tue Jan 30 11:31:23 2007 [Olin Library, St. Louis]
 *
 * Purpose:   Assures that a DP matrix <gx> is allocated
 *            for a model of size up to <M> and a sequence of
 *            length up to <L>; reallocates if necessary.
 *            
 *            This function does not respect the configured
 *            <RAMLIMIT>; it will allocate what it's told to
 *            allocate. 
 *
 * Returns:   <eslOK> on success, and <gx> may be reallocated upon
 *            return; any data that may have been in <gx> must be 
 *            assumed to be invalidated.
 *
 * Throws:    <eslEMEM> on allocation failure, and any data that may
 *            have been in <gx> must be assumed to be invalidated.
 */
int
p7_gmx_GrowTo(P7_GMX *gx, int M, int L)
{
  int      status;
  void    *p;
  int      i;
  uint64_t ncells;
  int      do_reset = FALSE;

  if (M < gx->allocW && L < gx->validR) return eslOK;
  
  /* must we realloc the 2D matrices? (or can we get away with just
   * jiggering the row pointers, if we are growing in one dimension
   * while shrinking in another?)
   */
  ncells = (uint64_t) (M+1) * (uint64_t) (L+1);
  if (ncells > gx->ncells) 
    {
      ESL_RALLOC(gx->dp_mem, p, sizeof(float) * ncells * p7G_NSCELLS);
      gx->ncells = ncells;
      do_reset   = TRUE;
    }

  /* must we reallocate the row pointers? */
  if (L >= gx->allocR)
    {
      ESL_RALLOC(gx->xmx, p, sizeof(float)   * (L+1) * p7G_NXCELLS);
      ESL_RALLOC(gx->dp,  p, sizeof(float *) * (L+1));
      gx->allocR = L+1;
      gx->allocW = M+1;
      do_reset   = TRUE;
    }

  /* must we widen the rows? */
  if (M >= gx->allocW)
    {
      gx->allocW = M+1;
      do_reset   = TRUE;
    }

  /* must we set some more valid row pointers? */
  if (L >= gx->validR)
    do_reset   = TRUE;

  /* reset all the row pointers.*/
  if (do_reset)
    {
      gx->validR = ESL_MIN(gx->ncells / gx->allocW, gx->allocR);
      for (i = 0; i < gx->validR; i++) 
	gx->dp[i] = gx->dp_mem + i * (gx->allocW) * p7G_NSCELLS;
    }

  gx->M      = 0;
  gx->L      = 0;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_gmx_Destroy()
 * Synopsis:  Frees a DP matrix.
 * Incept:    SRE, Tue Jan 30 11:17:36 2007 [Einstein's, in St. Louis]
 *
 * Purpose:   Frees a <P7_GMX>.
 *
 * Returns:   (void)
 */
void
p7_gmx_Destroy(P7_GMX *gx)
{
  if (gx == NULL) return;

  if (gx->dp      != NULL)  free(gx->dp);
  if (gx->xmx     != NULL)  free(gx->xmx);
  if (gx->dp_mem  != NULL)  free(gx->dp_mem);
  free(gx);
  return;
}

/*****************************************************************
 * 2. Debugging aids
 *****************************************************************/

/* Function:  p7_gmx_Dump()
 * Synopsis:  Dump a DP matrix to a stream, for diagnostics.
 * Incept:    SRE, Fri Jul 13 09:56:04 2007 [Janelia]
 *
 * Purpose:   Dump matrix <gx> to stream <fp> for diagnostics.
 */
int
p7_gmx_Dump(FILE *ofp, P7_GMX *gx)
{
  return p7_gmx_DumpWindow(ofp, gx, 0, gx->L, 0, gx->M, TRUE);
}


/* Function:  p7_gmx_DumpWindow()
 * Synopsis:  Dump a window of a DP matrix to a stream for diagnostics.
 * Incept:    SRE, Mon Apr 14 08:45:28 2008 [Janelia]
 *
 * Purpose:   Dump a window of matrix <gx> to stream <fp> for diagnostics,
 *            from row <istart> to <iend>, from column <kstart> to <kend>.
 *            
 *            If <show_specials> is <TRUE>, scores for the special
 *            <ENJBC> states are also displayed.
 *
 *            Asking for <0..L,0..M> with <show_specials=TRUE> is the
 *            same as <p7_gmx_Dump()>.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_gmx_DumpWindow(FILE *ofp, P7_GMX *gx, int istart, int iend, int kstart, int kend, int show_specials)
{
  int i, k, x;
  int width     = 9;
  int precision = 4;

  /* Header */
  fprintf(ofp, "     ");
  for (k = kstart; k <= kend;  k++) fprintf(ofp, "%*d ", width, k);
  if (show_specials)                fprintf(ofp, "%*s %*s %*s %*s %*s\n", width, "E", width, "N", width, "J", width, "B", width, "C");
  fprintf(ofp, "      ");
  for (k = kstart; k <= kend; k++) fprintf(ofp, "%*.*s ", width, width, "----------");
  if (show_specials)               fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");
  
  /* DP matrix data */
  for (i = istart; i <= iend; i++)
    {
      fprintf(ofp, "%3d M ", i);
      for (k = kstart; k <= kend;        k++) fprintf(ofp, "%*.*f ", width, precision, gx->dp[i][k * p7G_NSCELLS + p7G_M]);
      if (show_specials) 
	for (x = 0;    x <  p7G_NXCELLS; x++) fprintf(ofp, "%*.*f ", width, precision, gx->xmx[  i * p7G_NXCELLS + x]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d I ", i);
      for (k = kstart; k <= kend;        k++) fprintf(ofp, "%*.*f ", width, precision, gx->dp[i][k * p7G_NSCELLS + p7G_I]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d D ", i);
      for (k = kstart; k <= kend;        k++) fprintf(ofp, "%*.*f ", width, precision, gx->dp[i][k * p7G_NSCELLS + p7G_D]);
      fprintf(ofp, "\n\n");
    }
  return eslOK;
}


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7GMX_TESTDRIVE

static void
gmx_testpattern(P7_GMX *gx)
{
  int i,k,s,n, n2;

  /* Write a test pattern, via the dp[i] pointers */
  n = 0;
  for (i = 0; i < gx->validR; i++)
    for (k = 0; k < gx->allocW; k++)
      for (s = 0; s < p7G_NSCELLS; s++)
	gx->dp[i][k*p7G_NSCELLS+s] = n++;

  /* Read it back, via the dp[i] pointers */
  n = 0;
  for (i = 0; i < gx->validR; i++)
    for (k = 0; k < gx->allocW; k++)
      for (s = 0; s < p7G_NSCELLS; s++)
	{
	  if (gx->dp[i][k*p7G_NSCELLS+s] != n) esl_fatal("gmx unit test failed: test pattern corrupted");
	  n++;
	}
  
  /* Reading it back via the dp_mem vector itself ought to be the same */
  if (gx->allocR == gx->validR && gx->ncells == gx->validR*gx->allocW)
    {
      n2 = 0;
      for (i = 0; i < gx->ncells; i++)
	for (s = 0; s < p7G_NSCELLS; s++)
	  {
	    if (gx->dp_mem[i*p7G_NSCELLS+s] != n2) esl_fatal("gmx unit test failed: test pattern corrupted (2nd test)");
	    n2++;
	  }
      /* and the number of cells ought to match too */
      if (n != n2) esl_fatal("gmx unit test failed: unexpected # of cells");
    }
}


static void
utest_GrowTo(void)
{
  P7_GMX *gx = p7_gmx_Create(20, 20);
  gmx_testpattern(gx);

  p7_gmx_GrowTo(gx,  40,  20);  gmx_testpattern(gx);	/* grow in M, not L */
  p7_gmx_GrowTo(gx,  40,  40);  gmx_testpattern(gx);	/* grow in L, not M */
  p7_gmx_GrowTo(gx,  80,  10);  gmx_testpattern(gx);	/* grow in M, but with enough ncells */
  p7_gmx_GrowTo(gx,  10,  80);  gmx_testpattern(gx);	/* grow in L, but with enough ncells */
  p7_gmx_GrowTo(gx, 160, 160);  gmx_testpattern(gx);	/* grow in both L and M */
}
#endif /*p7GMX_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/


/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7GMX_TESTDRIVE
/*
  gcc -o gmx_utest -msse2 -g -Wall -I. -L. -I../easel -L../easel -Dp7GMX_TESTDRIVE p7_gmx.c -leasel -lm
  ./gmx_utest
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",               0},
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for p7_gmx.c";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);

  utest_GrowTo();

  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7GMX_TESTDRIVE*/
/*------------------ end, test driver ---------------------------*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
