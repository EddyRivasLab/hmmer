/* Using P7_SPARSEMX for AEC (anchor/envelope constrained) alignments.
 * 
 * p7_spaecmx does not provide an independent data structure. Rather,
 * we use the sparse DP matrix P7_SPARSEMX in three different ways:
 *     normal sparse                    : p7_sparsemx.[ch]
 *     ASC, anchor set constrained      : p7_spascmx.[ch]
 *     AEC, anchor/envelope constrained : p7_spaecmx.[ch]
 * and this module provides stuff for the AEC part.
 */
#include "p7_config.h"

#include <stdio.h>

#include "base/p7_envelopes.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/p7_spaecmx.h"



/* Function:  p7_spaecmx_Dump()
 * Synopsis:  Dump a sparse AEC DP matrix.
 *
 * Purpose:   Dump sparse AEC dynamic programming matrix <aec>,
 *            which was constrained by the envelopes in <env>,
 *            to open output stream <fp>.       
 *            
 *            The AEC matrix only uses L cells, not G. It works with
 *            posteriors marginalized over L/G choice.  Local/glocal
 *            decision is made already, in envelope determination.
 *            Therefore we only dump three rows per supercell, labeled
 *            M/I/D.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_spaecmx_Dump(FILE *fp, const P7_SPARSEMX *aec, const P7_ENVELOPES *env)
{
  const P7_SPARSEMASK *sm  = aec->sm;
  const float         *dpc = aec->dp;
  const float         *last_dpc;
  int   width     = 9;
  int   precision = 4;
  int   d,i,k,z,s;
  int   cells[3]  = { p7S_ML, p7S_IL, p7S_DL };
  char  labels[4] = "MID";

  ESL_DASSERT1(( aec->type == p7S_AEC_ALIGN ));
  ESL_DASSERT1(( sm->M == env->M ));
  ESL_DASSERT1(( sm->L == env->L ));


  fprintf(fp, "Sparse AEC matrix dump. L=%d, M=%d\n", sm->L, sm->M);
  fprintf(fp, "       ");
  for (k = 1; k <= sm->M; k++) fprintf(fp, "%*d ", width, k);
  fprintf(fp, "\n");
  
  fprintf(fp, "       ");
  for (k = 1; k <= sm->M; k++) fprintf(fp, "%*.*s ", width, width, "---------");
  fprintf(fp, "\n");

  for (d = 1; d <= env->D; d++)
    {
      /* UP sector:   ia(d)..i0(d)-1; 1..k0(d)-1 */
      for (i = env->arr[d].ia; i < env->arr[d].i0; i++)
	{
	  for (s = 0; s < 3; s++) 
	    {
	      fprintf(fp, "%3d %c ", i, labels[s]);
	      for (z = 0, k = 1; k < env->arr[d].k0; k++) {
		while (z < sm->n[i] && sm->k[i][z]  < k) z++;
		if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(fp, "%*.*f ", width, precision, *(dpc + z*p7S_NSCELLS + cells[s]));
		else                                     fprintf(fp, "%*s ",   width, "......");
	      }
	      for (; k <= sm->M; k++) {
		while (z < sm->n[i] && sm->k[i][z]  < k) z++;
		if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(fp, "%*s ", width, "...xx.");
		else                                     fprintf(fp, "%*s ", width, "......");
	      }
	      fprintf(fp, "\n");
	    }
	  fprintf(fp, "\n");

	  /* Now advance dpc */
	  for (z = 0; z < sm->n[i] && sm->k[i][z] < env->arr[d].k0; z++) dpc += p7S_NSCELLS;
	}


      /* DOWN sector: i0(d)..ib(d);   k0(d)..M   */
      for (i = env->arr[d].i0; i <= env->arr[d].ib; i++)
	{
	  last_dpc = dpc;
	  for (s = 0; s < 3; s++)
	    {
	      dpc = last_dpc;
	      fprintf(fp, "%3d %c ", i, labels[s]);
	      for (z = 0, k = 1; k < env->arr[d].k0; k++) {
		while (z < sm->n[i] && sm->k[i][z]  < k) z++;
		if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(fp, "%*s ", width, "...xx.");
		else                                     fprintf(fp, "%*s ", width, "......");
	      }
	      for (; k <= sm->M; k++) {
		while (z < sm->n[i] && sm->k[i][z]  < k)   z++;
		if    (z < sm->n[i] && sm->k[i][z] == k) { fprintf(fp, "%*.*f ", width, precision, *(dpc + cells[s])); dpc += p7S_NSCELLS; }
		else                                       fprintf(fp, "%*s ", width, "......");
	      }
	      fprintf(fp, "\n");
	    }
	  fprintf(fp, "\n");
	}

      fprintf(fp, "...\n\n");
    }
  return eslOK;
}


/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
