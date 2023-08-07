/* Using P7_SPARSEMX for AEC (anchor/envelope constrained) alignments.
 * 
 * p7_spaecmx does not provide an independent data structure. Rather,
 * we use the sparse DP matrix P7_SPARSEMX in three different ways:
 *     normal sparse                    : p7_sparsemx.[ch]
 *     ASC, anchor set constrained      : p7_spascmx.[ch]
 *     AEC, anchor/envelope constrained : p7_spaecmx.[ch]
 * and this module provides stuff for the AEC part.
 */
#include <p7_config.h>

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
 * 2. AEC matrix validation.
 *****************************************************************/

/* Only L (local) cells in main supercells are used in an AEC
 * matrix. G (glocal) cells and specials are unused, untouched, may be
 * arbitrary values.
 * 
 * For ML/IL/DL cells, AEC score values are either -inf or >= 0.0.
 */
static int
aec_supercell(const float *dpc, int M_used, int I_used, int D_used)
{
  int aectype[3] = { p7S_ML, p7S_IL, p7S_DL };
  int s;
  for (s = 0; s < 3; s++)
    if      (isfinite(dpc[aectype[s]]) && dpc[aectype[s]] >= 0.0) continue;
    else if (dpc[aectype[s]] == -eslINFINITY)                     continue;
    else  ESL_FAIL(eslFAIL, NULL, NULL);

  if (! M_used && (dpc[p7S_ML] != -eslINFINITY)) ESL_FAIL(eslFAIL, NULL, NULL);
  if (! I_used && (dpc[p7S_IL] != -eslINFINITY)) ESL_FAIL(eslFAIL, NULL, NULL);
  if (! D_used && (dpc[p7S_DL] != -eslINFINITY)) ESL_FAIL(eslFAIL, NULL, NULL);
  return eslOK;
}

/* Function:  p7_spaecmx_Validate()
 * Synopsis:  Validate a sparse AEC DP matrix.
 *
 * Purpose:   Validate sparse AEC DP alignment matrix <aec>, which was 
 *            constrained by envelopes <env>. Return <eslOK> if it's ok.
 *            Else if <aec> fails validation, return <eslFAIL>; if
 *            <errbuf> is non-NULL, leave an error message in it.
 */
int
p7_spaecmx_Validate(const P7_SPARSEMX *aec, const P7_ENVELOPES *env, char *errbuf)
{
  char                 msg[] = "sparse AEC matrix validation failed";
  const P7_SPARSEMASK *sm    = aec->sm;
  const float         *dpc   = aec->dp;
  int d,z,i;
  int ia,i0,k0,ib;


  if ( aec->type != p7S_AEC_ALIGN ) ESL_FAIL(eslFAIL, errbuf, msg);

  for (d = 1; d <= env->D; d++)
    {
      ia      =  env->arr[d].ia;
      i0      =  env->arr[d].i0;
      k0      =  env->arr[d].k0;
      ib      =  env->arr[d].ib;

      /* First row ia of UP sector, if there's an UP sector */
      if (ia < i0)
	{
	  for (z = 0; z < sm->n[ia] && sm->k[ia][z] < k0; z++, dpc += p7S_NSCELLS)
	    {
	      if (sm->k[ia][z] == 1) { if (aec_supercell(dpc, 1, 0, 0) != eslOK) ESL_FAIL(eslFAIL, errbuf, msg); }
	      else                   { if (aec_supercell(dpc, 1, 0, 1) != eslOK) ESL_FAIL(eslFAIL, errbuf, msg); }
	    }
	}

      /* Remaining rows ia+1..i0-1 of UP sector, if there's an UP sector */
      for (i = ia+1; i < i0; i++)
	{
	  for (z = 0; z < sm->n[i] && sm->k[i][z] < k0; z++, dpc += p7S_NSCELLS)
	    {
	      if (sm->k[i][z] == 1) { if (aec_supercell(dpc, 1, 1, 0) != eslOK) ESL_FAIL(eslFAIL, errbuf, msg); }
	      else                  { if (aec_supercell(dpc, 1, 1, 1) != eslOK) ESL_FAIL(eslFAIL, errbuf, msg); }
	    }
	}
	  
      /* Anchor row i0, top of DOWN sector */
      z = 0; while (z < sm->n[i0] && sm->k[i0][z] < k0) z++;
      for (;        z < sm->n[i0];                      z++, dpc += p7S_NSCELLS)
	{
	  if   (sm->k[i0][z] == k0) { if (aec_supercell(dpc, 1, 0, 0) != eslOK) ESL_FAIL(eslFAIL, errbuf, msg); }
	  else                      { if (aec_supercell(dpc, 0, 0, 1) != eslOK) ESL_FAIL(eslFAIL, errbuf, msg); }
	}

      /* Remaining rows i0+1..ib of DOWN sector, if any */
      for (i = i0+1; i <= ib; i++)
	{
	  z = 0; while (z < sm->n[i] && sm->k[i][z] < k0) z++;
	  for (;        z < sm->n[i];                     z++, dpc += p7S_NSCELLS)
	    {
	      if      (sm->k[i][z] == k0)     { if (aec_supercell(dpc, 0, 1, 0) != eslOK) ESL_FAIL(eslFAIL, errbuf, msg); }  // k0   == M impossible for i>i0
	      else if (sm->k[i][z] == k0+1)   { if (aec_supercell(dpc, 1, 1, 0) != eslOK) ESL_FAIL(eslFAIL, errbuf, msg); }  // k0+1 == M possible for i0+1 only; could check I=-inf but we don't
	      else if (sm->k[i][z] == env->M) { if (aec_supercell(dpc, 1, 0, 1) != eslOK) ESL_FAIL(eslFAIL, errbuf, msg); }  
	      else                            { if (aec_supercell(dpc, 1, 1, 1) != eslOK) ESL_FAIL(eslFAIL, errbuf, msg); }  
	    }
	}
    }

  return eslOK;
}
