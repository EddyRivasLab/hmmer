/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* emulation.c
 * SRE, Wed Jan 21 07:50:01 1998
 * 
 * Interfaces between HMMER and other software packages.
 * 
 * RCS $Id$
 */

#include <stdio.h>
#include <string.h>

#include "squid.h"
#include "config.h"
#include "structs.h"
#include "funcs.h"
#include "version.h"


/* Function: WriteProfile()
 * Date:     SRE, Wed Jan 21 07:58:09 1998 [St. Louis]
 *
 * Purpose:  Given an HMM, write a GCG profile .prf file as
 *           output. Based on examination of Michael Gribskov's Fortran
 *           source in GCG 9.0; on reverse engineering
 *           by examination of GCG 9.0 output from "profilemake"
 *           and how the .prf file is used by "profilesearch";
 *           and on the GCG 9.0 documentation.
 *           
 *           The profile for symbol i at model position k is:
 *
 *           PROF(i,k)   = log P(i,k) + TMM
 *           m_open(k)   = (TMI + TIM - TMM) - TII     
 *           m_extend(k) = TII
 *           
 *           Note that GCG affine gaps are gop + n * gex;
 *           HMMER affine gaps are gop + (n-1) * gex, thus an
 *           extra TII gets taken away from gop, since GCG will charge it.
 *           
 *           This is a simple mapping of HMM to
 *           profile. It ignores the delete transitions entirely,
 *           fitting only the match and inserts. Note the somewhat
 *           tricky method by which the M->M transition score is
 *           smuggled into the profile. 
 *           
 *           To get a close approximation of hmmsw scores, call
 *           profilesearch as
 *                profilesearch -noave -nonor -gap 10 -len 1
 *
 * Args:     fp   - open FILE to write to (or stdout, possibly)
 *           hmm  - the HMM to write   
 *
 * Returns:  (void)
 */
void
WriteProfile(FILE *fp, struct plan7_s *hmm)
{
  int k;			/* position in model      */
  int x;			/* symbol index           */
  int sc;			/* a score to print       */
  float nx;			/* expected # of symbol x */
  
  Plan7Logoddsify(hmm);

  /* GCG can't deal with long profiles. Their limit is 1000
   * positions.
   */
  if (hmm->M > 1000)
    Die("Can't convert HMM %s; M=%d, and GCG imposes a max limit of 1000 consensus positions", 
	hmm->name, hmm->M);

  /* Header information.
   * GCG 9.0 will look for sequence type and length of model.
   * Other than this, nothing is parsed until we get to the 
   * Cons line that has a ".." on it.
   * Lines that begin with "!" are comments.
   */
  if (Alphabet_type == hmmAMINO)        fprintf(fp, "!!AA_PROFILE 1.0\n");
  else if (Alphabet_type == hmmNUCLEIC) fprintf(fp, "!!NA_PROFILE 1.0\n");
  else Die("No support for profiles with non-biological alphabets");

  if (Alphabet_type == hmmAMINO)        fprintf(fp, "(Peptide) ");
  else if (Alphabet_type == hmmNUCLEIC) fprintf(fp, "(Nucleotide) ");
  fprintf(fp, "HMMCONVERT v%s of: %s  Length: %d\n",
	  RELEASE, hmm->name, hmm->M);
  
  /* Insert some HMMER-specific commentary
   */
  fprintf(fp, "   Profile converted from a profile HMM using HMMER v%s emulation.\n", RELEASE);
  fprintf(fp, "   Use -nonor -noave -gap 10.0 -len 1.0 with profilesearch and friends\n");
  fprintf(fp, "      to get the closest approximation to HMMER bit scores.\n");
  fprintf(fp, "   WARNING: There is a loss of information in this conversion.\n");
  fprintf(fp, "      Neither the scores nor even the rank order of hits will be precisely\n");
  fprintf(fp, "      preserved in a comparison of HMMER hmmsearch to GCG profilesearch.\n\n");


  /* Do the CONS line, which gives the valid IUPAC symbols and their order
   */
  fprintf(fp, "Cons");
  for (x = 0; x < Alphabet_iupac; x++)
    fprintf(fp, "    %c ", Alphabet[x]);
  fprintf(fp, "  Gap   Len ..\n");
 
  /* Now, the profile; for each position in the HMM, write a line of profile.
   */
  for (k = 1; k <= hmm->M; k++)
    {
				/* GCG adds some indexing as comments */
      if ((k-1)%10 == 0 && k > 10)
	fprintf(fp, "! %d\n", k);

				/* find consensus residue by max prob */
      x = FMax(hmm->mat[k], Alphabet_size);
      fprintf(fp, " %c  ", Alphabet[x]);
				/* generate emission score profile;
				 * Profiles are scaled by a factor of 100 
				 */
      for (x = 0; x < Alphabet_iupac; x++)
	{
	  sc = hmm->msc[x][k];
	  if (k < hmm->M) sc += hmm->tsc[k][TMM];
	  sc = sc * 100 / INTSCALE;
	  fprintf(fp, "%5d ", sc);
	}
				/* Generate gap open, gap extend penalties;
				   note profilesearch defaults of 10, 1,
				   and that GCG profile values are percentages
				   of these base penalties, 0..100.*/
				/* gap open */
      if (k < hmm->M)
	{
	  sc = -1 * (hmm->tsc[k][TMI] + hmm->tsc[k][TIM] - hmm->tsc[k][TMM] - hmm->tsc[k][TII]);
	  sc = sc * 100 / (10.0 * INTSCALE);
	}
      else sc = 100;
      fprintf(fp, "%5d ", sc);
				/* gap extend */
      if (k < hmm->M)
	{
	  sc = -1 * hmm->tsc[k][TII];
	  sc = sc * 100 / (1.0 * INTSCALE);
	}
      else sc = 100;
      fprintf(fp, "%5d\n", sc);
    }

  /* The final line of the profile is a count of the observed
   * residues in the training sequences. This information is not
   * available in an HMM, and I'm not sure that GCG ever uses it.
   * Approximate it by calculating a /very/ rough expectation.
   */
  fprintf(fp, " *  ");
  for (x = 0; x < Alphabet_size; x++)
    {
      nx = 0.0;
      for (k = 1; k <= hmm->M; k++)
	nx += hmm->mat[k][x];
      nx *= hmm->nseq;
      fprintf(fp, "%5d ", (int) nx);
    }
  for (; x < Alphabet_iupac; x++)
      fprintf(fp, "%5d ", 0);
  fprintf(fp, "\n");
  return;
}
      
