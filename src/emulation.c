/************************************************************
 * @LICENSE@
 ************************************************************/

/* emulation.c
 * SRE, Wed Jan 21 07:50:01 1998
 * 
 * Interfaces between HMMER and other software packages.
 * 
 * CVS $Id$
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <string.h>

#include "squid.h"
#include "structs.h"
#include "funcs.h"


/* Function: WriteProfile()
 * Date:     SRE, Wed Jan 21 07:58:09 1998 [St. Louis]
 *
 * Purpose:  Given an HMM, write a GCG profile .prf file as
 *           output. Based on examination of Michael Gribskov's Fortran
 *           source in GCG 9.1; on reverse engineering
 *           by examination of GCG 9.1 output from "profilemake"
 *           and how the .prf file is used by "profilesearch";
 *           and on the GCG 9.0 documentation.
 *           
 *           See notes 28 Jan 98 for detail; in brief, the conversion goes like:
 *           
 *           PROF(i,k) = match score         =  msc(i,k) + TMM(k-1)
 *           
 *           GAP(k)    = cost per insertion  =  TMI(k-1) + TIM(k-1) - TMM(k-1) - TII(k-1)
 *           LEN(k)    = cost per inserted x =  TII(k-1)
 *           
 *           QGAP(k)   = cost per deletion   =  TDM(k-1) + TMD(unknown) - TMM(k-1) - TDD(k-1)
 *           QLEN(k)   = cost per deleted k  =  TDD(k-1)
 *           
 *           Note that GCG affine gaps are GAP + n * LEN;
 *           HMMER affine gaps count (n-1) * gap-extend, thus an
 *           extra TII gets taken away from GAP (and TDD from QGAP),
 *           since GCG will charge it.
 *           
 *           Also note how the TMM transitions, which have no equivalent
 *           in a profile, get smuggled in OK.
 *           
 *           Also note that GCG charges gaps using the profile position
 *           /after/ the gap, not preceding the gap as HMMER does.
 *           
 *           Also note the TMD(unknown) in the QGAP calculation. HMMER
 *           distinguishes between gap-open and gap-close, but GCG does not,
 *           so there is a fundamental incompatibility here. Here
 *           we use an upper (best-scoring, minimum-cost) bound. 
 *           
 *           And finally note that GCG's implementation forces GAP=QGAP and
 *           LEN=QLEN. Here, we upper bound again. Compugen's implementation
 *           allows an "extended profile" format which distinguishes between
 *           the two. 
 *           
 *           The upper bound approach to these scores means that a 
 *           score given by an emulated profile is an upper bound: the HMMER
 *           score (for a single Smith/Waterman style local alignment)
 *           cannot be better than this. This is intentional, so that
 *           the Compugen BIC can be used for rapid prefiltering of
 *           the database.
 *           
 *           To get a close approximation of hmmsw scores, call
 *           profilesearch as
 *                profilesearch -noave -nonor -gap 10 -len 1
 *           On the Compugen BIC, using extended profiles, you want:
 *                om -model=xsw.model -gapop=10 -gapext=1 -qgapop=10 -qgapext=1 -noave -nonor
 *
 * Args:     fp      - open FILE to write to (or stdout, possibly)
 *           hmm     - the HMM to write   
 *           do_xsw  - TRUE to write Compugen's experimental extended profile format
 *
 * Returns:  (void)
 */
void
WriteProfile(FILE *fp, struct plan7_s *hmm, int do_xsw)
{
  int k;			/* position in model      */
  int x;			/* symbol index           */
  int sc;			/* a score to print       */
  float nx;			/* expected # of symbol x */
  int gap, len, qgap, qlen;	/* penalties to charge    */
  
  P7Logoddsify(hmm, TRUE);

  /* GCG can't deal with long profiles. Their limit is 1000
   * positions. However, Compugen can. Therefore we warn,
   * but don't die.
   */
  if (hmm->M > 1000 && !do_xsw)
    Warn("Profile %s will have more than 1000 positions. GCG won't read it; Compugen will.",
	 hmm->name);

  /* Header information.
   * GCG will look for sequence type and length of model.
   * Other than this, nothing is parsed until we get to the 
   * Cons line that has a ".." on it.
   * Lines that begin with "!" are comments.
   */
  if (Alphabet_type == hmmAMINO)        fprintf(fp, "!!AA_PROFILE 1.0\n");
  else if (Alphabet_type == hmmNUCLEIC) fprintf(fp, "!!NA_PROFILE 1.0\n");
  else Die("No support for profiles with non-biological alphabets");

  if (Alphabet_type == hmmAMINO)        fprintf(fp, "(Peptide) ");
  else if (Alphabet_type == hmmNUCLEIC) fprintf(fp, "(Nucleotide) ");
  fprintf(fp, "HMMCONVERT v%s Length: %d %s|%s|%s\n",
	  PACKAGE_VERSION, hmm->M, hmm->name,
	  hmm->flags & PLAN7_ACC ? hmm->acc : "",
	  hmm->flags & PLAN7_DESC ? hmm->desc : "");
  
  /* Insert some HMMER-specific commentary
   */
  if (do_xsw)
    {
      fprintf(fp, "   Profile converted from a profile HMM using HMMER v%s emulation.\n", PACKAGE_VERSION);
      fprintf(fp, "   Compugen XSW extended profile format.\n");
      fprintf(fp, "   Use -model=xsw.model -nonor -noave -gapop=10 -gapext=1 -qgapop=10 -qgapext=1\n");
      fprintf(fp, "      with om on the Compugen BIC to get the closest approximation to HMMER bit scores.\n");
      fprintf(fp, "   WARNING: There is a loss of information in this conversion.\n");
      fprintf(fp, "      Neither the scores nor even the rank order of hits will be precisely\n");
      fprintf(fp, "      preserved in a comparison of HMMER hmmsearch to GCG profilesearch.\n");
      fprintf(fp, "      The profile score is an approximation of the (single-hit) HMMER score.\n\n");
    }
  else
    {
      fprintf(fp, "   Profile converted from a profile HMM using HMMER v%s emulation.\n", PACKAGE_VERSION);
      fprintf(fp, "   Use -nonor -noave -gap=10 -len=1 with profilesearch and friends\n");
      fprintf(fp, "      to get the closest approximation to HMMER bit scores.\n");
      fprintf(fp, "   WARNING: There is a loss of information in this conversion.\n");
      fprintf(fp, "      Neither the scores nor even the rank order of hits will be precisely\n");
      fprintf(fp, "      preserved in a comparison of HMMER hmmsearch to GCG profilesearch.\n");
      fprintf(fp, "      The profile score is an approximation of the (single-hit) HMMER score.\n\n");
    }


  /* Do the CONS line, which gives the valid IUPAC symbols and their order
   */
  fprintf(fp, "Cons");
  for (x = 0; x < Alphabet_iupac; x++)
    fprintf(fp, "    %c ", Alphabet[x]);
  if (do_xsw)
    fprintf(fp, "  Gap   Len  QGap  Qlen ..\n"); 
  else
    fprintf(fp, "  Gap   Len ..\n");
 
  /* Now, the profile; for each position in the HMM, write a line of profile.
   */
  for (k = 1; k <= hmm->M; k++)
    {
				/* GCG adds some indexing as comments */
      if ((k-1)%10 == 0 && k > 10)
	fprintf(fp, "! %d\n", k);

				/* find consensus residue by max prob */
      x = FArgMax(hmm->mat[k], Alphabet_size);
      fprintf(fp, " %c  ", Alphabet[x]);
				/* generate emission score profile;
				 * Profiles are scaled by a factor of 100 
				 */
      for (x = 0; x < Alphabet_iupac; x++)
	{
	  sc = hmm->msc[x][k];
	  if (k < hmm->M) sc += hmm->tsc[TMM][k];
	  sc = sc * 100 / INTSCALE;
	  fprintf(fp, "%5d ", sc);
	}
				/* Generate gap open, gap extend penalties;
				   note we will force profilesearch to weights of 10, 1,
				   and that GCG profile values are percentages
				   of these base penalties, 0..100.*/
				/* gap open (insertion)*/
      if (k > 1)
	{
	  gap = -1 * (hmm->tsc[TMI][k-1] + hmm->tsc[TIM][k-1] - hmm->tsc[TMM][k-1] - hmm->tsc[TII][k-1]);
	  gap = gap * 100 / (10.0 * INTSCALE);
	}
      else gap = 100;		/* doesn't matter because GAP_1 is never used */

				/* gap extend (insertion)*/
      if (k > 1)
	{
	  len = -1 * hmm->tsc[TII][k-1];
	  len = len * 100 / (1.0 * INTSCALE);
	}
      else len = 100;		/* again, doesn't matter because LEN_1 is never used */

				/* gap open (deletion) */
      if (k > 1)
	{
	  qgap = -1 * (hmm->tsc[TDM][k-1] + hmm->tsc[TMD][k-1] - hmm->tsc[TMM][k-1] - hmm->tsc[TDD][k-1]);
	  qgap = qgap * 100 / (10.0 * INTSCALE);
	}
      else qgap = 100;
				/* gap extend (deletion) */
      if (k > 1)
	{
	  qlen = -1 * hmm->tsc[TDD][k-1];
	  qlen = qlen * 100 / (1.0 * INTSCALE);
	}
      else qlen = 100;

      
      if (do_xsw)
	fprintf(fp, "%5d %5d %5d %5d\n", gap, len, qgap, qlen);
      else
	fprintf(fp, "%5d %5d\n", gap, len); /* assume insertions >= deletions */
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
      
