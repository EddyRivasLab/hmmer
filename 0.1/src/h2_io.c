/* Input/output of HMMER3 HMMs in HMMER2 save file formats:
 * for backwards compatibility.
 * 
 * Contents:
 * 
 * 
 * SRE, Sun Dec 21 10:26:22 2008 [Zaragoza]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "hmmer.h"
#include "easel.h"

static void h2_multiline(FILE *fp, const char *pfx, char *s);
static void printprob(FILE *fp, int fieldwidth, float p, float null);

/* Function:  p7_h2io_WriteASCII()
 * Synopsis:  Write an H3 HMM in HMMER2 compatible format
 * Incept:    SRE, Sun Dec 21 10:26:53 2008 [Janelia]
 *
 * Purpose:   Write HMM <hmm> to stream <fp> in HMMER2 ASCII save
 *            file format.
 *            
 *            HMMER2 saved the null model and the search configuration
 *            (local vs. glocal, for example) as part of its HMM file;
 *            H3 only saves the core HMM. The HMMER2 file is created
 *            for HMMER2's default ``ls mode'' (glocal) with default
 *            null model transitions and default special state
 *            transitions (NECJ).
 *            
 *            Optional statistical calibration and alignment checksum
 *            are not written, because for these H3 and H2 differ too
 *            much.
 *
 * Args:      fp   - stream to write save file format to
 *            hmm  - HMM to save
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error. <eslEINVAL> if <hmm> can't
 *            be converted; for example, if it is not in a protein 
 *            or nucleic acid alphabet (H2 requires biosequence in its
 *            save files). 
 */
int
p7_h2io_WriteASCII(FILE *fp, P7_HMM *hmm)
{
  P7_BG *bg;			/* H2 saves null model in HMM file   */
  int    k;                     /* counter for nodes                 */
  int    x;                     /* counter for symbols               */
  int    ts;			/* counter for state transitions     */
  float  pmove,ploop;		/* default H2 null model transitions */
  int    status;

  if ((bg = p7_bg_Create(hmm->abc)) == NULL) { status = eslEMEM; goto ERROR; }


  fprintf(fp, "HMMER2.0  [converted from %s]\n", HMMER_VERSION);  /* magic header */
  fprintf(fp, "NAME  %s\n", hmm->name);
  if (hmm->acc)  fprintf(fp, "ACC   %s\n", hmm->acc);
  if (hmm->desc) fprintf(fp, "DESC  %s\n", hmm->desc);
  fprintf(fp, "LENG  %d\n", hmm->M);

  if      (hmm->abc->type == eslAMINO)   fprintf(fp, "ALPH  Amino\n");  
  else if (hmm->abc->type == eslDNA)     fprintf(fp, "ALPH  Nucleic\n");  
  else if (hmm->abc->type == eslRNA)     fprintf(fp, "ALPH  Nucleic\n");  
  else ESL_XEXCEPTION(eslEINVAL, "Only protein, DNA, RNA HMMs can be saved in H2 format");

  fprintf(fp, "RF    %s\n", (hmm->flags & p7H_RF)  ? "yes" : "no");
  fprintf(fp, "CS    %s\n", (hmm->flags & p7H_CS)  ? "yes" : "no");
  fprintf(fp, "MAP   %s\n", (hmm->flags & p7H_MAP) ? "yes" : "no");

  if (hmm->comlog != NULL)  h2_multiline(fp, "COM   ",     hmm->comlog);
  if (hmm->nseq   != -1)    fprintf     (fp, "NSEQ  %d\n", hmm->nseq);
  if (hmm->ctime  != NULL)  fprintf     (fp, "DATE  %s\n", hmm->ctime); 

  /* Checksum is not written; H2 and H3 use different checksum algorithms */

  if (hmm->flags & p7H_GA) fprintf(fp, "GA    %.1f %.1f\n", hmm->cutoff[p7_GA1], hmm->cutoff[p7_GA2]);
  if (hmm->flags & p7H_TC) fprintf(fp, "TC    %.1f %.1f\n", hmm->cutoff[p7_TC1], hmm->cutoff[p7_TC2]);
  if (hmm->flags & p7H_NC) fprintf(fp, "NC    %.1f %.1f\n", hmm->cutoff[p7_NC1], hmm->cutoff[p7_NC2]);

  /* in H3, the HMM does not include NECJ; these are part of the profile.
   * for emulating H2 output, assume default LS config */
  fputs("XT     ", fp);
  pmove = ( (hmm->abc->type == eslAMINO) ?   1./351. :    1./1001.); 
  ploop = ( (hmm->abc->type == eslAMINO) ? 350./351. : 1000./1001.); 
  printprob(fp, 6, pmove, 1.0);	/* NB */
  printprob(fp, 6, ploop, 1.0);	/* NN */
  printprob(fp, 6,   0.5, 1.0);	/* EC */
  printprob(fp, 6,   0.5, 1.0);	/* EJ */
  printprob(fp, 6, pmove, 1.0);	/* CT */
  printprob(fp, 6, ploop, 1.0);	/* CC */
  printprob(fp, 6, pmove, 1.0);	/* JB */
  printprob(fp, 6, ploop, 1.0);	/* JJ */
  fputc('\n', fp);

  /* Save the default H2 null model transitions, not H3's null model transitions  */
  fprintf(fp, "NULT   ");
  printprob(fp, 6, ploop, 1.0);	/* 1-p1 */
  printprob(fp, 6, pmove, 1.0);	/* p1   */
  fputc('\n', fp);
  
  /* but null emissions really are the H3 null model emissions */
  fputs("NULE   ", fp);
  for (x = 0; x < hmm->abc->K; x++)
    printprob(fp, 6, bg->f[x], 1./(float)hmm->abc->K);
  fputc('\n', fp);

  /* Don't save stats; H3 local alignment stats are different from H2 calibration */
  
  /* The main model section */
  fprintf(fp, "HMM      ");
  for (x = 0; x < hmm->abc->K; x++) fprintf(fp, "  %c    ", hmm->abc->sym[x]);
  fprintf(fp, "\n");
  fprintf(fp, "       %6s %6s %6s %6s %6s %6s %6s %6s %6s\n",
          "m->m", "m->i", "m->d", "i->m", "i->i", "d->m", "d->d", "b->m", "m->e");

  /* Print HMM parameters (main section of the save file) */
  fprintf(fp, "      ");
  printprob(fp, 6, 1.-hmm->t[0][p7H_MD], 1.0);
  fprintf(fp, " %6s", "*");
  printprob(fp, 6, hmm->t[0][p7H_MD], 1.0);
  fputc('\n', fp);

  for (k = 1; k <= hmm->M; k++)
    {
				/* Line 1: k, match emissions, map */
      fprintf(fp, " %5d ", k);
      for (x = 0; x < hmm->abc->K; x++) 
	printprob(fp, 6, hmm->mat[k][x], bg->f[x]);
      if (hmm->flags & p7H_MAP) fprintf(fp, " %5d", hmm->map[k]);
      fputc('\n', fp);
				/* Line 2: RF and insert emissions */
      fprintf(fp, " %5c ", hmm->flags & p7H_RF ? hmm->rf[k] : '-');
      for (x = 0; x < hmm->abc->K; x++) 
	printprob(fp, 6, ((k < hmm->M) ? hmm->ins[k][x] : 0.0), bg->f[x]);
      fputc('\n', fp);
				/* Line 3: CS and transition probs */
      fprintf(fp, " %5c ", hmm->flags & p7H_CS ? hmm->cs[k] : '-');
      for (ts = 0; ts < 7; ts++)
	printprob(fp, 6, ((k < hmm->M) ? hmm->t[k][ts] : 0.0), 1.0);
      printprob(fp, 6, ((k==1)     ? hmm->t[0][p7H_MM] : 0.0), 1.0);
      printprob(fp, 6, ((k<hmm->M) ? 0.0: 1.0), 1.0);
      fputc('\n', fp);

    }
  fputs("//\n", fp);

  p7_bg_Destroy(bg);
  return eslOK;

 ERROR:
  p7_bg_Destroy(bg);
  return status;
}

/* h2_multiline()
 * SRE, Sun Dec 21 10:40:23 2008 [Zaragoza]
 * 
 * Used to print the command log to HMMER2 ASCII save files.
 * H3 records command numbers in brackets, as in "COM [1] hmmbuild ..."
 * H2 just records commands, as in "COM   hmmbuild ...".
 * Compare p7_hmmfile.c::multiline().
 *
 * Given a record (like the comlog) that contains 
 * multiple lines, print it as multiple lines with
 * a given prefix. e.g.:
 *           
 * given:   "COM   ", "foo\nbar\nbaz"
 * print:   COM   foo
 *          COM   bar
 *          COM   baz
 *
 * If <s> is NULL, no-op. Otherwise <s> must be a <NUL>-terminated
 * string.  It does not matter if it ends in <\n> or not. <pfx>
 * must be a valid <NUL>-terminated string; it may be empty.
 *           
 * Args:     fp:   FILE to print to
 *           pfx:  prefix for each line
 *           s:    line to break up and print; tolerates a NULL
 */
static void
h2_multiline(FILE *fp, const char *pfx, char *s)
{
  char *sptr  = s;
  char *end   = NULL;
  int   n     = 0;

  do {
    end = strchr(sptr, '\n');

    if (end != NULL) 		             /* if there's no \n left, end == NULL */
      {
	n = end - sptr;	                     /* n chars exclusive of \n */
	fprintf(fp, "%s ", pfx);
	fwrite(sptr, sizeof(char), n, fp);   /* using fwrite lets us write fixed # of chars   */
	fprintf(fp, "\n");                   /* while writing \n w/ printf allows newline conversion */
	sptr += n + 1;	                     /* +1 to get past \n */
      } 
    else 
      {
	fprintf(fp, "%s %s\n", pfx, sptr); /* last line */
      }
  } while (end != NULL  && *sptr != '\0');   /* *sptr == 0 if <s> terminates with a \n */
}


/* printprob()
 * Print a probability (with a leading space), formatted
 * for an H2 ASCII save file.
 */
static void
printprob(FILE *fp, int fieldwidth, float p, float null)
{
  if      (p == 0.0)                fprintf(fp, " %*s", fieldwidth, "*");
  else if (null == 1.0 && p == 1.0) fprintf(fp, " %*d", fieldwidth, 0);
  else                              fprintf(fp, " %*d", fieldwidth, (int) floor(0.5 + 1442.695 * log(p/null)));
}

