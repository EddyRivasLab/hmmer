/* General routines used throughout HMMER.
 * 
 * Internally, HMMER profiles work in scaled integer log-odds (SILO)
 * scores.  The dynamic range of SILO scores is controlled by
 * p7_INTSCALE at compile time (p7_INTSCALE defaults to
 * 1000). Externally, HMMER reports real-numbered bit scores.
 * The code refers to "SILO score" and "bit score" to differentiate.
 * 
 * SRE, Fri Jan 12 13:19:38 2007 [Janelia] [Franz Ferdinand, eponymous]
 * SVN $Id$
 */

#include "p7_config.h"

#include <math.h>
#include <float.h>

#include "easel.h"
#include "hmmer.h"

/* Function:  p7_banner()
 * Synopsis:  print standard HMMER application output header
 * Incept:    SRE, Wed May 23 10:45:53 2007 [Janelia]
 *
 * Purpose:   Print the standard HMMER command line application banner
 *            to <fp>, constructing it from <progname> (the name of the
 *            program) and a short one-line description <banner>.
 *            For example, 
 *            <p7_banner(stdout, "hmmsim", "collect profile HMM score distributions");>
 *            might result in:
 *            
 *            \begin{cchunk}
 *            # hmmsim :: collect profile HMM score distributions
 *            # HMMER 3.0 (May 2007)
 *            # Copyright (C) 2004-2007 HHMI Janelia Farm Research Campus
 *            # Freely licensed under the Janelia Software License.
 *            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *            \end{cchunk}
 *              
 *            <progname> would typically be an application's
 *            <argv[0]>, rather than a fixed string. This allows the
 *            program to be renamed, or called under different names
 *            via symlinks. Any path in the <progname> is discarded;
 *            for instance, if <progname> is "/usr/local/bin/hmmsim",
 *            "hmmsim" is used as the program name.
 *            
 * Note:    
 *    Needs to pick up preprocessor #define's from p7_config.h,
 *    as set by ./configure:
 *            
 *    symbol          example
 *    ------          ----------------
 *    HMMER_VERSION   "3.0"
 *    HMMER_DATE      "May 2007"
 *    HMMER_COPYRIGHT "Copyright (C) 2004-2007 HHMI Janelia Farm Research Campus"
 *    HMMER_LICENSE   "Freely licensed under the Janelia Software License."
 *
 * Returns:   (void)
 */
void
p7_banner(FILE *fp, char *progname, char *banner)
{
  char *appname = NULL;

  if (esl_FileTail(progname, FALSE, &appname) != eslOK) appname = progname;

  fprintf(fp, "# %s :: %s\n", appname, banner);
  fprintf(fp, "# HMMER %s (%s)\n", HMMER_VERSION, HMMER_DATE);
  fprintf(fp, "# %s\n", HMMER_COPYRIGHT);
  fprintf(fp, "# %s\n", HMMER_LICENSE);
  fprintf(fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  if (appname != NULL) free(appname);
  return;
}



/* Function:  p7_SILO2Lod()
 * Synopsis:  Convert a scaled integer score to a lod score.
 * Incept:    SRE, Wed Jul 18 08:52:19 2007 [Janelia]
 *
 * Purpose:   Convert a SILO (scaled integer log-odds) score
 *            to a lod score.
 *            
 *            HMMER3 normally works internally in floating-point log
 *            odds probabilities, in nats. Some alternate/optimized
 *            implementations of alignment algorithms (Viterbi,
 *            Forward) work in scaled integer scores called SILO
 *            scores. These routines must convert their final answer
 *            back to lod score form before returning it to H3.
 *            
 *            Beside simply descaling and casting to a float, the
 *            conversion needs to be careful to convert the integer
 *            <p7_IMPOSSIBLE> to $\-infty$; and indeed, anything
 *            within shouting range of <p7_IMPOSSIBLE>, because we may
 *            have added a positive score to an impossible value.
 *
 */
float
p7_SILO2Lod(int silo)
{
  if (silo <= p7_IMPOSSIBLE + 10 * p7_INTSCALE) return -eslINFINITY; /* anything within 10 nats of impossible is impossible */
  else return (float) silo / p7_INTSCALE;
}


/* Function:  p7_AminoFrequencies()
 * Incept:    SRE, Fri Jan 12 13:46:41 2007 [Janelia]
 *
 * Purpose:   Fills a vector <f> with amino acid background frequencies,
 *            in [A..Y] alphabetic order, same order that Easel digital
 *            alphabet uses. Caller must provide <f> allocated for at
 *            least 20 floats.
 *            
 *            The frequencies here were counted over SwissProt 34,
 *            21.2M residues.
 *
 * Returns:   <eslOK> on success.
 *
 * Xref:      STL11/125
 */
int
p7_AminoFrequencies(float *f)
{
  f[0]  = 0.075520;			/* A */
  f[1]  = 0.016973;			/* C */
  f[2]  = 0.053029;			/* D */
  f[3]  = 0.063204;			/* E */
  f[4]  = 0.040762;			/* F */
  f[5]  = 0.068448;			/* G */
  f[6]  = 0.022406;			/* H */
  f[7]  = 0.057284;			/* I */
  f[8]  = 0.059398;			/* K */
  f[9]  = 0.093399;			/* L */
  f[10] = 0.023569;			/* M */
  f[11] = 0.045293;			/* N */
  f[12] = 0.049262;			/* P */
  f[13] = 0.040231;			/* Q */
  f[14] = 0.051573;			/* R */
  f[15] = 0.072214;			/* S */
  f[16] = 0.057454;			/* T */
  f[17] = 0.065252;			/* V */
  f[18] = 0.012513;			/* W */
  f[19] = 0.031985;			/* Y */
  return eslOK;
}



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
