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


/* Function: p7_Prob2SILO()
 * 
 * Purpose:  Convert a probability to a scaled integer log odds score. 
 *           Round to nearest integer (i.e. note use of +0.5 and floor())
 *           Return the score. 
 */
int
p7_Prob2SILO(float p, float null)
{
  if   (p == 0.0) return p7_IMPOSSIBLE;
  else            return (int) floor(0.5 + p7_INTSCALE * log(p/null));
}

/* Function:  p7_LL2SILO()
 * Incept:    SRE, Mon May  2 08:19:36 2005 [St. Louis]
 *
 * Purpose:   Convert a log likelihood to a scaled integer log odds score,
 *            rounded to nearest integer, given a <null> probability; 
 *            return the score. 
 *            
 *            Note that <ll> is a log(prob), but <null> is a probability.
 */
int
p7_LL2SILO(float ll, float null)
{
  int sc;
  sc = (int) floor(0.5 + p7_INTSCALE * (ll - log(null)));
  if (sc < p7_IMPOSSIBLE) sc = p7_IMPOSSIBLE;
  return sc;
}

/* Function: p7_SILO2Prob()
 * 
 * Purpose:  Convert a scaled integer lod score back to a probability;
 *           needs the null model probability (or 1.0) to do the conversion.
 */
float 
p7_SILO2Prob(int sc, float null)
{
  if (sc == p7_IMPOSSIBLE) return 0.;
  else                     return (null * exp((float) sc / p7_INTSCALE));
}

/* Function:  p7_SILO2Bitscore()
 * Incept:    SRE, Thu Feb  1 10:13:40 2007 [UA8018 St. Louis to Dulles]
 *
 * Purpose:   Convert an scaled integer lod score to a
 *            standard real-valued bit score, suitable for output.
 *
 */
float 
p7_SILO2Bitscore(int sc)
{
  return ((float) sc / p7_INTSCALE / eslCONST_LOG2);
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


static int ilogsum_lookup[p7_LOGSUM_TBL];

static void 
init_ilogsum(void)
{
  int i;
  for (i = 0; i < p7_LOGSUM_TBL; i++) 
    ilogsum_lookup[i] = (int) ((float) p7_INTSCALE * (log(1.+exp((float) -i/ (float) p7_INTSCALE))));
}

/* Function: p7_ILogsum()
 * 
 * Purpose:  Return the scaled integer log probability of
 *           the sum of two scores <s1> and <s2>, where
 *           <s1> and <s2> are also given as scaled log probabilities.
 *         
 *           $\log(\exp(s_1)+\exp(s_2)) = s_1 + \log(1 + \exp(s_2-s_1))$ for $s_1 > s_2$
 *           
 * Note:     For speed, builds a lookup table the first time it's
 *           called.  The table size <p7_LOGSUM_TBL> is set to 20000
 *           by default, in <p7_config.h>.
 *
 *           Because of the one-time initialization, we have to
 *           be careful in a multithreaded implementation... hence
 *           the use of <pthread_once()>, which forces us to put
 *           the initialization routine and the lookup table outside
 *           <p7_ILogsum()>. (Thanks to Henry Gabb at Intel for pointing
 *           out this problem.)
 *           
 * Args:     s1,s2 -- scaled integer log_2 probabilities to be summed
 *                    in probability space.
 *                    
 * Return:   scaled integer log_2 probability of the sum.
 */
int 
p7_ILogsum(int s1, int s2)
{
  int    diff;
#ifdef HMMER_THREADS
  static pthread_once_t firsttime = PTHREAD_ONCE_INIT;
  pthread_once(&firsttime, init_ilogsum);
#else
  static int firsttime = 1;
  if (firsttime) { init_ilogsum(); firsttime = 0; }
#endif

  diff = s1-s2;
  if      (diff >=  p7_LOGSUM_TBL) return s1;
  else if (diff <= -p7_LOGSUM_TBL) return s2;
  else if (diff > 0)               return s1 + ilogsum_lookup[diff];
  else                             return s2 + ilogsum_lookup[-diff];
} 


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
