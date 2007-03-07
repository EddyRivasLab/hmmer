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




/*****************************************************************
 * @LICENSE@
 *****************************************************************/
