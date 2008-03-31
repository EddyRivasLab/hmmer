/* p7_evalues.h
 * Calculation of E-values for HMMER scores.
 * 
 * A P7_EVINFO structure holds information about the expected score
 * distribution.
 * 
 * SRE, Mon Apr  3 09:46:18 2006
 * SVN $Id$
 */
#ifndef P7_EVALUES_INCLUDED
#define P7_EVALUES_INCLUDED

/* Possible distributions (for P7_EVINFO mode): */
#define p7_EXPTAIL 0
#define p7_GUMBEL  1

/* Structure: P7_EVINFO
 */
typedef struct {
  int    mode;		/* what type of distribution this is          */
  double mu;		/* location param (for gumbel or exponential) */
  double lambda;	/* scale param (gumbel or exponential)        */
  double tailmass;	/* 1.0, or fraction of tail that's fit        */

  int N;		/* this info was calibrated on <N> seqs...    */
  int L;		/*   ...of length <L>.                        */
} P7_EVINFO;


#endif /*P7_EVALUES_INCLUDED*/
/************************************************************
 * @LICENSE@
 ************************************************************/

