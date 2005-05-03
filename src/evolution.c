/* evolution.c
 * Routines for phylogenetic "extrapolation" of profile HMMs.
 * 
 * SRE, Sun Oct 31 08:13:33 2004 [St. Louis]
 * SVN $Id$
 */

#include <math.h>

#include "funcs.h"


/* Function:  EvolveOneTransitionVector()
 * Incept:    SRE, Sun Oct 31 08:15:02 2004 [St. Louis]
 *
 * Purpose:   Implementation of Rivas' evolutionary modeling of a model transition
 *            process [Rivas05].
 *            
 *            Given: at a nontrivial time t* (ts), a transition vector
 *            q* (qs). Also given, assertions of the transition vector (q0)
 *            at time 0, and the stationary vector q_\infty (qz) at time 
 *            infinity.
 *            
 *            Each component 0..i..n of the three vectors must obey one
 *            of these three conditions:
 *                q0[i] < qs[i] < qz[i]
 *                q0[i] > qs[i] > qz[i]
 *                q0[i] = qs[i] = qz[i]
 *            (because all we're really doing is fitting each component
 *            independently to an exponential decay/accumulation curve.)
 *            (Note that this is not an "evolutionary model" in the sense that
 *            multiplicativity will not hold.)
 *            
 *            Each probability vector must also be properly normalized
 *            (\sum = 1). The code does not check this.
 *            
 *            Calculate: for some given time t, calculate and return a new
 *            transition vector, q.
 *            
 * Method:    [Rivas05] expresses the calculation in matrix form; but can
 *            also think of it as fitting independent exponential curves
 *            starting from q0[i], asymptoting to qz[i], while passing through
 *            qs[i].
 *            
 * Args:      qs:   probability vector at time t*
 *            ts:   time t* (typically set to 1.0 by convention, but could
 *                  also be something in units of substitutions/site or PAMs)
 *            n:    # of elements in q vectors.
 *            q0:   probability vector at time 0
 *            qz:   probability vector at time infinity
 *            t:    desired time.
 *            q:    RETURN: probability vector at time t.
 *
 * Returns:   On success: returns 1, and q contains the vector at time t; 
 *
 *            If a component of the vector doesn't obey one of the necessary
 *            conditions, returns 0.
 *
 * Xref:      [Rivas05]; STL8 p.117
 */
int
EvolveOneTransitionVector(float *qs, float ts, int n, float *q0, float *qz, float t, float *q)
{
  int   i;
  float wt;

  /* Check that the conditions hold.
   */
  for (i = 0; i < n; i++)
    {
      if (! (q0[i] <  qs[i] <  qz[i] ||
	     q0[i] >  qs[i] >  qz[i] ||
	     q0[i] == qs[i] == qz[i]))
	return 0;
    }
  
  /* Temporarily use q[i] to hold the eigenvalues (decay constants).
   * [Rivas05], eqn 110.
   */
  for (i = 0; i < n; i++)
    q[i] = (1./ts) * log( (qs[i] - qz[i]) / (q0[i] - qz[i]));
  
  /* Calculate the normalization constant wt.
   * [Rivas05], eqn 113.
   */
  wt = 0.;
  for (i = 0; i < n; i++)
    wt += exp( t * q[i] ) * (q0[i]-qz[i]);
 
  /* Calculate q.
   * [Rivas05], eqn 114.
   */
  for (i = 0; i < n; i++)
    q[i] = (qz[i] + exp(t*q[i])*(q0[i]-qz[i])) / (1. + wt);
  return;
}

/************************************************************
 * @LICENSE@
 ************************************************************/

