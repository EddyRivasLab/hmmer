/* Reference implementation of envelope definition.
 * 
 * Envelope definition step comes after we define the anchor set, and
 * before we do the alignment.
 * 
 * Contents:
 *    1. Envelope definition API call.
 *    2. Internal functions: the steps of envelope definition.
 *    3. Copyright and license information.
 */
// TODO:
//       write an example driver
//       get reference anchor set definition coded as designed

#include "p7_config.h"

#include "easel.h"

#include "base/p7_envelopes.h" 
#include "base/p7_coords2.h"
#include "base/p7_profile.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_asc_fwdback.h"


static int envcoords(P7_ENVELOPE *env, int D, const P7_REFMX *apd, float threshold);
static int glocality(P7_ENVELOPE *env, int D, const P7_REFMX *apd);
static int outcoords(P7_ENVELOPE *env, int D, const P7_REFMX *apd, float epsilon);
static int approxsc (P7_ENVELOPE *env, int D, const P7_REFMX *afd, const P7_PROFILE *gm);


int
p7_reference_Envelopes(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_COORD2 *anch, int D,
		       const P7_REFMX *apu, const P7_REFMX *apd,
		       P7_REFMX *afu, P7_REFMX *afd, P7_ENVELOPE *env)
{
  int d;
  int status;

  for (d = 0; d < D; d++)
    {
      env[d].i0    = anch[d].n1;
      env[d].k0    = anch[d].n2;
      env[d].flags = 0;
    }
  
  if ((status = envcoords(env, D, apd, 0.5))   != eslOK) goto ERROR;
  if ((status = glocality(env, D, apd))        != eslOK) goto ERROR;
  if ((status = outcoords(env, D, apd, 0.005)) != eslOK) goto ERROR;

  if (D > 0)
    {
      if ((status = approxsc(env, D, afd, gm)) != eslOK) goto ERROR;

      p7_refmx_Reuse(afu);
      p7_refmx_Reuse(afd);
      
      for (d = 0; d < D; d++)
	if (! (env[d].flags & p7E_ENVSC_APPROX))
	  {
	    status = p7_ReferenceASCForward( dsq+env[d].oea-1, env[d].oeb-env[d].oea+1, gm,
					     &(anch[d]), 1, afu, afd, &(env[d].env_sc));
	    if (status != eslOK) goto ERROR;

	    p7_refmx_Reuse(afu);
	    p7_refmx_Reuse(afd);
	  }
    }
  else 
    env[0].env_sc = P7R_XMX(afd, L, p7R_C) + gm->xsc[p7P_C][p7P_MOVE];

  /* ka,kb are not set yet. alignment does that part...
   * even though we already know for glocal alis that ka=1,kb=M.
   */
  return eslOK;
  
 ERROR:
  return status;

}

/* all specials are in <apd> */
/* threshold = 0.5 */
static int
envcoords(P7_ENVELOPE *env, int D, const P7_REFMX *apd, float threshold)
{
  float  phomology;
  int    d, i, iz;
  int    L = apd->L;
  
  for (d = 0; d < D; d++)
    {
      /* Determine ia, envelope start position;
       *  min_{i <= i0} P(x_i in domain d) >= threshold
       */
      phomology = 1.0;
      iz        = (d == 0 ? 0 : env[d-1].i0);
      for (i = env[d].i0 - 1; i >= iz; i--)
	{
	  phomology -= P7R_XMX(apd, i, p7R_B);   // now phomology = P(res i in domain d)
	  if (phomology < threshold) break;      // if i is not in domain...
	}
      env[d].ia = i+1;                           // but i+1 was, so ia = i+1.

      /* Determine ib, envelope end position;
       *   max_{i >= i0} P(x_i in domain d) >= threshold
       */
      phomology = 1.0;
      iz        = (d == D-1 ? L+1 : env[d+1].i0);
      for (i = env[d].i0 + 1; i < iz; i++)
	{
	  phomology -= P7R_XMX(apd, i, p7R_E);    // now phomology = P(x_i in dom_d)
	  if (phomology < threshold) break;       // if i is not in domain...
	}
      env[d].ib = i-1;		                  // but i-1 was, so ib = i-1.
    }

  return eslOK;
}



/* 
 *  For each domain d,
 *    pL = \sum_i=i0(d-1)..i0(d)-1 P(i,L)
 *    pG = \sum_i=i0(d-1)..i0(d)-1 P(i,G)
 *  and if pG >= pL, it's a glocal alignment.
 *  
 *  We do need the i0(d)-1 row, though it may not look like it at
 *  first. From a B on i0(d)-1 row, we need to use B->G->DDDDk0-1 path
 *  to be able to reach the anchor Mk on row i0(d). In decoding, we
 *  wing unfold the G->DD->Mk paths, so this path exists in a decoding
 *  matrix. Only for G, though; P(L, i0(d)) is 0.
 * 
 */
static int
glocality(P7_ENVELOPE *env, int D, const P7_REFMX *apd)
{
  float pL = 0.0;
  float pG = 0.0;
  int   d, i, iz;

  for (d = 0; d < D; d++)
    {
      iz = (d == 0 ? 0 : env[d-1].i0);
      for (i = iz; i < env[d].i0; i++)
	{
	  pL += P7R_XMX(apd, i, p7R_L);
	  pG += P7R_XMX(apd, i, p7R_G);
	}
      if (pG >= pL) env[d].flags |= p7E_IS_GLOCAL;
    }
  return eslOK;
}


/* "outer envelope", oea..oeb, for env score purposes */
/* for epsilon = 0.005 (allowing 0.005 of the probability
 * mass to be in paths that extend beyond the envelope on 
 * each side), we lose a max of 0.01 of the mass, which
 * is a log-odds score of ~0.01 in nats; negligible.
 */
static int
outcoords(P7_ENVELOPE *env, int D, const P7_REFMX *apd, float epsilon)
{
  int d, i, iz;
  float phomology;
  int   s;
  int   L = apd->L;
  int   leftchoked; 
  int   rightchoked;

  for (d = 0; d < D; d++)
    {
      leftchoked = rightchoked = FALSE;

      phomology = 1.0;
      iz = (d == 0 ? 0 : env[d-1].i0);
      s  = (d == 0 ? p7R_N : p7R_J);
      for (i = env[d].i0 - 1; i >= iz; i--)
	{
	  if (P7R_XMX(apd, i, s) >= 1.0 - epsilon) { leftchoked = TRUE; break; }
	  phomology -= P7R_XMX(apd, i, p7R_B);
	  if (phomology < epsilon) break;
	}
      env[d].oea = i+1;

      phomology = 1.0;
      iz = (d == D-1 ? L : env[d+1].i0-1);
      s  = (d == D-1 ? p7R_C : p7R_J);
      for (i = env[d].i0; i < iz; i++)
	{
	  if (P7R_XMX(apd, i, s) >= 1.0 - epsilon) { rightchoked = TRUE; break; }
	  phomology -= P7R_XMX(apd, i+1, p7R_E); // now phomology = P(x_i+1 in dom_d)
	  if (phomology < epsilon) break;        // if i+1 is not in domain...
	}
      env[d].oeb = i;		                 // but i=1 was, so oeb = i.

      if (leftchoked && rightchoked) 
	env[d].flags |= p7E_ENVSC_APPROX;
    }
  return eslOK;
}



/* If t_NN = t_JJ = t_CC, then we can often calculate the envelope
 * score of a single domain by an accurate approximation from the ASC
 * Forward matrix for the whole sequence.
 * 
 * Assume that we know that Xi and Yj are "choke points" flanking our
 * domain, such that >= 1.0-epsilon of the posterior path probability
 * flows through each of them.
 * 
 * For d=0, X=N, Y=J; for d=D-1, X=J, Y=C; for internal d, X=Y=J.
 * 
 * Then F(j,X) - F(i,X) = the sum of all path suffixes that start in
 * Xi and end at Xj.
 * 
 * Then, because t_NN = t_JJ = tCC, we just need to add in the log
 * transitions to get us from S to Xi via N's, and get us from Xj to T
 * via C's to get the envelope score. That's t_SN + i * t_NN, and
 * (L-j) * t_CC + t_CT. (And t_SN = 0).
 */
static int
approxsc(P7_ENVELOPE *env, int D, const P7_REFMX *afd, const P7_PROFILE *gm) 
{
  int d;
  int L = afd->L;
  int s1, s2;


  for (d = 0; d < D; d++)
    if (env[d].flags & p7E_ENVSC_APPROX)
      {
	s1 = (d == 0   ? p7R_N : p7R_J);
	s2 = (d == D-1 ? p7R_C : p7R_J);

	env[d].env_sc = 
	  P7R_XMX(afd, env[d].oeb,   s2) -
	  P7R_XMX(afd, env[d].oea-1, s1);

	env[d].env_sc += 
	  gm->xsc[p7P_N][p7P_LOOP] * (env[d].oea-1) +
	  gm->xsc[p7P_C][p7P_LOOP] * (L-env[d].oeb) +
	  gm->xsc[p7P_C][p7P_MOVE];
      }
  return eslOK;
}


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
