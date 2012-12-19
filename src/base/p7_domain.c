/* Domain coordinates.
 * 
 * Each P7_HIT structure contains an array of 1 or more P7_DOMAIN
 * structures, with coordinates of each 'domain' alignment in the
 * target. In turn, each P7_DOMAIN contains a P7_ALIDISPLAY.
 * 
 * P7_TOPHITS => P7_HIT array => P7_DOMAIN array => P7_ALIDISPLAY
 * {--- p7_tophits.c -------}  {- p7_domain.c -} {- p7_alidisplay.c -}
 * 
 * Contents:
 *    1. P7_DOMAIN object (arrays of)
 *    2. Debugging and development tools
 *    3. Copyright and license information.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"

#include "base/p7_domain.h"
#include "base/p7_alidisplay.h"

/*****************************************************************
 * 1. P7_DOMAIN object (arrays of)
 *****************************************************************/

P7_DOMAIN *
p7_domain_Create(int ndom)
{
  P7_DOMAIN *dcl = NULL;
  int        d;
  int        status;

  ESL_ALLOC(dcl, sizeof(P7_DOMAIN) * ndom);
  for (d = 0; d < ndom; d++) dcl[d].ad = NULL;
  return dcl;

 ERROR:
  return NULL;
}

void
p7_domain_Destroy(P7_DOMAIN *dcl, int ndom)
{
  int d;
  if (dcl) {
    for (d = 0; d < ndom; d++)
      if (dcl[d].ad) p7_alidisplay_Destroy(dcl[d].ad);
    free(dcl);
  }
}

/*****************************************************************
 * 2. Debugging, development tools
 *****************************************************************/

int
p7_domain_TestSample(ESL_RANDOMNESS *rng, int alen, P7_DOMAIN *dcl)
{
  int status;

  dcl->iae           = 1+esl_rnd_Roll(rng, 100000);
  dcl->ibe           = 1+esl_rnd_Roll(rng, 100000);
  dcl->kae           = 1+esl_rnd_Roll(rng, 100000);
  dcl->kbe           = 1+esl_rnd_Roll(rng, 100000);
  dcl->ia            = 1+esl_rnd_Roll(rng, 100000);
  dcl->ib            = 1+esl_rnd_Roll(rng, 100000);
  dcl->ka            = 1+esl_rnd_Roll(rng, 100000);
  dcl->kb            = 1+esl_rnd_Roll(rng, 100000);
  dcl->envsc         = -1000. + 2000.*esl_random(rng);
  dcl->domcorrection = -1000. + 2000.*esl_random(rng);
  dcl->dombias       = -1000. + 2000.*esl_random(rng);
  dcl->oasc          = -1000. + 2000.*esl_random(rng);
  dcl->bitscore      = -1000. + 2000.*esl_random(rng);
  dcl->lnP           = -1000. + 2000.*esl_random(rng);
  dcl->is_reported   = esl_rnd_Roll(rng, 2);
  dcl->is_included   = esl_rnd_Roll(rng, 2);

  if ((status = p7_alidisplay_TestSample(rng, alen, &(dcl->ad))) != eslOK) return status;

  return eslOK;
}


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/


