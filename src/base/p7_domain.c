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
  for (d = 0; d < ndom; d++) {
    dcl[d].ad = NULL;
    dcl[d].scores_per_pos = NULL;
  }
  return dcl;

 ERROR:
  return NULL;
}

void
p7_domain_Destroy(P7_DOMAIN *dcl, int ndom)
{
  int d;
  if (dcl) {
    for (d = 0; d < ndom; d++) {
      if (dcl[d].ad) p7_alidisplay_Destroy(dcl[d].ad);
      if (dcl[d].scores_per_pos) free(dcl[d].scores_per_pos);
    }
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
  dcl->scores_per_pos = NULL; //fixme. TJW is using this; SRE is not testing it.

  if ((status = p7_alidisplay_TestSample(rng, alen, &(dcl->ad))) != eslOK) return status;
  if ((status = p7_alidisplay_Serialize (dcl->ad))               != eslOK) return status;

  return eslOK;
}

int
p7_domain_Validate(const P7_DOMAIN *dcl, char *errbuf)
{
  if (dcl->iae < 1 || dcl->ibe < 1) ESL_FAIL(eslFAIL, errbuf, "bad iae/ibe");
  if (dcl->kae < 1 || dcl->kbe < 1) ESL_FAIL(eslFAIL, errbuf, "bad kae/kbe");
  if (dcl->ia  < 1 || dcl->ib  < 1) ESL_FAIL(eslFAIL, errbuf, "bad ia/ib");
  if (dcl->ka  < 1 || dcl->kb  < 1) ESL_FAIL(eslFAIL, errbuf, "bad ka/kb");

  if ( isnan(dcl->envsc) ||
       isnan(dcl->domcorrection) ||
       isnan(dcl->dombias) ||
       isnan(dcl->oasc) ||
       isnan(dcl->bitscore) ||
       isnan(dcl->lnP)) ESL_FAIL(eslFAIL, errbuf, "NaN");

  if (dcl->is_reported != 0 && dcl->is_reported != 1) ESL_FAIL(eslFAIL, errbuf, "is_reported: 0|1");
  if (dcl->is_included != 0 && dcl->is_included != 1) ESL_FAIL(eslFAIL, errbuf, "is_included: 0|1");

  if (dcl->ad == NULL) ESL_FAIL(eslFAIL, errbuf, "bad (NULL) ad");

  return p7_alidisplay_Validate(dcl->ad, errbuf);
}

int
p7_domain_Compare(const P7_DOMAIN *dcl1, const P7_DOMAIN *dcl2, float tol)
{
  if ( dcl1->iae != dcl2->iae ||  dcl1->ibe != dcl2->ibe ) return eslFAIL;
  if ( dcl1->kae != dcl2->kae ||  dcl1->kbe != dcl2->kbe ) return eslFAIL;
  if ( dcl1->ia  != dcl2->ia  ||  dcl1->ib  != dcl2->ib  ) return eslFAIL;
  if ( dcl1->ka  != dcl2->ka  ||  dcl1->kb  != dcl2->kb  ) return eslFAIL;

  if (esl_FCompare_old(dcl1->envsc,         dcl2->envsc,         tol) != eslOK) return eslFAIL;
  if (esl_FCompare_old(dcl1->domcorrection, dcl2->domcorrection, tol) != eslOK) return eslFAIL;
  if (esl_FCompare_old(dcl1->dombias,       dcl2->dombias,       tol) != eslOK) return eslFAIL;
  if (esl_FCompare_old(dcl1->oasc,          dcl2->oasc,          tol) != eslOK) return eslFAIL;
  if (esl_FCompare_old(dcl1->bitscore,      dcl2->bitscore,      tol) != eslOK) return eslFAIL;
  if (esl_DCompare_old(dcl1->lnP,           dcl2->lnP,           tol) != eslOK) return eslFAIL;

  if (dcl1->is_reported != dcl2->is_reported) return eslFAIL;
  if (dcl1->is_included != dcl2->is_included) return eslFAIL;

  return p7_alidisplay_Compare(dcl1->ad, dcl2->ad);
}



