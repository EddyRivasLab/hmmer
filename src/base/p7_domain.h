#ifndef p7DOMAIN_INCLUDED
#define p7DOMAIN_INCLUDED

#include "p7_config.h"

#include "base/p7_alidisplay.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
  
typedef struct p7_dom_s { 
  int            iae, ibe;	/* Envelope coords on sequence: iae..ibe in 1..L                              */
  int            kae, kbe;	/* Envelope coords on model:    kae..kbe in 1..M                              */
  int            ia,  ib;	/* Alignment coords on sequence                                               */
  int            ka,  kb;	/* Alignment coords on model                                                  */
  float          envsc;  	/* Forward score in envelope ienv..jenv; NATS; without null2 correction       */
  float          domcorrection;	/* null2 score when calculating a per-domain score; NATS                      */
  float          dombias;	/* FLogsum(0, log(bg->omega) + domcorrection): null2 score contribution; NATS */
  float          oasc;		/* optimal accuracy score (units: expected # residues correctly aligned)      */
  float          bitscore;	/* overall score in BITS, null corrected, if this were the only domain in seq */
  double         lnP;	        /* log(P-value) of the bitscore                                               */
  int            is_reported;	/* TRUE if domain meets reporting thresholds                                  */
  int            is_included;	/* TRUE if domain meets inclusion thresholds                                  */
  float         *scores_per_pos; /* score in BITS that each position in the alignment contributes to an overall viterbi score */
  P7_ALIDISPLAY *ad; 
} P7_DOMAIN;


/* 1. P7_DOMAIN object (arrays of) */
extern P7_DOMAIN *p7_domain_Create(int ndom);
extern void       p7_domain_Destroy(P7_DOMAIN *dcl, int ndom);

/* 2. Debugging, development tools */
extern int        p7_domain_TestSample(ESL_RANDOMNESS *rng, int alen, P7_DOMAIN *dcl);
extern int        p7_domain_Validate(const P7_DOMAIN *dcl, char *errbuf);
extern int        p7_domain_Compare(const P7_DOMAIN *dcl1, const P7_DOMAIN *dcl2, float tol);
#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif

#endif /*p7DOMAIN_INCLUDED*/

