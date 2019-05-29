#ifndef h4BUILD_INCLUDED
#define h4BUILD_INCLUDED

#include "esl_alphabet.h"
#include "esl_msa.h"

#include "h4_prior.h"
#include "h4_profile.h"


typedef struct h4_build_config_s {
  int       arch_strategy;            // h4_BUILD_ARCH_{GIVEN | RULES}; default RULES
  float     fragthresh;               // seq fragment if aspan/alen < fragthresh
  float     symfrac;                  // consensus column if nres/(nres+ngap) >= symfrac
  int       nsamp;                    // if nseq > nsamp, determine cons by sampling
  float     maxfragdens;              // ... unless nfrag/nseq > maxfragdens, then use all <nseq>
  uint64_t  seed;                     // ... use this RNG64 seed in taking reproducible sample

  int       wgt_strategy;             // h4BUILD_WGT_{NONE | GIVEN | PB | GSC | BLOSUM}; default PB
  float     wid;                      // for _BLOSUM: single linkage clustering threshold

  int       effn_strategy;            // h4BUILD_EFFN_{NONE | GIVEN | EWGT}; default EWGT
  float     re_target;                // for _EWGT: relative entropy target for eweights
  float     re_sigma;                 // for _EWGT: minimum total score param for eweighting short models
  float     effn_set;                 // for _NONE: effective seq number set by caller; or -1.

  H4_PRIOR *pri;

  int       stop_early;                // TRUE to stop after counting, return observed counts
  const ESL_ALPHABET  *abc;                
} H4_BUILD_CONFIG;


#define h4BUILD_FRAGTHRESH    0.5
#define h4BUILD_SYMFRAC       0.5
#define h4BUILD_NSAMP         10000    
#define h4BUILD_MAXFRAGDENS   0.5      
#define h4BUILD_RNGSEED       42       
#define h4BUILD_WID           0.62

                                       // Entropy weighting target defaults (in bits):
#define h4BUILD_ETARG_PRT     0.59     // .. for protein; from work of Steve Johnson.
#define h4BUILD_ETARG_NUC     0.45     // .. for DNA/RNA; from work of Travis Wheeler.
#define h4BUILD_ETARG_OTH     1.0      // .. for other alphabets; arbitrary choice.
#define h4BUILD_ESIGMA        45.0     // Entropy weighting parameter for short MSAs

                                       // Exclusive choices for <arch_strategy>: 
#define h4BUILD_ARCH_GIVEN    1        // .. use consensus columns specified by user's input MSA
#define h4BUILD_ARCH_RULES    2        // .. use fragthresh/symfrac rules (default)

                                       // Exclusive choices for <wgt_strategy>:
#define h4BUILD_WGT_NONE      1        // .. no relative weights (all seqs = 1)
#define h4BUILD_WGT_GIVEN     2        // .. use weights specified by user's input MSA
#define h4BUILD_WGT_PB        3        // .. Henikoff position-based weights (default)
#define h4BUILD_WGT_GSC       4        // .. Gerstein-Sonnhammer-Chothia tree weights
#define h4BUILD_WGT_BLOSUM    5        // .. BLOSUM cluster weights

                                       // Exclusive choices for <effn_strategy>:
#define h4BUILD_EFFN_NONE     1        // .. just use <nseq>
#define h4BUILD_EFFN_GIVEN    2        // .. use value specified by user
#define h4BUILD_EFFN_EWGT     3        // .. use entropy weighting algorithm (default)


extern int h4_Build(const H4_BUILD_CONFIG *cfg, ESL_MSA *msa, H4_PROFILE **ret_hmm, char *errbuf);

extern H4_BUILD_CONFIG *h4_build_config_Create(const ESL_ALPHABET *abc);
extern void             h4_build_config_Destroy(H4_BUILD_CONFIG *cfg);

#endif /*h4BUILD_INCLUDED*/
