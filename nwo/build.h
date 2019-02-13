#ifndef h4BUILD_INCLUDED
#define h4BUILD_INCLUDED

#include "esl_alphabet.h"
#include "esl_msa.h"

#include "h4_prior.h"
#include "h4_profile.h"

enum h4_archchoice_e {
  h4_ARCH_SYMFRAC = 0,
  h4_ARCH_GIVEN   = 1
};

enum h4_wgtchoice_e {
  h4_WGT_PB     = 0,
  h4_WGT_GIVEN  = 1,
  h4_WGT_GSC    = 2,
  h4_WGT_BLOSUM = 3,
  h4_WGT_NONE   = 4
};

typedef struct h4_build_config_s {
  /* Model architecture */
  enum h4_archchoice_e arch_strategy;
  enum h4_wgtchoice_e  wgt_strategy;
  float                symfrac;
  float                fragthresh;
  float                wid;

  H4_PRIOR            *pri;

  int                  stop_after_counting;
  const ESL_ALPHABET  *abc;                
} H4_BUILD_CONFIG;


extern int h4_Build(H4_BUILD_CONFIG *cfg, ESL_MSA *msa, H4_PROFILE **ret_hmm, char *errbuf);

extern H4_BUILD_CONFIG *h4_build_config_Create(const ESL_ALPHABET *abc);
extern void             h4_build_config_Destroy(H4_BUILD_CONFIG *cfg);

#endif /*h4BUILD_INCLUDED*/
