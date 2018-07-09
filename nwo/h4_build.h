#ifndef h4_BUILD_INCLUDED
#define h4_BUILD_INCLUDED

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
  float                alen_occthresh;
  float                fragthresh;
  float                wid;
} H4_BUILD_CONFIG;

extern H4_BUILD_CONFIG *h4_build_config_Create(void);
extern void             h4_build_config_Destroy(H4_BUILD_CONFIG *cfg);

#endif /*h4_BUILD_INCLUDED*/
