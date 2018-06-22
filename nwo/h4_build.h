#ifndef h4_BUILD_INCLUDED
#define h4_BUILD_INCLUDED

enum h4_archchoice_e {  h4_ARCH_SYMFRAC = 0, h4_ARCH_GIVEN = 1 };

typedef struct h4_build_config_s {
  /* Model architecture */
  enum h4_archchoice_e arch_strategy;
  float                symfrac;
  float                fragthresh;
} H4_BUILD_CONFIG;

extern H4_BUILD_CONFIG *h4_build_config_Create(void);
extern void             h4_build_config_Destroy(H4_BUILD_CONFIG *cfg);

#endif /*h4_BUILD_INCLUDED*/
