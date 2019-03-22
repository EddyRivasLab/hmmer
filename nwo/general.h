#ifndef h4GENERAL_INCLUDED
#define h4GENERAL_INCLUDED

#include "esl_getopts.h"

extern ESL_GETOPTS *h4_CreateDefaultApp(ESL_OPTIONS *options, int nargs, int argc, char **argv, char *banner, char *usage);
extern int          h4_AminoFrequencies(float *f);

#endif // h4GENERAL_INCLUDED
