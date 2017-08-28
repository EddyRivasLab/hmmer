#ifndef p7H2_IO_INCLUDED
#define p7H2_IO_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#include "base/p7_hmm.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif

extern int   p7_h2io_WriteASCII(FILE *fp, P7_HMM *hmm);

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*p7H2_IO_INCLUDED*/

