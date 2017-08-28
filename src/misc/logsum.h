#ifndef p7LOGSUM_INCLUDED
#define p7LOGSUM_INCLUDED

#include "p7_config.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif

extern int   p7_FLogsumInit(void);
extern float p7_FLogsum(float a, float b);

extern int   p7_logsum_InitMax(void);
extern int   p7_logsum_Reinit(void);
extern int   p7_logsum_IsSlowExact(void);
#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*p7LOGSUM_INCLUDED*/

