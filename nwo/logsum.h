#ifndef h4LOGSUM_INCLUDED
#define h4LOGSUM_INCLUDED
#include "h4_config.h"

extern int   h4_logsum_Init(void);
extern int   h4_logsum_Reinit(void);
extern float h4_logsum(float a, float b);

extern int   h4_logsum_IsSlowExact(void);
extern int   h4_logsum_InitMax(void);

#endif /*h4LOGSUM_INCLUDED*/

