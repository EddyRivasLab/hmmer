#ifndef h4ZIGAR_INCLUDED
#define h4ZIGAR_INCLUDED

#include "h4_path.h"

extern int h4_zigar_Encode(const H4_PATH *pi, int z1, char **ret_zali);
extern int h4_zigar_Decode(const char *zali, H4_PATH *pi, int *opt_seqlen, int *opt_hmmlen);

#endif // h4ZIGAR_INCLUDED
