#ifndef CUSTOMSTRUCTSH_INCLUDED
#define CUSTOMSTRUCTSH_INCLUDED

#ifdef SLOW
#include "defaultstructs.h"
#endif
#ifdef FAST
#include "defaultstructs.h"
#endif
#ifdef ALTIVEC
#include "altivecstructs.h"
#endif
#ifdef LOCAL
#include "localstructs.h"
#endif

#endif /*LOCALSTRUCTSH_INCLUDED*/
