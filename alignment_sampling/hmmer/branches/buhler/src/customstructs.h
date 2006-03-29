#ifndef CUSTOMSTRUCTSH_INCLUDED
#define CUSTOMSTRUCTSH_INCLUDED

#ifdef IMPL_SLOW
#include "defaultstructs.h"
#endif
#ifdef IMPL_FAST
#include "defaultstructs.h"
#endif
#ifdef IMPL_ALTIVEC
#include "altivecstructs.h"
#endif
#ifdef IMPL_LOCAL
#include "localstructs.h"
#endif
#ifdef IMPL_JDB
#include "jdbstructs.h"
#endif

#endif /*LOCALSTRUCTSH_INCLUDED*/
