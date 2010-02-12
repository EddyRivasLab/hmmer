/* The "dummy" implementation: providing an "optimized" interface when
 * no optimization is available.
 * 
 * On Intel/AMD platforms, you should be using the SSE implementation
 * (impl_sse.c); on PowerPC platforms, you should be using the VMX
 * implementation (impl_vmx.c). The "dummy" implementation is provided
 * solely for portability to the minority of platforms that lack
 * either of the two industry SIMD instruction sets.
 * 
 * SRE, Sun Feb 17 10:02:22 2008 [Janelia]
 * SVN $Id$
 */
#ifdef  p7_IMPL_DUMMY
#ifndef P7_IMPL_DUMMY_INCLUDED
#define P7_IMPL_DUMMY_INCLUDED
#include "p7_config.h"
#include "esl_alphabet.h"
#include "hmmer.h"

/* An impl_ must provide P7_OPROFILE and P7_OMX.
 * The dummy impl_ just uses the generic forms.
 */
typedef P7_PROFILE P7_OPROFILE;
typedef P7_GMX     P7_OMX;


/* This is the API that an impl_ must provide: */
extern P7_OPROFILE *p7_oprofile_Create(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy(P7_OPROFILE *om);

extern P7_OMX      *p7_omx_Create(int allocM);
extern int          p7_omx_GrowTo(P7_OMX *ox, int allocM);
extern void         p7_omx_Destroy(P7_OMX *ox);

extern int          p7_oprofile_Dump(FILE *fp, P7_OPROFILE *om);
extern int          p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox, int truefalse);

extern int          p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_ReconfigLength(P7_OPROFILE *om, int L);

extern int p7_MSVFilter    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_ForwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_ViterbiScore (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);

#endif /*P7_IMPL_DUMMY_INCLUDED*/
#endif /*p7_IMPL_DUMMY*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
