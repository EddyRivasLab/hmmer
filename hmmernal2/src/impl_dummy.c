/* The "dummy" implementation provides an "optimized" interface when
 * none is available, for portability to all ANSI C99 platforms.
 * 
 * SRE, Sun Feb 17 10:10:45 2008 [Janelia]
 * SVN $Id$
 */

#include "p7_config.h"
#ifdef p7_IMPL_DUMMY
#include "hmmer.h"
#include "easel.h"

P7_OPROFILE *
p7_oprofile_Create(int allocM, const ESL_ALPHABET *abc)
{
  return (P7_OPROFILE *) p7_profile_Create(allocM, abc);
}

void
p7_oprofile_Destroy(P7_OPROFILE *om)
{
  p7_profile_Destroy((P7_PROFILE *) om);
}

P7_OMX *
p7_omx_Create(int allocM)
{
  return (P7_OMX *) p7_gmx_Create(allocM, 1); /* L=1 makes it a two-row matrix */
}

int
p7_omx_GrowTo(P7_OMX *ox, int allocM)
{
  return p7_gmx_GrowTo((P7_GMX *) ox, allocM, 1);
}

void
p7_omx_Destroy(P7_OMX *ox)
{
  p7_gmx_Destroy((P7_GMX *) ox);
}

int
p7_oprofile_Dump(FILE *fp, P7_OPROFILE *om)
{
  return eslOK; 		/* FIXME: currently a no-op! */
}

int
p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox, int truefalse)
{
  return eslOK;			/* FIXME: currently a no-op! */
}

int
p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om)
{
  return p7_profile_Copy(gm, (P7_PROFILE *) om);
}

int
p7_oprofile_ReconfigLength(P7_OPROFILE *om, int L)
{
  return p7_ReconfigLength((P7_PROFILE *) om, L);
}

int
p7_MSVFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  return p7_GMSV(dsq, L, (const P7_PROFILE *) om, (P7_GMX *) ox, ret_sc);
}


int
p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  return p7_GViterbi(dsq, L, (const P7_PROFILE *) om, (P7_GMX *) ox, ret_sc);
}

int 
p7_ForwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  return p7_GForward(dsq, L, (const P7_PROFILE *) om, (P7_GMX *) ox, ret_sc);
}

int
p7_ViterbiScore(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  return p7_GViterbi(dsq, L, (const P7_PROFILE *) om, (P7_GMX *) ox, ret_sc);
}




#endif /* p7_IMPL_DUMMY */

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
