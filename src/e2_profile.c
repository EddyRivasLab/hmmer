/* Routines for the E2_PROFILE structure - e2 profile with two evolutionary models
 *                                         
 *    1. The E2_PROFILE object: allocation, initialization, destruction.
 *    2. Access methods.
 *    3. Debugging and development code.
 *    4. Unit tests.
 *    5. Test driver.
 *
 * See also: 
 *   modelconfig.c : routines that configure a profile given two evolutionary models
 */
#include "p7_config.h"
#include "hmmer.h"

#include "esl_vectorops.h"

#include "e1_model.h"
#include "e2.h"
#include "e2_profile.h"
#include "e2hmmer_profile.h"

/*****************************************************************
 * 1. The E2_PROFILE object: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  e2_profile_Create()
 * Synopsis:  Allocates a profile.
 *
 * Purpose:   Allocates for a profile of up to <M> nodes, for digital
 *            alphabet <abc>.
 *            
 *            Because this function might be in the critical path (in
 *            hmmscan, for example), we leave much of the model
 *            unintialized, including scores and length model
 *            probabilities. The <e2_ProfileConfig()> call is what
 *            sets these. 
 *            
 *            The alignment mode is set to <e2_NO_MODE>.  The
 *            reference pointer <gm->abc> is set to <abc>.
 *
 * Returns:   a pointer to the new profile.
 *
 * Throws:    <NULL> on allocation error.
 *
 * Xref:      STL11/125.
 */
E2_PROFILE *
e2_profile_Create(const ESL_ALPHABET *abc)
{
  E2_PROFILE *gm = NULL;
  int         status;

  /* level 0 */
  ESL_ALLOC(gm, sizeof(E2_PROFILE));

  /* Set remaining info  */
  gm->mode             = e2_NOMODE;
  gm->name             = NULL;
  gm->acc              = NULL;
  gm->desc             = NULL;  
  gm->abc              = abc;
  return gm;

 ERROR:
  e2_profile_Destroy(gm);
  return NULL;
}


/* Function:  e2_profile_Copy()
 * Synopsis:  Copy a profile.
 *
 * Purpose:   Copies profile <src> to profile <dst>, where <dst>
 *            has already been allocated to be of sufficient size.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if <dst> is too small 
 *            to fit <src>.
 */
int
e2_profile_Copy(const E2_PROFILE *src, E2_PROFILE *dst)
{
  int x;
  int status;

  esl_vec_FCopy(src->tsc, e2P_NTRANS, dst->tsc);
  for (x = 0; x < src->abc->K;   x++) esl_vec_FCopy(src->sssc[x], e2_MAXABET, dst->sssc[x]);
  for (x = 0; x < src->abc->K;   x++) esl_vec_FCopy(src->ssc[x],  e2P_NS,     dst->ssc[x]);
  for (x = 0; x < src->abc->K;   x++) esl_vec_FCopy(src->isc[x],  e2P_NS,     dst->isc[x]);
  for (x = 0; x < src->abc->K;   x++) esl_vec_FCopy(src->fsc[x],  e2P_NS,     dst->fsc[x]);
  
  dst->mode = src->mode;

  if (dst->name != NULL) free(dst->name);
  if (dst->acc  != NULL) free(dst->acc);
  if (dst->desc != NULL) free(dst->desc);

  if ((status = esl_strdup(src->name,      -1, &(dst->name)))      != eslOK) return status;
  if ((status = esl_strdup(src->acc,       -1, &(dst->acc)))       != eslOK) return status;
  if ((status = esl_strdup(src->desc,      -1, &(dst->desc)))      != eslOK) return status;

  return eslOK;
}

/* Function:  e2_profile_Clone()
 * Synopsis:  Duplicates a profile.
 *
 * Purpose:   Duplicate profile <gm>; return a pointer
 *            to the newly allocated copy.
 */
E2_PROFILE *
e2_profile_Clone(const E2_PROFILE *gm)
{
  E2_PROFILE *g2 = NULL;
  int         status;

  if ((g2 = e2_profile_Create(gm->abc)) == NULL) return NULL;
  if ((status = e2_profile_Copy(gm, g2)) != eslOK) goto ERROR;
  return g2;
  
 ERROR:
  e2_profile_Destroy(g2);
  return NULL;
}


/* Function:  e2_profile_Reuse()
 * Synopsis:  Prepare profile to be re-used for a new HMM.
 *
 * Purpose:   Prepare profile <gm>'s memory to be re-used
 *            for a new HMM.
 */
int
e2_profile_Reuse(E2_PROFILE *gm)
{
  /* name, acc, desc annotation is dynamically allocated for each HMM */
  if (gm->name != NULL) { free(gm->name); gm->name = NULL; }
  if (gm->acc  != NULL) { free(gm->acc);  gm->acc  = NULL; }
  if (gm->desc != NULL) { free(gm->desc); gm->desc = NULL; }

  /* reset some other things, but leave the rest alone. */
  gm->mode = e2_NOMODE;

  return eslOK;
}


/* Function:  e2_profile_Destroy()
 * Synopsis:  Frees a profile.
 *
 * Purpose:   Frees a profile <gm>.
 *
 * Returns:   (void).
 *
 * Xref:      STL11/125.
 */
void
e2_profile_Destroy(E2_PROFILE *gm)
{
  if (gm != NULL) {
    if (gm->name      != NULL) free(gm->name);
    if (gm->acc       != NULL) free(gm->acc);
    if (gm->desc      != NULL) free(gm->desc);
    free(gm);
  }
  return;
}

/*****************************************************************
 * 2. Access methods.
 *****************************************************************/
 
int
e2_profile_GetT(const E2_PROFILE *gm, char st1, char st2, float *ret_tsc)
{
  float tsc = 0.0f;
  int   status;

  switch (st1) {
  case e2T_S:   break;
  case e2T_T:   break;
  case e2T_EE:  
    switch(st2) {
    case e2T_J1:   tsc = e2P_XSCLOOP(gm, e2P_EE); break;
    case e2T_C1:   tsc = e2P_XSCMOVE(gm, e2P_EE); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
  case e2T_N1:  
    switch(st2) {
    case e2T_N1:   tsc = e2P_XSCLOOP(gm, e2P_N1); break;
    case e2T_N2:   tsc = e2P_XSCMOVE(gm, e2P_N1); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
 case e2T_N2:  
    switch(st2) {
    case e2T_N2:  tsc = e2P_XSCLOOP(gm, e2P_N2); break;
    case e2T_BB:  tsc = e2P_XSCMOVE(gm, e2P_N2); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;

 case e2T_J1:  
    switch(st2) {
    case e2T_J1:   tsc = e2P_XSCLOOP(gm, e2P_J1); break;
    case e2T_J2:   tsc = e2P_XSCMOVE(gm, e2P_J1); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
 case e2T_J2:  
    switch(st2) {
    case e2T_J2:  tsc = e2P_XSCLOOP(gm, e2P_J2); break;
    case e2T_BB:  tsc = e2P_XSCMOVE(gm, e2P_J2); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;

 case e2T_C1:  
    switch(st2) {
    case e2T_C1:  tsc = e2P_XSCLOOP(gm, e2P_C1); break;
    case e2T_C2:  tsc = e2P_XSCMOVE(gm, e2P_C1); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
 case e2T_C2:  
    switch(st2) {
    case e2T_C2:  tsc = e2P_XSCLOOP(gm, e2P_C2); break;
    case e2T_T:   tsc = e2P_XSCMOVE(gm, e2P_C2); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;

  case e2T_BB:  
    switch (st2) {
    case e2T_IB:  tsc = e2P_TSC(gm, e2P_BB_IB); break;
    case e2T_SS:  tsc = e2P_TSC(gm, e2P_BB_SS); break;
    case e2T_DS:  tsc = e2P_TSC(gm, e2P_BB_DS); break;
    case e2T_SD:  tsc = e2P_TSC(gm, e2P_BB_SD); break;
    case e2T_BI:  tsc = e2P_TSC(gm, e2P_BB_BI); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
  case e2T_IB:  
    switch (st2) {
    case e2T_BB:  break;
    case e2T_IB:  tsc = e2P_TSC(gm, e2P_IB_IB); break;
    case e2T_SS:  tsc = e2P_TSC(gm, e2P_IB_SS); break;
    case e2T_DS:  tsc = e2P_TSC(gm, e2P_IB_DS); break;
    case e2T_SD:  tsc = e2P_TSC(gm, e2P_IB_SD); break;
    case e2T_DD:  tsc = e2P_TSC(gm, e2P_IB_DD); break;
    case e2T_II:  tsc = e2P_TSC(gm, e2P_IB_II); break;
    case e2T_EE:  tsc = e2P_TSC(gm, e2P_IB_EE); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
  case e2T_SS:  
    switch (st2) {
    case e2T_SS:  tsc = e2P_TSC(gm, e2P_SS_SS); break;
    case e2T_DS:  tsc = e2P_TSC(gm, e2P_SS_DS); break;
    case e2T_IS:  tsc = e2P_TSC(gm, e2P_SS_IS); break;
    case e2T_SD:  tsc = e2P_TSC(gm, e2P_SS_SD); break;
    case e2T_DD:  tsc = e2P_TSC(gm, e2P_SS_DD); break;
    case e2T_SI:  tsc = e2P_TSC(gm, e2P_SS_SI); break;
    case e2T_EE:  tsc = e2P_TSC(gm, e2P_SS_EE); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
  case e2T_DS:  
    switch (st2) {
    case e2T_SS:  tsc = e2P_TSC(gm, e2P_DS_SS); break;
    case e2T_DS:  tsc = e2P_TSC(gm, e2P_DS_DS); break;
    case e2T_IS:  tsc = e2P_TSC(gm, e2P_DS_IS); break;
    case e2T_SD:  tsc = e2P_TSC(gm, e2P_DS_SD); break;
    case e2T_DD:  tsc = e2P_TSC(gm, e2P_DS_DD); break;
    case e2T_DI:  tsc = e2P_TSC(gm, e2P_DS_DI); break;
    case e2T_EE:  tsc = e2P_TSC(gm, e2P_DS_EE); break;
     default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
  case e2T_IS:  
    switch (st2) {
    case e2T_SS:  tsc = e2P_TSC(gm, e2P_IS_SS); break;
    case e2T_DS:  tsc = e2P_TSC(gm, e2P_IS_DS); break;
    case e2T_IS:  tsc = e2P_TSC(gm, e2P_IS_IS); break;
    case e2T_SD:  tsc = e2P_TSC(gm, e2P_IS_SD); break;
    case e2T_DD:  tsc = e2P_TSC(gm, e2P_IS_DD); break;
    case e2T_II:  tsc = e2P_TSC(gm, e2P_IS_II); break;
    case e2T_EE:  tsc = e2P_TSC(gm, e2P_IS_EE); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
  case e2T_SD:  
    switch (st2) {
    case e2T_SS:  tsc = e2P_TSC(gm, e2P_SD_SS); break;
    case e2T_DS:  tsc = e2P_TSC(gm, e2P_SD_DS); break;
    case e2T_SD:  tsc = e2P_TSC(gm, e2P_SD_SD); break;
    case e2T_DD:  tsc = e2P_TSC(gm, e2P_SD_DD); break;
    case e2T_ID:  tsc = e2P_TSC(gm, e2P_SD_ID); break;
    case e2T_SI:  tsc = e2P_TSC(gm, e2P_SD_SI); break;
    case e2T_EE:  tsc = e2P_TSC(gm, e2P_SD_EE); break;
     default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
  case e2T_DD:  
    switch (st2) {
    case e2T_SS:  tsc = e2P_TSC(gm, e2P_DD_SS); break;
    case e2T_DS:  tsc = e2P_TSC(gm, e2P_DD_DS); break;
    case e2T_SD:  tsc = e2P_TSC(gm, e2P_DD_DS); break;
    case e2T_DD:  tsc = e2P_TSC(gm, e2P_DD_DD); break;
    case e2T_ID:  tsc = e2P_TSC(gm, e2P_DD_ID); break;
    case e2T_DI:  tsc = e2P_TSC(gm, e2P_DD_DI); break;
    case e2T_EE:  tsc = e2P_TSC(gm, e2P_DD_EE); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
  case e2T_ID:  
    switch (st2) {
    case e2T_SS:  tsc = e2P_TSC(gm, e2P_ID_SS); break;
    case e2T_DS:  tsc = e2P_TSC(gm, e2P_ID_DS); break;
     case e2T_SD:  tsc = e2P_TSC(gm, e2P_ID_SD); break;
    case e2T_DD:  tsc = e2P_TSC(gm, e2P_ID_DD); break;
    case e2T_ID:  tsc = e2P_TSC(gm, e2P_ID_ID); break;
    case e2T_II:  tsc = e2P_TSC(gm, e2P_ID_II); break;
    case e2T_EE:  tsc = e2P_TSC(gm, e2P_ID_EE); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
  case e2T_BI:  
    switch (st2) {
    case e2T_SS:  tsc = e2P_TSC(gm, e2P_BI_SS); break;
    case e2T_DS:  tsc = e2P_TSC(gm, e2P_BI_DS); break;
    case e2T_SD:  tsc = e2P_TSC(gm, e2P_BI_DS); break;
    case e2T_DD:  tsc = e2P_TSC(gm, e2P_BI_EE); break;
    case e2T_BI:  tsc = e2P_TSC(gm, e2P_BI_DD); break;
    case e2T_EE:  tsc = e2P_TSC(gm, e2P_BI_EE); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
  case e2T_SI:  
    switch (st2) {
    case e2T_SS:  tsc = e2P_TSC(gm, e2P_SI_SS); break;
    case e2T_DS:  tsc = e2P_TSC(gm, e2P_SI_DS); break;
    case e2T_SD:  tsc = e2P_TSC(gm, e2P_SI_SD); break;
    case e2T_DD:  tsc = e2P_TSC(gm, e2P_SI_DD); break;
    case e2T_SI:  tsc = e2P_TSC(gm, e2P_SI_SI); break;
    case e2T_EE:  tsc = e2P_TSC(gm, e2P_SI_EE); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;
  case e2T_DI:  
    switch (st2) {
    case e2T_SS:  tsc = e2P_TSC(gm, e2P_DI_SS); break;
    case e2T_DS:  tsc = e2P_TSC(gm, e2P_DI_DS); break;
    case e2T_SD:  tsc = e2P_TSC(gm, e2P_DI_SD); break;
    case e2T_DD:  tsc = e2P_TSC(gm, e2P_DI_DD); break;
    case e2T_DI:  tsc = e2P_TSC(gm, e2P_DI_DI); break;
    case e2T_EE:  tsc = e2P_TSC(gm, e2P_DI_EE); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
   break;
  case e2T_II:  
    switch (st2) {
    case e2T_SS:  tsc = e2P_TSC(gm, e2P_II_SS); break;
    case e2T_DS:  tsc = e2P_TSC(gm, e2P_II_DS); break;
    case e2T_SD:  tsc = e2P_TSC(gm, e2P_II_SD); break;
    case e2T_DD:  tsc = e2P_TSC(gm, e2P_II_DD); break;
    case e2T_II:  tsc = e2P_TSC(gm, e2P_II_II); break;
    case e2T_EE:  tsc = e2P_TSC(gm, e2P_II_EE); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", e2_model_DecodeStatetype(st1), e2_model_DecodeStatetype(st2));
    }
    break;

  default: ESL_XEXCEPTION(eslEINVAL, "bad state type %d in traceback", st1);
  }
 
  *ret_tsc = tsc;
  return eslOK;

 ERROR:
  *ret_tsc = -eslINFINITY;
  return status;
}

