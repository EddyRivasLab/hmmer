/* e1_rate
 *
 *   
*/
#ifndef E2_PROFILE_INCLUDED
#define E2_PROFILE_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"
#include "esl_alphabet.h"	/* ESL_DSQ, ESL_ALPHABET */
#include "esl_dmatrix.h"	/* ESL_DMATRIX           */
#include "esl_getopts.h"	/* ESL_GETOPTS           */
#include "esl_tree.h"

#include "e2_config.h"
#include "e2.h"

/*******************************************************************************
 * 3. E2_PROFILE: a scoring profile for two descendant from a common ancestral
 *******************************************************************************/
/* Indices for six special state types x that have transition parameters
 * in gm->xsc[x][y], where all of them have two choices y.
 */
#define e2P_NXSTATES 7
#define e2P_NXTRANS  2
#define e2P_LOOP 0              /* gm->xsc[x][e2P_LOOP]    gm->xsc[x][e2P_MOVE] */        
#define e2P_MOVE 1              /* -------------------     -------------------- */
#define e2P_EE  0            	/*     EE->J1                     EE->C1        */
#define e2P_N1  1		/*     N1->N1                     N1->N2        */
#define e2P_N2  2		/*     N2->N2                     N2->BB        */
#define e2P_J1  3               /*     J1->J1                     J1->J2        */
#define e2P_J2  4               /*     J2->J2                     J2->BB        */
#define e2P_C1  5		/*     C1->C1                     C1->C2        */
#define e2P_C2  6		/*     C2->C2                     C2->T         */

/* Indices for state types
 */
enum e2p_states_e { 
  e2P_BB = 0,
  e2P_IB = 1,
  e2P_SS = 2,
  e2P_DS = 3,
  e2P_IS = 4,
  e2P_SD = 5,
  e2P_DD = 6,
  e2P_ID = 7,
  e2P_BI = 8,
  e2P_SI = 9,
  e2P_DI = 10,
  e2P_II = 11
};
#define e2P_NSTATES 12

/* Indices for transition scores gm->tsc[] */
/* order is optimized for dynamic programming */
enum e2p_tsc_e {
  e2P_BB_IB = 0,  // 2 to IB
  e2P_IB_IB = 1,
  e2P_BB_SS = 2,  // 12 to SS
  e2P_IB_SS = 3,
  e2P_SS_SS = 4,
  e2P_DS_SS = 5,
  e2P_IS_SS = 6,
  e2P_SD_SS = 7,
  e2P_DD_SS = 8,
  e2P_ID_SS = 9,
  e2P_BI_SS = 10,
  e2P_SI_SS = 11,
  e2P_DI_SS = 12,
  e2P_II_SS = 13,
  e2P_BB_DS = 14,  // 12 to DS
  e2P_IB_DS = 15,
  e2P_SS_DS = 16,
  e2P_DS_DS = 17,
  e2P_IS_DS = 18,
  e2P_SD_DS = 19,
  e2P_DD_DS = 20,
  e2P_ID_DS = 21,
  e2P_BI_DS = 22,
  e2P_SI_DS = 23,
  e2P_DI_DS = 24,
  e2P_II_DS = 25,
  e2P_SS_IS = 26,  // 3 to IS
  e2P_DS_IS = 27,  
  e2P_IS_IS = 28, 
  e2P_BB_SD = 29,  // 12 to SD
  e2P_IB_SD = 30,
  e2P_SS_SD = 31,
  e2P_DS_SD = 32,
  e2P_IS_SD = 33,
  e2P_SD_SD = 34,
  e2P_DD_SD = 35,
  e2P_ID_SD = 36,
  e2P_BI_SD = 37,
  e2P_SI_SD = 38,
  e2P_DI_SD = 39,
  e2P_II_SD = 40,
  e2P_IB_DD = 41, // 11 (but 10 used) to DD but DD_DD is used to normalize 
  e2P_SS_DD = 42,
  e2P_DS_DD = 43,
  e2P_IS_DD = 44,
  e2P_SD_DD = 45,
  e2P_DD_DD = 46,
  e2P_ID_DD = 47,
  e2P_BI_DD = 48,
  e2P_SI_DD = 49,
  e2P_DI_DD = 50,
  e2P_II_DD = 51,
  e2P_SD_ID = 52,  // 3 to ID
  e2P_DD_ID = 53,  
  e2P_ID_ID = 54, 
  e2P_BB_BI = 55,  // 2 to BI
  e2P_BI_BI = 56,
  e2P_SS_SI = 57,  // 3 to SI
  e2P_SD_SI = 58,  
  e2P_SI_SI = 59, 
  e2P_DS_DI = 60,  // 3 to DI
  e2P_DD_DI = 61,  
  e2P_DI_DI = 62,  
  e2P_IB_II = 63,  // 4 to II
  e2P_IS_II = 64,  
  e2P_ID_II = 65,  
  e2P_II_II = 66, 
  e2P_IB_EE = 67,  // 11 to E 
  e2P_SS_EE = 68,  
  e2P_DS_EE = 69,  
  e2P_IS_EE = 70,
  e2P_SD_EE = 71,
  e2P_DD_EE = 72,
  e2P_ID_EE = 73,
  e2P_BI_EE = 74,
  e2P_SI_EE = 75,
  e2P_DI_EE = 76,
  e2P_II_EE = 77
};
#define e2P_NTRANS 78

/* Indices for sequences
 */
enum e2p_sq_e {
  e2P_SL = 0, 
  e2P_SR = 1
};
#define e2P_NS 2

/* Accessing transition, emission scores */
#define e2P_TSC(gm, s)     ((gm)->tsc[(s)])
#define e2P_XSCLOOP(gm, s) ((gm)->xsc[(s)][e2P_LOOP])
#define e2P_XSCMOVE(gm, s) ((gm)->xsc[(s)][e2P_MOVE])
#define e2P_SSSC(gm, x, y) ((gm)->sssc[(x)][(y)])
#define e2P_SLSC(gm, x)    ((gm)->ssc[(x)][e2P_SL])
#define e2P_SRSC(gm, y)    ((gm)->ssc[(y)][e2P_SR])
#define e2P_ILSC(gm, x)    ((gm)->isc[(x)][e2P_SL])
#define e2P_IRSC(gm, y)    ((gm)->isc[(y)][e2P_SR])
#define e2P_FLSC(gm, x)    ((gm)->fsc[(x)][e2P_SL])
#define e2P_FRSC(gm, y)    ((gm)->fsc[(y)][e2P_SR])

typedef struct e2_profile_s {
  float  tsc[e2P_NTRANS];                 /* transitions  [0..e2P_NTRANS-1], hand-indexed                                */
  float  sssc[e2_MAXABET][e2_MAXABET];    /* seq1/seq2 joint subtitution emissions [0..K-1][0..K-1], hand-indexed        */
  float  ssc[e2_MAXABET][e2P_NS];         /* orphan subtitution emissions [0..K-1][0..e2P_NR-1], hand-indexed            */
  float  isc[e2_MAXABET][e2P_NS];         /* insertion emissions [0..K-1][0..e2P_NR-1], hand-indexed                     */
  float  xsc[e2P_NXSTATES][e2P_NXTRANS];  /* special transitions [ENJCBG][LOOP,MOVE] */
  float  fsc[e2_MAXABET][e2P_NS];         /* flanking emissions [0..K-1][0..e2P_NR-1], hand-indexed                      */

  int    mode;                            /* configured algorithm mode (e.g. e2_LOCAL) */ 

  /* Info, most of which is a copy from parent HMM:                                       */
  char  *name;			/* unique name of model                                   */
  char  *acc;			/* unique accession of model, or NULL                     */
  char  *desc;                  /* brief (1-line) description of model, or NULL           */

  const ESL_ALPHABET *abc;	/* copy of pointer to appropriate alphabet                */
} E2_PROFILE;

/* e2_profile.c */
extern E2_PROFILE *e2_profile_Create(const ESL_ALPHABET *abc);
extern int         e2_profile_Copy(const E2_PROFILE *src, E2_PROFILE *dst);
extern E2_PROFILE *e2_profile_Clone(const E2_PROFILE *gm);
extern int         e2_profile_Reuse(E2_PROFILE *gm);
extern size_t      e2_profile_Sizeof(E2_PROFILE *gm);
extern void        e2_profile_Destroy(E2_PROFILE *gm);
extern int         e2_profile_GetT(const E2_PROFILE *gm, char st1, char st2, float *ret_tsc);

#endif
