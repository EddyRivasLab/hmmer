/* p7_profile.h
 * The "common representation" of a HMMER score profiles and DP matrices,
 * P7_PROFILE and P7_GMX.
 * 
 * SVN $Id$
 * SRE, Tue Mar 21 14:54:43 2006 [St. Louis]
 */
#ifndef P7_PROFILE_INCLUDED
#define P7_PROFILE_INCLUDED


/* Search modes.
 */
enum p7_searchmode {
  P7_NO_MODE    = 0,
  P7_LS_MODE    = 1,
  P7_FS_MODE    = 2,
  P7_SW_MODE    = 3,
  P7_S_MODE     = 4,
};

/* The P7_PROFILE structure: generic score profile.
 */
typedef struct {
  enum p7_searchmode mode;	/* configured algorithm mode                */
  ESL_ALPHABET *abc;		/* copy of pointer to appropriate alphabet  */
  int    M;
  int  **tsc;                   /* transition scores     [0.6][1.M-1]       */
  int  **msc;                   /* match emission scores [0.Kp-1][1.M]      */
  int  **isc;                   /* ins emission scores   [0.Kp-1][1.M-1]    */
  int    xsc[4][2];             /* N,E,C,J transitions   [][LOOP,MOVE]      */
  int   *bsc;                   /* begin transitions     [1.M]              */
  int   *esc;			/* end transitions       [1.M]              */
} P7_PROFILE;


/* The P7_GMX structure: the generic DP matrix.
 */
typedef struct {
  int **xmx;			/* special scores [0.1..N][BECJN]     */
  int **mmx;			/* match scores [0.1..N][0.1..M]      */
  int **imx;			/* insert scores [0.1..N][0.1..M-1.M] */
  int **dmx;			/* delete scores [0.1..N][0.1..M-1.M] */
  
  int N;			/* alloc'ed for seq of length N; N+1 rows */
  int M;			/* alloc'ed for HMM of length M; M+1 cols */

  /* If either pad is 0, we're not growable in that direction.       */
  int padN;			/* extra pad in sequence length/rows */
  int padM;			/* extra pad in HMM length/columns   */
} P7_GMX;


extern P7_PROFILE *p7_profile_Create(int M);
extern void        p7_profile_Destroy(P7_PROFILE *gm);

#endif /*P7_PROFILE_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/


