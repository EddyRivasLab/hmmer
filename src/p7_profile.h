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
#define p7_LOCAL     1		/* multi-hit local:  "fs" mode   */
#define p7_GLOCAL    2		/* multi-hit glocal: "ls" mode   */
#define p7_UNILOCAL  3		/* one-hit local: "sw" mode      */
#define p7_UNIGLOCAL 4		/* one-hit glocal: "s" mode      */


/* The P7_PROFILE structure: generic score profile.
 */
typedef struct {
  int    mode;         	/* configured algorithm mode (e.g. p7_LOCAL)   */ 
  ESL_ALPHABET *abc;	/* copy of pointer to appropriate alphabet     */
  int    M;		/* number of nodes in the model                */
  int  **tsc;           /* transition scores     [0.6][1.M-1]          */
  int  **msc;           /* match emission scores [0.Kp-1][1.M]         */
  int  **isc;           /* ins emission scores   [0.Kp-1][1.M-1]       */
  int    xsc[4][2];     /* N,E,C,J transitions   [][LOOP,MOVE]         */
  int   *bsc;           /* begin transitions     [1.M]                 */
  int   *esc;		/* end transitions       [1.M]                 */
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
extern int         p7_profile_GetTScore(P7_PROFILE *gm, 
					char st1, int k1, char st2, int k2,
					int *ret_tsc);


#endif /*P7_PROFILE_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/


