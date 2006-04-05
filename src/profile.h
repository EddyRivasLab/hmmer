/* p7_profile.h
 * The "common representation" of a HMMER score profile. 
 * Usually built from a profile HMM's probabilities.
 * 
 * SRE, Tue Mar 21 14:54:43 2006 [St. Louis]
 * SVN $Id$
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

/* The P7_PROFILE structure
 */
typedef struct {
  enum p7_searchmode mode;	/* configured algorithm mode                */
  int    M;
  int  **tsc;                   /* transition scores     [0.6][1.M-1]       */
  int  **msc;                   /* match emission scores [0.MAXCODE-1][1.M] */
  int  **isc;                   /* ins emission scores [0.MAXCODE-1][1.M-1] */
  int    xsc[4][2];             /* N,E,C,J transitions   [][LOOP,MOVE]      */
  int   *bsc;                   /* begin transitions     [1.M]              */
  int   *esc;			/* end transitions       [1.M]              */
} P7_PROFILE;



extern P7_PROFILE *p7_profile_Create(int M);
extern void        p7_profile_Destroy(P7_PROFILE *gm);

#endif /*P7_PROFILE_INCLUDED*/
