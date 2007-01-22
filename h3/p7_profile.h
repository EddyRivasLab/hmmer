/* The P7_PROFILE structure - a search profile.
 *                            
 * SRE, Thu Jan 11 15:05:51 2007 [Janelia] [Holst, The Planets]
 * SVN $Id$
 */
#ifndef P7_PROFILEH_INCLUDED
#define P7_PROFILEH_INCLUDED

#include "esl_alphabet.h"

/*----------------------------------------------------------------
 * P7_PROFILE 
 * A Plan 7 search profile.
 */
typedef struct {
  int    M;		/* number of nodes in the model                */
  int  **tsc;           /* transition scores     [0.6][1.M-1]          */
  int  **msc;           /* match emission scores [0.Kp-1][1.M]         */
  int  **isc;           /* ins emission scores   [0.Kp-1][1.M-1]       */
  int    xsc[4][2];     /* N,E,C,J transitions   [][LOOP,MOVE]         */
  int   *bsc;           /* begin transitions     [1.M]                 */
  int   *esc;		/* end transitions       [1.M]                 */

  /* We also have some probabilities relevant to the search profile but
   * not to the core model.
   */
  float *xt[4][2];	/* [NECJ][MOVE,LOOP] transitions               */
  float *begin;		/* 1..M begin "probabilities"                  */
  float *end;		/* 1..M end "probabilities"                    */

  /* Objects we keep references to
   */
  ESL_ALPHABET *abc;	/* copy of pointer to appropriate alphabet     */
  P7_HMM       *hmm;	/* who's your daddy                            */
  P7_NULL      *bg;	/* background null model                       */
  
  int    mode;         	/* configured algorithm mode (e.g. p7_LOCAL)   */ 
} P7_PROFILE;

/* Indices for special state types, I: used for dynamic programming xmx[][]
 * mnemonic: eXtra Matrix for B state = XMB
 */
#define p7_XMB 0
#define p7_XME 1
#define p7_XMC 2
#define p7_XMJ 3
#define p7_XMN 4
p
/* Indices for special state types, II: used for hmm->xt[] indexing
 * mnemonic: eXtra Transition for N state = XTN
 */
#define p7_XTN  0
#define p7_XTE  1
#define p7_XTC  2
#define p7_XTJ  3

/* Indices for extra state transitions
 * Used for indexing hmm->xt[][].
 */
#define p7_MOVE 0          /* trNB, trEC, trCT, trJB */
#define p7_LOOP 1          /* trNN, trEJ, trCC, trJJ */

/* Search modes.
 */
#define p7_LOCAL     1		/* multi-hit local:  "fs" mode   */
#define p7_GLOCAL    2		/* multi-hit glocal: "ls" mode   */
#define p7_UNILOCAL  3		/* one-hit local: "sw" mode      */
#define p7_UNIGLOCAL 4		/* one-hit glocal: "s" mode      */


#endif /*P7_PROFILEH_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/ 
