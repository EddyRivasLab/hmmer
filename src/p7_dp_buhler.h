/* p7_dp_buhler.h
 * 
 * Definitions of optimized profile and DP matrix (P7_OPROFILE,
 * P7_OMX) structures for Jeremy Buhler's DP implementation.
 *
 * For information, please contact jbuhler@cse.wustl.edu.
 * Copyright (C) 2005 Washington University School of Medicine
 *
 * SVN $Id$
 */

#ifndef P7_DPBUHLER_INCLUDED
#define P7_DPBUHLER_INCLUDED

/* Model constants for JB implementation of Viterbi algo */

#define N_XPARMS 4
enum xParms { _XTE = 0, _XTN = 1, _XTJ = 2, _XTC = 3 };

#define N_XMOVES 2
enum xMoves { _LOOP = 0, _MOVE = 1 };

#define N_TMOVES 9
enum tMoves { _TMM = 0, _TIM = 1, _TDM = 2, 
              _TMD = 3, _TDD = 4, 
              _TMI = 5, _TII = 6, 
              _BSC = 7, _ESC = 8 };

#define N_MISC 2
enum miscTypes { _MSC = 0, _ISC = 1 };

typedef struct {
  int *tsc;    /*                         (0 .. M) x N_TMOVES */
  int **misc;  /* (0..AlphabetSize - 1) x (1 .. M) x N_MISC */
  
  int xsc[N_XPARMS][N_XMOVES];
  int M;
} P7_OPROFILE;


/* DP matrix constants for JB implementation of Viterbi algo */

#define N_MSTATES 3
enum mStates { _MMX = 0, _IMX = 1, _DMX = 2 };

#define N_XSTATES 5
enum xStates { _XME = 0, _XMN = 1, 
               _XMJ = 2, _XMB = 3, 
               _XMC = 4 };

typedef struct {
  int **dp;    /* (0 .. SeqLength) x (0 .. M) x N_MSTATES */
  int *xmx;    /* (0 .. SeqLength) x            N_XSTATES */
  
  unsigned int maxM, maxN;
  unsigned int padM, padN;
} P7_OMX;

#endif /*P7_DPBUHLER_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
