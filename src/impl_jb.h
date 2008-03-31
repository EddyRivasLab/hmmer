/* Optimized implementation of HMMER2 from Jeremy Buhler (Washington
 * University, St. Louis); adapted to HMMER3.
 * 
 * SRE, Wed Jun 27 07:59:50 2007 [Janelia]
 * SVN $Id$
 */
#ifndef P7_IMPL_JB_INCLUDED
#define P7_IMPL_JB_INCLUDED

/*****************************************************************
 * 1. P7_OPROFILE: a scoring profile
 *****************************************************************/

/* Indices for special state types in the length model, om->xsc[x][]
 */
enum p7o_xstates_e {
  p7O_E = 0,
  p7O_N = 1, 
  p7O_J = 2, 
  p7O_C = 3 
};
#define p7O_NXSTATES  4

/* Indices for transitions from the length modeling scores om->xsc[][x]
 */
enum p7o_xtransitions_e { 
    p7O_LOOP = 0, 
    p7O_MOVE = 1 
};
#define p7O_NXTRANS  2

/* Indices for transition scores gm->tsc[k][] */
/* order is optimized for dynamic programming */
enum p7o_tsc_e { 
  p7O_MM = 0,
  p7O_IM = 1, 
  p7O_DM = 2, 
  p7O_MD = 3,
  p7O_DD = 4, 
  p7O_MI = 5, 
  p7O_II = 6, 
  p7O_BM = 7 
};
#define p7O_NTRANS 8

/* Indices for residue emission score vectors
 */
enum p7o_rsc_e {
  p7O_MSC = 0,
  p7O_ISC = 1 
};
#define p7O_NR 2

typedef struct p7_oprofile_s {
  int     mode;
  int     M;
  int    *tsc;	/* [0.1..M-1][0..p7X_NTSC-1] */
  int   **rsc;	/* [0..Kp-1][0.1..M][p7X_NR] */
  int     xsc[p7O_NXSTATES][p7O_NXTRANS];

  const ESL_ALPHABET *abc;
} P7_OPROFILE;


/*****************************************************************
 * 2. P7_OMX: a dynamic programming matrix
 *****************************************************************/


enum p7x_scells_e {
  p7X_M = 0, 
  p7X_I = 1,
  p7X_D = 2 
};
#define p7X_NSCELLS 3

enum p7x_xcells_e { 
  p7X_E = 0, 
  p7X_N = 1, 
  p7X_J = 2,
  p7X_B = 3, 
  p7X_C = 4 
};
#define p7X_NXCELLS 5

typedef struct p7_omx_s {
  int M;
  int L;

  size_t ncells;
  size_t nrows;
  
  int **dp;			/* [0.1..L][0.1..M][0..p7X_NSCELLS-1] */
  int  *xmx;              	/* [0.1..L][0..p7X_NXCELLS-1] */
} P7_OMX;


extern P7_OPROFILE *p7_oprofile_Create(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy(P7_OPROFILE *om);
extern P7_OMX      *p7_omx_Create(int allocM, int allocL);
extern void         p7_omx_Destroy(P7_OMX *ox);
extern int          p7_omx_GrowTo(P7_OMX *ox, int allocM, int allocL);
extern int          p7_omx_Dump(FILE *ofp, P7_OMX *ox);
extern int          p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_Viterbi(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int          p7_Forward(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);



#endif /*P7_IMPL_JB_INCLUDED*/
