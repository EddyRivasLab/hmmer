/* Original HMMER2 implementation of dynamic programming algorithms.
 * 
 * This may be useful for reference or for regression testing.
 * 
 * Implementation plugins provide a standard API, including support
 * for two objects, a P7_OPROFILE score profile and a P7_OMX dynamic
 * programming matrix.
 * 
 * SRE, Fri Jul 13 07:51:24 2007
 * SVN $Id$
 */

#ifndef P7_IMPL_H2_INCLUDED
#define P7_IMPL_H2_INCLUDED

/*****************************************************************
 * 1. P7_OPROFILE: a scoring profile
 *****************************************************************/

/* Indices for special state types in the length model, om->xsc[x][]
 */
enum p7o_xstates_e { 
  p7O_N = 0,
  p7O_E = 1,
  p7O_C = 2,
  p7O_J = 3
};
#define p7O_NXSTATES 4

/* Indices for transitions from the length modeling scores om->xsc[][x]
 */
enum p7o_xtransitions_e {
  p7O_MOVE = 0,
  p7O_LOOP = 1
};
#define p7O_NXTRANS 2

/* Indices x for transition scores gm->tsc[x][k] */
/* order is unoptimized, same as in core HMM */
/* also note that H2 erroneously has transitions stored in bad cache stride */
enum p7o_tsc_e {
  p7O_MM = 0, 
  p7O_MI = 1, 
  p7O_MD = 2, 
  p7O_IM = 3, 
  p7O_II = 4, 
  p7O_DM = 5, 
  p7O_DD = 6, 
};
#define p7O_NTRANS 7

typedef struct p7_oprofile_s {
  int    mode;         	/* configured algorithm mode (e.g. p7_LOCAL)   */ 
  int    M;		/* number of nodes in the model                */
  int  **tsc;           /* transition scores     [0.6][1.M-1]          */
  int  **msc;           /* match emission scores [0.Kp-1][1.M]         */
  int  **isc;           /* ins emission scores   [0.Kp-1][1.M-1]       */
  int    xsc[p7O_NXSTATES][p7O_NXTRANS]; /* specials [NECJ][LOOP,MOVE] */
  int   *bsc;           /* begin transitions     [1.M]                 */
  int   *esc;		/* end transitions       [1.M]                 */

  const ESL_ALPHABET    *abc;	/* copy of ptr to appropriate alphabet */
} P7_OPROFILE;


/*****************************************************************
 * 2. P7_OMX: a dynamic programming matrix
 *****************************************************************/

enum p7x_xcells_e {
  p7X_B  = 0,
  p7X_E  = 1,
  p7X_C  = 2,
  p7X_J  = 3,
  p7X_N  = 4
};
#define p7X_NXCELLS 5

typedef struct p7_omx_s {
  int  M;		/* actual model dimension (model 1..M)    */
  int  L;		/* actual sequence dimension (seq 1..L)   */
  
  size_t ncells;	/* current cell allocation limit: >= (M+1)*(L+1) */
  size_t nrows;    	/* current row allocation limit:  >= L+1  */

  int **xmx;		/* special scores [0.1..L][BECJN]      */
  int **mmx;		/* match scores   [0.1..L][0.1..M]     */
  int **imx;		/* insert scores  [0.1..L][0.1..M-1.M] */
  int **dmx;		/* delete scores  [0.1..L][0.1..M-1.M] */

  int  *xmx_mem;	
  int  *mmx_mem;
  int  *imx_mem;
  int  *dmx_mem;
} P7_OMX;



extern P7_OPROFILE *p7_oprofile_Create(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy(P7_OPROFILE *om);

extern P7_OMX      *p7_omx_Create(int allocM, int allocL);
extern int          p7_omx_GrowTo(P7_OMX *ox, int allocM, int allocL);
extern void         p7_omx_Destroy(P7_OMX *ox);
extern int          p7_omx_Dump(FILE *ofp, P7_OMX *ox);

extern int          p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om);

extern int          p7_Viterbi(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int          p7_Forward(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);


#endif /*P7_IMPL_H2_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
