/* Optimized implementation of HMMER2 from  Jeremy Buhler (Washington
 * University, St. Louis); adapted to HMMER3 and converted to floating
 * point operation.
 * 
 * SRE, Wed Jun 27 10:08:48 2007 [Janelia]
 * SVN $Id$
 */
#ifndef P7_IMPL_FP_INCLUDED
#define P7_IMPL_FP_INCLUDED

/* special transition scores are reordered relative to generic profiles */
#define p7X_NX   4
enum p7X_xsc_e { p7X_XTE = 0, p7X_XTN = 1, p7X_XTJ = 2, p7X_XTC = 3 };

#define p7X_NXT  2
enum p7X_xmove_e { p7X_LOOP = 0, p7X_MOVE = 1 };

/* transition scores reordered relative to generic profiles. */
#define p7X_NT 8
enum p7X_tsc_e { p7X_TMM = 0, p7X_TIM = 1, p7X_TDM = 2, 
               p7X_TMD = 3, p7X_TDD = 4, 
               p7X_TMI = 5, p7X_TII = 6, 
               p7X_BSC = 7 };

#define p7X_NR 2
enum p7X_rsc_e { p7X_MSC = 0, p7X_ISC = 1 };


typedef struct p7_oprofile_s {
  int     M;
  float  *tsc;	/* [0.1..M-1]0..p7X_NTSC-1] */
  float **rsc;	/* [0..Kp-1][0.1..M][p7X_NR] */
  float   xsc[p7X_NX][p7X_NXT];

  const ESL_ALPHABET *abc_r;
} P7_OPROFILE;

#define p7X_NSCELLS 3
enum p7X_scells_e { p7X_XMM = 0, p7X_XMI = 1, p7X_XMD = 2 };

#define p7X_NXCELLS 5
enum p7X_xcells_e { p7X_XME = 0, p7X_XMN = 1, 
		    p7X_XMJ = 2, p7X_XMB = 3, 
		    p7X_XMC = 4 };

typedef struct p7_omx_s {
  int M;
  int L;

  size_t ncells;
  size_t nrows;
  
  float **dp;			/* [0.1..L][0.1..M][0..p7X_NSCELLS-1] */
  float *xmx;              	/* [0.1..L][0..p7X_NXCELLS-1] */
} P7_OMX;


extern P7_OPROFILE *p7_oprofile_Create(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy(P7_OPROFILE *om);
extern P7_OMX      *p7_omx_Create(int allocM, int allocL);
extern void         p7_omx_Destroy(P7_OMX *ox);
extern int          p7_oprofile_Config(const P7_HMM *hmm, const P7_BG *bg, P7_OPROFILE *om, int mode);
extern int          p7_oprofile_ReconfigLength(P7_OPROFILE *om, int L);
extern int          p7_Viterbi(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);

#endif /*P7_IMPL_FP_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
