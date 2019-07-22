/* H4_CHECKPTMX: checkpointed, striped DP matrix for F/B filter.
 * 
 * Independent of vector ISA. Do not add any ISA-specific code.
 * 
 * See h4_checkptmx.md for notes.
 */
#ifndef h4CHECKPTMX_INCLUDED
#define h4CHECKPTMX_INCLUDED
#include "h4_config.h"


#if eslDEBUGLEVEL > 0
#include "h4_refmx.h"
#endif


/* H4_CHECKPTMX
 * Checkpointed DP matrix used by the F/B filter
 */
typedef struct {
  /* Aligned memory allocation (for all rows): */
  int M;               // current actual query model dimension (consensus positions)
  int L;               // current actual target seq dimension (residues)
  int Vf;              // width of striped vectors, in # floats. 0: no striping set yet
  int Q;               // current actual number of fb vectors per row = H4_Q(M,Vf)
  int R;               // current actual number of rows (<=Ra+Rb+Rc), excluding R0 
  float   *dp_mem;     // raw memory allocation, that dp[] rows point into
  int64_t  allocW;     // alloced width/row, cells (floats); multiple of h4_VMAX_FB 
  int64_t  allocN;     // total # of alloc'ed cells: allocN >= (validR)(allocW) 
  int64_t  redline;    // recommended allocN limit (floats); can temporarily exceed it

  /* Forward/Backward matrix rows: */
  float **dpf;         // row ptrs, dpf[0.R0-1,R0..R0+R-1], aligned memory
  int     allocR;      // allocated dpf[]. R+R0 <= R0+Ra+Rb+Rc <= validR <= allocR
  int     validR;      // # dpf[] rows on valid dp_mem; maybe < allocR after Reinit()

  /* Checkpointed layout, mapping DP rows 1..R to residues 1..L: */
  int R0;              // # of extra rows: one for fwd[0] boundary, two for bck[prv,cur]
  int Ra;              // # of rows used in "all" region (uncheckpointed)
  int Rb;              // # of rows in "between" region (one incomplete checkpoint segment)
  int Rc;              // # of rows in "checkpointed" region
  int La;              // residues 1..La are in "all" region
  int Lb;              // residues La+1..La+Lb are in "between" region
  int Lc;              // residues La+Lb+1..La+Lb+Lc=L are in "checkpointed" region

#if eslDEBUGLEVEL > 0
  /* Info for dumping debugging info, conditionally compiled */
  int       do_dumping;         // TRUE if matrix is in dumping mode
  FILE     *dfp;                // open output stream for debug dumps
  int       dump_maxpfx;        // each line prefixed by tag of up to this # chars
  int       dump_width;         // cell values in diagnostic output are fprintf'ed:
  int       dump_precision;     //   dfp, "%*.*f", dbg_width, dbg_precision, val

  H4_REFMX *fwd;                // full Forward matrix, saved for unit test diffs
  H4_REFMX *bck;                // ... full Backward matrix, ditto
  H4_REFMX *pp;                 // ... full posterior probability matrix, ditto
  float     bcksc;              // Backwards score: which we check against Forward
#endif // eslDEBUGLEVEL
} H4_CHECKPTMX;


/* indices for main MID states in a DP row */
#define h4C_M       0
#define h4C_I       1
#define h4C_D       2
#define h4C_NSCELLS 3

/* indices for special states/values at end of a DP row */
#define h4C_E       0
#define h4C_N       1
#define h4C_JJ      2
#define h4C_J       3
#define h4C_B       4
#define h4C_CC      5
#define h4C_C       6
#define h4C_SCALE   7
#define h4C_NXCELLS 8

extern H4_CHECKPTMX *h4_checkptmx_Create(int M, int L, int64_t redline);
extern int           h4_checkptmx_Reinit(H4_CHECKPTMX *cpx, int M, int L);
extern size_t        h4_checkptmx_Sizeof(const H4_CHECKPTMX *cpx);
extern size_t        h4_checkptmx_MinSizeof(int M, int L);
extern void          h4_checkptmx_Destroy(H4_CHECKPTMX *cpx);


#endif // h4CHECKPTMX_INCLUDED
