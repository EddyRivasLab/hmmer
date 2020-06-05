/* H4_FILTERMX: one row of DP memory for Viterbi filter.
 * 
 * Independent of vector ISA. Do not add any ISA-specific code.
 */
#ifndef h4FILTERMX_INCLUDED
#define h4FILTERMX_INCLUDED
#include "h4_config.h"

#include <stdio.h>

#include "simdvec.h"


/* H4_FILTERMX
 * One row vectorized DP matrix used by Viterbi filter.
 */
typedef struct {  
  int      M;           // current profile size: determines width of <dp> row    
  int      Vw;          // width of striped vectors, in int16_t's.  0 if no striped data set yet

  int16_t *dp;          // aligned, one row of DP memory: >= 3*Qf vectors. 
  int      allocM;      // <dp_mem> is allocated to hold up to M=allocM

#if eslDEBUGLEVEL > 0
  FILE    *dfp;         // if non-NULL, dump matrix rows to stream <dfp>
#endif
} H4_FILTERMX;


/* VF matrix has three cells per vector supercell. */
#define h4F_M      0
#define h4F_I      1
#define h4F_D      2
#define h4F_NCELLS 3

/* VF implementations use these macros to improve clarity of accesses.
 * However, a variable named <dp> must be set to <dp = fx->dp> to use them,
 * for an H4_FILTERMX <fx>.
 */
#define MMXf(q)   (dp[(q) * h4F_NCELLS + h4F_M])
#define IMXf(q)   (dp[(q) * h4F_NCELLS + h4F_I])
#define DMXf(q)   (dp[(q) * h4F_NCELLS + h4F_D])


extern H4_FILTERMX *h4_filtermx_Create  (int allocM);
extern int          h4_filtermx_Reinit  (H4_FILTERMX *fx, int allocM);
extern int          h4_filtermx_GetScore(const H4_FILTERMX *fx, int k, int s);
extern int          h4_filtermx_Reuse    (H4_FILTERMX *fx);
extern void         h4_filtermx_Destroy  (H4_FILTERMX *fx);

#if eslDEBUGLEVEL > 0
extern int h4_filtermx_SetDumpMode(H4_FILTERMX *fx, FILE *dfp);
extern int h4_filtermx_DumpRow(const H4_FILTERMX *fx, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
#endif // eslDEBUGLEVEL > 0

#endif // h4FILTERMX_INCLUDED
