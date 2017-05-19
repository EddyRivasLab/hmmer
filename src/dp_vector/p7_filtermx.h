/* P7_FILTERMX: one-row of DP memory for MSV and Viterbi filters.
 * 
 * Independent of vector ISA. (Do not add any ISA-specific code.)
 */
#ifndef p7FILTERMX_INCLUDED
#define p7FILTERMX_INCLUDED
#include "p7_config.h"

#include <stdio.h>

#include "simdvec.h"

enum p7f_mxtype_e { p7F_NONE = 0, p7F_SSVFILTER = 1, p7F_MSVFILTER = 2, p7F_VITFILTER = 3 };
                                
typedef struct p7_filtermx_s {  
  int      M;			// current profile size: determines width of <dp> row    
  int      V;                   // width of striped vectors. 0 if no striped data set yet

                                // SRE REVISIT: char? only because it's both int8_t, int16_t
  char    *dp;			// aligned, one row of DP memory: >= 3*Qf vectors 
  int      allocM;		// <dp_mem> is allocated to hold up to M=allocM

  enum p7f_mxtype_e type;	// p7F_NONE | p7F_SSVFILTER | p7F_MSVFILTER | p7F_VITFILTER 

#if eslDEBUGLEVEL > 0
  int      do_dumping;		// TRUE if we're dumping each row for diagnostics 
  FILE    *dfp;			// where we're dumping it (or NULL, if debugging=FALSE) 
#endif
} P7_FILTERMX;


/* Viterbi filter has three cells per vector supercell. MSV/SSV only have match states. */
enum p7f_scells_e { p7F_M = 0, p7F_D = 1, p7F_I = 2 };
#define p7F_NSCELLS 3

/* Viterbi filter uses these macros to improve clarity of accesses */
#define MMXf(q)   (dp[(q) * p7F_NSCELLS + p7F_M])
#define DMXf(q)   (dp[(q) * p7F_NSCELLS + p7F_D])
#define IMXf(q)   (dp[(q) * p7F_NSCELLS + p7F_I])

extern P7_FILTERMX *p7_filtermx_Create(int allocM);
extern int          p7_filtermx_Reinit(P7_FILTERMX *fx, int allocM);
extern int          p7_filtermx_GetScore(const P7_FILTERMX *fx, int k, enum p7f_scells_e s);
extern size_t       p7_filtermx_Sizeof(const P7_FILTERMX *fx);
extern size_t       p7_filtermx_MinSizeof(int M);
extern void         p7_filtermx_Destroy(P7_FILTERMX *fx);

#if eslDEBUGLEVEL > 0
extern int p7_filtermx_SetDumpMode(P7_FILTERMX *fx, FILE *dfp, int truefalse);
extern int p7_filtermx_DumpMFRow(const P7_FILTERMX *fx, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
extern int p7_filtermx_DumpVFRow(const P7_FILTERMX *fx, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
#endif

#endif /*p7FILTERMX_INCLUDED*/




