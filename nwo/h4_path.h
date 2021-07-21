#ifndef h4PATH_INCLUDED
#define h4PATH_INCLUDED
#include "h4_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"

#include "h4_counts.h"
#include "h4_mode.h"
#include "h4_profile.h"


typedef struct {
  int     Z;        // length of <st>, <rle> (actual; i.e. inclusive of run length compression)
  int8_t *st;       // state codes
  int    *rle;      // run lengths of st[z];  or, for st[z] == L, k for L->Mk entry. (for G, rle[z]=1, and some code depends on this)

  int     Zalloc;   // current allocation for st, rle
  int     Zredline; // Reuse() will downalloc to this, if exceeded.
} H4_PATH;

/* Codes for states, esp. states used in H4_PATH pi->st[] */
/* Do not change order. Some routines make assumptions about it: for
 * example, h4_path_TestSample() and h4_emit() assume that MG-IG-DG
 * and ML-IL-DL are contiguous; tracebacks can
 * assume G-MG-IG-DG, L-ML-IL-DL contiguity and order.
 */
#define h4P_NONE  0
#define h4P_S     1   // unused in paths
#define h4P_N     2   
#define h4P_B     3   // unused
#define h4P_G     4
#define h4P_MG    5
#define h4P_IG    6
#define h4P_DG    7
#define h4P_L     8
#define h4P_ML    9
#define h4P_IL    10
#define h4P_DL    11
#define h4P_E     12  // unused
#define h4P_J     13
#define h4P_C     14  
#define h4P_T     15  // unused
#define h4P_NST   16  

#define h4_path_IsB(s)     ( (s) == h4P_G  || (s) == h4P_L  )
#define h4_path_IsMain(s)  ( (s) >= h4P_MG || (s) <= h4P_DL )
#define h4_path_IsM(s)     ( (s) == h4P_MG || (s) == h4P_ML )
#define h4_path_IsI(s)     ( (s) == h4P_IG || (s) == h4P_IL )
#define h4_path_IsD(s)     ( (s) == h4P_DG || (s) == h4P_DL )
#define h4_path_IsX(s)     ( (s) == h4P_N  || (s) == h4P_J || (s) == h4P_C )

/* 1. H4_PATH structure */
extern H4_PATH *h4_path_Create (void);
extern H4_PATH *h4_path_Clone  (const H4_PATH *pi);
extern int      h4_path_Copy   (const H4_PATH *src, H4_PATH *dst);
extern int      h4_path_Resize (H4_PATH *pi, int Z);
extern int      h4_path_Grow   (H4_PATH *pi);
extern int      h4_path_Append (H4_PATH *pi, int8_t st);
extern int      h4_path_AppendElement(H4_PATH *pi, int8_t st, int r);
extern int      h4_path_AppendSeveral(H4_PATH *pi, int8_t st, int n);
extern int      h4_path_Reverse(H4_PATH *pi);
extern int      h4_path_Reuse  (H4_PATH *pi);
extern void     h4_path_Destroy(H4_PATH *pi);

/* 2. Getting info from H4_PATH */
extern int h4_path_GetSeqlen(const H4_PATH *pi);
extern int h4_path_GetDomainCount(const H4_PATH *pi);
extern int h4_path_FetchDomainBounds(const H4_PATH *pi, int whichd, int *opt_ia, int *opt_ib, int *opt_ka, int *opt_kb);

/* 3. Inferring paths from existing alignments */
extern int h4_path_InferLocal (const ESL_ALPHABET *abc, const ESL_DSQ *ax, int alen, const int8_t *matassign, int lcol, int rcol, H4_PATH *pi);
extern int h4_path_InferGlocal(const ESL_ALPHABET *abc, const ESL_DSQ *ax, int alen, const int8_t *matassign, int lcol, int rcol, H4_PATH *pi);

/* 4. Counting paths into new HMMs */
extern int h4_path_Count(const H4_PATH *pi, const ESL_DSQ *dsq, float wgt, H4_COUNTS *ctm);

/* 5. Calculating lod scores for paths */
extern int h4_path_Score(const H4_PATH *pi, const ESL_DSQ *dsq, const H4_PROFILE *hmm, const H4_MODE *mo, float *ret_sc);

/* 6. Debugging and development tools */
extern int   h4_path_Example(H4_PATH **ret_pi);
extern int   h4_path_TestSample(ESL_RANDOMNESS *rng, H4_PATH **ret_pi);
extern char *h4_path_DecodeStatetype(int8_t st);
extern int   h4_path_Validate(const H4_PATH *pi, int M, int L, char *errbuf);
extern int   h4_path_Compare(const H4_PATH *pi1, const H4_PATH *pi2);
extern int   h4_path_Dump     (FILE *fp, const H4_PATH *pi);
extern int   h4_path_DumpCigar(FILE *fp, const H4_PATH *pi);






#endif // h4PATH_INCLUDED


  
