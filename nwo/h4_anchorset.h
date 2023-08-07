#ifndef h4ANCHORSET_INCLUDED
#define h4ANCHORSET_INCLUDED
#include <h4_config.h>

#include <stdio.h>

#include "esl_random.h"

#include "h4_path.h"


/* H4_ANCHOR
 *   Just an integer pair (i0,k0), for H4_ANCHORSET to make an array of.
 */
typedef struct {
  int32_t i0;
  int32_t k0;
} H4_ANCHOR;


/* H4_ANCHORSET
 *   Is a memory-managed wrapper around an array of H4_ANCHOR's.
 *   MPAS algorithm needs to rapidly create and reuse many trial
 *   anchorsets, and the sparse segmental MPAS implementation needs to 
 *   concatenate anchorsets.
 *
 *   Anchors are indexed 1..D, with a[0] and a[D+1] used as sentinels.
 *     a[0].i0,k0   = (0,M+1)
 *     a[D+1].i0,k0 = (L+1,0)
 *
 *   Therefore, for D domains, anch->a[] is always allocated for at
 *   least D+2 <H4_ANCHOR> structures, (0) 1..D (D+1).
 *
 *   These particular sentinel values are required for two access
 *   patterns we use in ASC matrices. See h4_anchorset.md for
 *   explication.
 * 
 *   In an empty anchorset (D=0), sentinels are allowed to be (0,0),
 *   for cases where L,M are not yet known. Sentinels must be set
 *   whenever D>0.
 */ 
typedef struct {
  H4_ANCHOR *a;         // array of anchor (i0,k0) coords: (0) 1..D (D+1)
  int32_t    D;         // number of domains in <a>, exclusive of 2 sentinels

  int32_t    nalloc;    // current allocation size for <a>; >= D+2 because of 2 sentinels
  int32_t    nredline;  // _Reuse() call pulls allocation back down to this
} H4_ANCHORSET;


extern H4_ANCHORSET *h4_anchorset_Create (int D, int L, int M);
extern int           h4_anchorset_Add    (H4_ANCHORSET *anch, int i0, int k0);
extern int           h4_anchorset_GetSentinels(const H4_ANCHORSET *anch, int *ret_L, int *ret_M);
extern int           h4_anchorset_SetSentinels(H4_ANCHORSET *anch, int L, int M);
extern int           h4_anchorset_GrowFor(H4_ANCHORSET *anch, int D);
extern int           h4_anchorset_Copy(const H4_ANCHORSET *src, H4_ANCHORSET *dst);
extern int           h4_anchorset_Reuse  (H4_ANCHORSET *anch);
extern void          h4_anchorset_Destroy(H4_ANCHORSET *anch);


extern int h4_anchorset_Dump          (FILE *fp, const H4_ANCHORSET *anch);
extern int h4_anchorset_DumpOneLine   (FILE *fp, const H4_ANCHORSET *anch);
extern int h4_anchorset_Validate      (const H4_ANCHORSET *anch, char *errmsg);
extern int h4_anchorset_Compare       (const H4_ANCHORSET *anch1, const H4_ANCHORSET *anch2);
extern int h4_anchorset_Sample        (ESL_RANDOMNESS *rng, int L, int M, int maxD, H4_ANCHORSET *anch);
extern int h4_anchorset_SampleFromPath(ESL_RANDOMNESS *rng, const H4_PATH *pi, H4_ANCHORSET *anch);

#endif // h4ANCHORSET_INCLUDED
