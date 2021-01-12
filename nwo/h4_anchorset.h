
#ifndef h4ANCHORSET_INCLUDED
#define h4ANCHORSET_INCLUDED

#include "h4_config.h"


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
 *     a[0].i0,k0   = (0,0)
 *     a[D+1].i0,k0 = (L+1,M+1)
 *
 *   Therefore, for D domains, anch->a[] is always allocated for at
 *   least D+2 <H4_ANCHOR> structures, (0) 1..D (D+1).
 */
typedef struct {
  H4_ANCHOR *a;         // array of anchor (i0,k0) coords: (0) 1..D (D+1)
  int32_t    D;         // number of domains in <a>, exclusive of 2 sentinels

  int32_t    nalloc;    // current allocation size for <a>; >= D+2 because of 2 sentinels
  int32_t    nredline;  // _Reuse() call pulls allocation back down to this
} H4_ANCHORSET;


extern H4_ANCHORSET *h4_anchorset_Create (int D, int L, int M);
extern int           h4_anchorset_Add    (H4_ANCHORSET *anch, int i0, int k0);
extern int           h4_anchorset_GrowFor(H4_ANCHORSET *anch, int D);
extern int           h4_anchorset_Reuse  (H4_ANCHORSET *anch);
extern void          h4_anchorset_Destroy(H4_ANCHORSET *anch);

extern int           h4_anchorset_Dump(FILE *fp, H4_ANCHORSET *anch);
extern int           h4_anchorset_Validate(H4_ANCHORSET *anch, char *errmsg);

#endif // h4ANCHORSET_INCLUDED
