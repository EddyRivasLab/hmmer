#ifndef P7_COORDS2_INCLUDED
#define P7_COORDS2_INCLUDED

#include "p7_config.h"

typedef struct {
  int start;
  int end;
} P7_COORD2;

typedef struct {
  P7_COORD2 *seg;
  int        nseg;
  int        L;

  int        nalloc;		/* current allocation size for <seg> */
  int        nredline;		/* Reuse() will pull a large allocation back down to this */
} P7_COORDS2;

#define p7COORDS2_INITALLOC_DEFAULT 8
#define p7COORDS2_REDLINE_DEFAULT   64

extern P7_COORDS2 *p7_coords2_Create (int nalloc, int nredline);
extern int         p7_coords2_Grow   (P7_COORDS2 *c2);
extern int         p7_coords2_GrowTo (P7_COORDS2 *c2, int nalloc);
extern int         p7_coords2_Reuse  (P7_COORDS2 *c2);
extern void        p7_coords2_Destroy(P7_COORDS2 *c2);

#endif /* P7_COORDS2_INCLUDED */

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

  
  


