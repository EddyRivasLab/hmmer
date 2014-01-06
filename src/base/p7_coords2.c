/* P7_COORDS2 maintains a resizeable array of start/end coord pairs. 
 * Used for domain segment coords in a target sequence, for example.
 * 
 * Contents:
 *   x. P7_COORDS2 object.
 *   x. Copyright and license information.
 */
#include "p7_config.h"

#include "base/p7_coords2.h"

P7_COORDS2 *
p7_coords2_Create(int nalloc, int nredline)
{
  P7_COORDS2 *c = NULL;
  int         status;

  if (nalloc   <= 0) nalloc   = p7COORDS2_INITALLOC_DEFAULT;
  if (nredline <= 0) nredline = p7COORDS2_REDLINE_DEFAULT;

  ESL_ALLOC(c2, sizeof(P7_COORDS2));
  c2->seg = NULL;
  
  ESL_ALLOC(c2->seg, sizeof(P7_COORD2) * nalloc);

  c2->nalloc   = nalloc;
  c2->nredline = nredline;

  c2->nseg = 0;
  c2->L    = 0;
  return c2;

 ERROR:
  p7_coords2_Destroy(c2);
  return NULL;
}

int
p7_coords2_Grow(P7_COORDS2 *c2)
{
  int status;

  if (c2->nseg < c2->nalloc) return eslOK;

  ESL_REALLOC(c2->seg, sizeof(P7_COORD2) * c2->nalloc * 2);
  c2->nalloc = c2->nalloc * 2;
  return eslOK;

 ERROR:
  return status;
}

int
p7_coords2_Reuse(P7_COORDS2 *c2)
{
  int status;

  if (c2->nalloc > c2->nredline) 
    {
      ESL_REALLOC(c2->seg, sizeof(P7_COORD2) * c2->nredline);
      c2->nalloc = c2->nredline;
    }

  c2->nseg = 0;
  c2->L    = 0;
  return eslOK;

 ERROR:
  return status;
}

void
p7_coords2_Destroy(P7_COORDS2 *c2)
{
  if (c2) 
    {
      if (c2->seg) free(c2->seg);
      free(c2);
    }
  return;
}


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

