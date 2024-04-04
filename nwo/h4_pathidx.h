/* H4_PATHIDX - index of the locations of domains in a path
 *
 * See h4_pathidx.md for notes.
 */
#ifndef h4PATHIDX_INCLUDED
#define h4PATHIDX_INCLUDED

#include <h4_config.h>

#include "easel.h"

#include "h4_path.h"

/* H4_PATHIDX
 */
typedef struct {
  int *ia, *ib;     // ia[d],ib[d]: domain d start,end on target seq; 1..L  (d = 1..D)
  int *ka, *kb;     // ka[d],kb[d]: domain d start,end on profile;    1..M
  int *za, *zb;     // za[d],zb[d]: domain d start|end on path;       1..Z (start on L|G; end on last Mk|Dk)
  int *is_glocal;   // is_glocal[d] = TRUE if G, FALSE if L
  int  D;           // number of domains in path (d=1..D)
  int  L;           // target sequence length deduced from path
} H4_PATHIDX;


extern int  h4_pathidx_Build(const H4_PATH *pi, H4_PATHIDX **ret_pidx);
extern void h4_pathidx_Destroy(H4_PATHIDX *pidx);

extern int  h4_pathidx_Validate(const H4_PATHIDX *pidx, const H4_PATH *pi, char *errbuf);


#endif // h4PATHIDX_INCLUDED
