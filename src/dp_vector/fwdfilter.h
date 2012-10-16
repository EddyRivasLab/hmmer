#ifndef P7_FWDFILTER_INCLUDED
#define P7_FWDFILTER_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_sparse/p7_sparsemx.h"

extern int p7_ForwardFilter (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int p7_BackwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm);


#endif /*P7_FWDFILTER_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
