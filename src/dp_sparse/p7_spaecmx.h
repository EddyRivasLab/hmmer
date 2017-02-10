#ifndef p7SPAECMX_INCLUDED
#define p7SPAECMX_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#include "base/p7_envelopes.h"

#include "dp_sparse/p7_sparsemx.h"

extern int p7_spaecmx_Dump(FILE *fp, const P7_SPARSEMX *aec, const P7_ENVELOPES *env);
extern int p7_spaecmx_Validate(const P7_SPARSEMX *aec, const P7_ENVELOPES *env, char *errbuf);


/*****************************************************************
 * 2. Footnotes
 *****************************************************************
 *
 * [1] ON SPECIAL CASE VALUES IN SPARSE AEC MATRIX
 * 
 * '*' marks cells prohibited by construction, which must have values
 * of -inf. Compare notes in p7_sparsemx.h, p7_spascmx.h.
 * 
 *          [M I D]..[M I D]..[M I D]  [M I D]  [M I D]..[M I D]..[M I D]
 *  k:         1        k      k0-1      k0      k0+1       k        M        
 *    ------------------------------------------------------------------
 *    ia   | . * *    . * .    . * . |  <= UP sector
 *    ..   | . . *   {. . .}   . . . |
 *    i0-1 | . . *    . . .    . . . |
 *    ------------------------------------------------------------------
 *    i0   |                         |  . * *    * * .    * * .    * * .
 *    ..   |         DOWN sector =>  |  * . *    . . *   {. . .}   . * .
 *    ib   |                         |  * . *    . . *    . . .    . * . 
 *    
 */
#endif /*p7SPAECMX_INCLUDED*/

