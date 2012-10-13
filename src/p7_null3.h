#ifndef P7_NULL3_INCLUDED
#define P7_NULL3_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"

#include "base/p7_bg.h"
#include "base/p7_trace.h"

extern void p7_null3_score(const ESL_ALPHABET *abc, const ESL_DSQ *dsq, P7_TRACE *tr, int start, int stop, P7_BG *bg, float *ret_sc);
extern void p7_null3_windowed_score(const ESL_ALPHABET *abc, const ESL_DSQ *dsq, int start, int stop, P7_BG *bg, float *ret_sc);

#endif /*P7_NULL3_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
