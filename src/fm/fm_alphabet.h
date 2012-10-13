#ifndef P7_FM_ALPHABET_INCLUDED
#define P7_FM_ALPHABET_INCLUDED

#include "p7_config.h"

#include "fm/fm.h"

extern int fm_createAlphabet (FM_METADATA *meta, uint8_t *alph_bits);
extern int fm_reverseString (char* str, int N);
extern int fm_getComplement (char c, uint8_t alph_type);

#endif /*P7_FM_ALPHABET_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

