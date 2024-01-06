#ifndef p7IO_INCLUDED
#define p7IO_INCLUDED
#include <p7_config.h>

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"

#include "base/p7_hmmfile.h"
#include "dp_vector/p7_oprofile.h"

extern int p7_oprofile_Write       (FILE *ffp, FILE *pfp, P7_OPROFILE *om);
extern int p7_oprofile_ReadSSV     (P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
extern int p7_oprofile_ReadInfoSSV (P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
extern int p7_oprofile_ReadRest    (P7_HMMFILE *hfp, P7_OPROFILE *om);
extern int p7_oprofile_ReadBlockSSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OM_BLOCK *hmmBlock);

extern P7_OM_BLOCK *p7_oprofile_CreateBlock(int size);
extern void         p7_oprofile_DestroyBlock(P7_OM_BLOCK *block);
extern int          p7_oprofile_Position(P7_HMMFILE *hfp, off_t offset);

#endif /*p7IO_INCLUDED*/


