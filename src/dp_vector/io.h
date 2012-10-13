#ifndef P7_IO_INCLUDED
#define P7_IO_INCLUDED

#include "p7_config.h"
#include <stdio.h"

#include "esl_alphabet.h"

#include "dp_vector/p7_oprofile.h"
#include "p7_hmmfile.h"

extern int p7_oprofile_Write(FILE *ffp, FILE *pfp, P7_OPROFILE *om);
extern int p7_oprofile_ReadMSV (P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
extern int p7_oprofile_ReadInfoMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
extern int p7_oprofile_ReadBlockMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OM_BLOCK *hmmBlock);
extern int p7_oprofile_ReadRest(P7_HMMFILE *hfp, P7_OPROFILE *om);
extern int p7_oprofile_Position(P7_HMMFILE *hfp, off_t offset);

extern P7_OM_BLOCK *p7_oprofile_CreateBlock(int size);
extern void p7_oprofile_DestroyBlock(P7_OM_BLOCK *block);

#endif /*P7_IO_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

