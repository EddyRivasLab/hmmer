#ifndef p7IO_INCLUDED
#define p7IO_INCLUDED

#include "p7_config.h"
#include <stdio.h>

#include "esl_alphabet.h"

#include "base/p7_hmmfile.h"
#include "hardware/hardware.h"
#include "dp_vector/p7_oprofile.h"


extern int p7_oprofile_Write(FILE *ffp, FILE *pfp, P7_OPROFILE *om);
extern int p7_oprofile_Write_sse(FILE *ffp, FILE *pfp, P7_OPROFILE *om);
extern int p7_oprofile_Write_avx(FILE *ffp, FILE *pfp, P7_OPROFILE *om);
extern int p7_oprofile_Write_avx512(FILE *ffp, FILE *pfp, P7_OPROFILE *om);

extern int p7_oprofile_ReadMSV (P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om, SIMD_TYPE simd);
extern int p7_oprofile_ReadMSV_sse (P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
extern int p7_oprofile_ReadMSV_avx (P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
extern int p7_oprofile_ReadMSV_avx512 (P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
extern int p7_oprofile_ReadMSV_neon (P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);

extern int p7_oprofile_ReadInfoMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om, SIMD_TYPE simd);
extern int p7_oprofile_ReadInfoMSV_sse(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
extern int p7_oprofile_ReadInfoMSV_avx(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
extern int p7_oprofile_ReadInfoMSV_avx512(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
extern int p7_oprofile_ReadInfoMSV_neon(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
extern int p7_oprofile_ReadBlockMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OM_BLOCK *hmmBlock, SIMD_TYPE simd);
extern int p7_oprofile_ReadRest(P7_HMMFILE *hfp, P7_OPROFILE *om);
extern int p7_oprofile_ReadRest_sse(P7_HMMFILE *hfp, P7_OPROFILE *om);
extern int p7_oprofile_ReadRest_avx(P7_HMMFILE *hfp, P7_OPROFILE *om);
extern int p7_oprofile_ReadRest_avx512(P7_HMMFILE *hfp, P7_OPROFILE *om);
extern int p7_oprofile_ReadRest_neon(P7_HMMFILE *hfp, P7_OPROFILE *om);
extern int p7_oprofile_Position(P7_HMMFILE *hfp, off_t offset);

extern P7_OM_BLOCK *p7_oprofile_CreateBlock(int size);
extern void p7_oprofile_DestroyBlock(P7_OM_BLOCK *block);

#endif /*p7IO_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

