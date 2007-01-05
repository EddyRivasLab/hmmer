/* P7_HMMFILE: an HMM save file opened for reading.
 * 
 * SRE, Wed Jan  3 18:27:13 2007 [Casa de Gatos] [Leftover Salmon, Euphoria]
 * SVN $Id$
 */
#ifndef P7_HMMFILEH_INCLUDED
#define P7_HMMFILEH_INCLUDED

#include "p7_config.h"
#include <stdio.h>
#include "esl_alphabet.h"
#include "esl_ssi.h"
#include "p7_hmm.h"

/* P7_HMMFILE
 * An HMM save file, opened for reading.
 */
typedef struct p7_hmmfile_s {
  FILE         *f;		 /* pointer to stream for reading                */
  ESL_ALPHABET *abc;   		 /* ptr to alphabet in use for these HMMs        */
  int (*parser)(struct p7_hmmfile_s *, ESL_ALPHABET **, P7_HMM **);  /* parsing function */
} P7_HMMFILE;


extern int  p7_hmmfile_Open(char *filename, char *env, P7_HMMFILE **ret_hfp);
extern void p7_hmmfile_Close(P7_HMMFILE *hfp);

extern int  p7_hmmfile_Write(FILE *fp, P7_HMM *hmm);
extern int  p7_hmmfile_Read(P7_HMMFILE *hfp, ESL_ALPHABET **ret_abc,  P7_HMM **ret_hmm);

#endif /* P7_HMMFILE_INCLUDED */
/************************************************************
 * @LICENSE@
 ************************************************************/
