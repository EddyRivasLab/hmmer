/* h4_hmmfile : profile HMM input and output
 */
#ifndef h4HMMFILE_INCLUDED
#define h4HMMFILE_INCLUDED
#include "h4_config.h"

#include <stdio.h>

#include "esl_buffer.h"
#include "esl_json.h"
#include "esl_keyhash.h"

#include "h4_profile.h"


typedef struct {
  ESL_BUFFER      *bf;    // open input stream
  ESL_JSON_PARSER *jprs;  // precise state of JSON parse after last byte read
  ESL_JSON        *pi;    // JSON parse tree
  ESL_KEYHASH     *kh;    // keyhash of top level keys ("length", "alphabet", "match", "t", etc.)

  int             *tokmap;    // map keyhash indices to token indices. For <keyi> from <kh>, tokmap[keyi] = idx in pi->tok[idx]:
  int              nmap;      // = the number of unique keys in an H4 profile file. 
  int              nmapalloc;

  char  errmsg[eslERRBUFSIZE]; // user-directed error message on parse failure
} H4_HMMFILE;


extern int  h4_hmmfile_Open(const char *filename, const char *envvar, H4_HMMFILE **ret_hfp);
extern void h4_hmmfile_Close(H4_HMMFILE *hfp);

extern int  h4_hmmfile_Read(H4_HMMFILE *hfp, ESL_ALPHABET **ret_abc, H4_PROFILE **opt_hmm);

extern int  h4_hmmfile_Write(FILE *fp, const H4_PROFILE *hmm);


#endif /* h4HMMFILE_INCLUDED */

