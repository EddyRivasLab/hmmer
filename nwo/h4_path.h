#ifndef h4PATH_INCLUDED
#define h4PATH_INCLUDED

#include "easel.h"
#include "esl_alphabet.h"

/* Codes used in H4_PATH pi->st[] */
enum h4_statetypes_e {
  h4_NONE = 0,
  h4_S    = 1,
  h4_N    = 2,
  h4_B    = 3,
  h4_G    = 4, 
  h4_MG   = 5,
  h4_IG   = 6,
  h4_DG   = 7,
  h4_L    = 8,
  h4_ML   = 9,
  h4_IL   = 10,
  h4_DL   = 11,
  h4_E    = 12,
  h4_J    = 13,
  h4_C    = 14,
  h4_T    = 15
};

typedef struct {
  int     Z;      // length of <st>, <r> (actual; i.e. inclusive of run length compression)
  int8_t *st;     // state codes
  int    *rle;    // run lengths of st[z]

  int     Zalloc;   // current allocation for st, r
  int     Zredline; // Reuse() will downalloc to this, if exceeded.
} H4_PATH;

extern H4_PATH *h4_path_Create (void);
extern int      h4_path_Grow   (H4_PATH *pi);
extern int      h4_path_Append (H4_PATH *pi, int8_t st);
extern int      h4_path_Reuse  (H4_PATH *pi);
extern void     h4_path_Destroy(H4_PATH *pi);

extern int h4_path_InferLocal (ESL_ALPHABET *abc, ESL_DSQ *ax, int alen, int8_t *matassign, H4_PATH *pi);
extern int h4_path_InferGlocal(ESL_ALPHABET *abc, ESL_DSQ *ax, int alen, int8_t *matassign, H4_PATH *pi);

extern char *h4_path_DecodeStatetype(int8_t st);
extern int   h4_path_Dump(FILE *fp, const H4_PATH *pi);







#endif // h4PATH_INCLUDED


  
