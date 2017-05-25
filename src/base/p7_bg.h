/* P7_BG: a null (background) model.
 */

/* This really contains three different things: 
 *     
 *   - the "null1" model, a one-state HMM consisting of background
 *     frequencies <f> and a parameter <p1> for a target-length
 *     dependent geometric;
 *     
 *   - the "bias filter" <fhmm> a two-state HMM composed from null1's
 *     background <f> and the model's mean composition <compo>. This
 *     model is constructed dynamically, every time a new profile is 
 *     considered;
 *     
 *   - a single term <omega> that's needed by the "null2" model to set
 *     a balance between the null1 and null2 scoring terms.  The null2
 *     model is otherwise defined by construction, in p7_domaindef.c.
 *
 * Someday we might pull this apart into two or three separate
 * objects.
 */
#ifndef p7BG_INCLUDED
#define p7BG_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_hmm.h"


typedef struct p7_bg_s {
  float   *f;		/* null1 background residue frequencies [0..K-1]: set at initialization    */
  float    p1;		/* null1's transition prob: p7_bg_SetLength() sets this from target seq L  */

  ESL_HMM *fhmm;	/* bias filter: p7_bg_SetFilter() sets this, from model's mean composition */

  float    omega;	/* the "prior" on null2/null3: set at initialization (one omega for both null types)  */
  int      use_null3;  /* use null3 in addition to null2 ?*/

  const ESL_ALPHABET *abc;	/* reference to alphabet in use: set at initialization             */
} P7_BG;


extern P7_BG *p7_bg_Create(const ESL_ALPHABET *abc);
extern P7_BG *p7_bg_CreateUniform(const ESL_ALPHABET *abc);
extern P7_BG *p7_bg_Clone(const P7_BG *bg);
extern int    p7_bg_Dump(FILE *ofp, const P7_BG *bg);
extern void   p7_bg_Destroy(P7_BG *bg);
extern int    p7_bg_SetLength(P7_BG *bg, int L);
extern int    p7_bg_NullOne(const P7_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc);

extern int    p7_bg_Read(char *bgfile, P7_BG *bg, char *errbuf);
extern int    p7_bg_Write(FILE *fp, P7_BG *bg);

extern int    p7_bg_SetFilter  (P7_BG *bg, int M, const float *compo);
extern int    p7_bg_FilterScore(P7_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc);


#endif /*p7BG_INCLUDED*/

