/* Structure: P7_ALIDISPLAY
 * 
 * Alignment of a sequence domain to a profile, formatted for printing.
 * 
 * A homology domain produces a chunk of P7_TRACE:
 *     ... B -> {GL} -> {MD}k1 -> ... {MD}k2 -> E ...
 * There may be more than one per trace. Number them d, d=0..ndom-1.
 * For a given d, a trace's index tells us:
 *    tfrom[d]   = position of B in the trace's arrays (0..tr->N-1)
 *    tto[d]     = position of E (0..tr->N-1)
 *    sqfrom[d]  = position of first seq residue accounted for (1..L)
 *    sqto[d]    = position of last seq residue (1..L)
 *    hmmfrom[d] = k1
 *    hmmto[d]   = k2
 *    
 * A P7_ALIDISPLAY is an annotated text representation of such a chunk
 * of trace. From tfrom[d]+2 to tto[d]-1 (skipping the B, {GL}, and
 * {E}), we have a series of {MDI} states: each is converted to an
 * aligned symbol pair. These strings are tto[d]-1 - (tfrom[d]+2) + 1 = 
 * tto[d] - tfrom[d] - 2 symbols long. (That's ad->N.)
 * 
 * The alidisplay also records whether state tfrom[d]+1 was G or L
 * by setting the is_glocal flag. We need this, for instance, to backconvert
 * an alidisplay to a trace.
 * 
 * Memory in this structure may either be serialized or deserialized.
 * If serialized, ad->mem is non-NULL, ad->memsize is >0, and all the ptrs
 * point into that memory. If not, ad->mem is NULL, ad->memsize is 0,
 * and all the ptrs have their own allocation as NUL-terminated strings.
 *
 * For an alignment of L residues, names/accessions C chars long, descriptions
 * D characters long; requires 7L + 4C + 2D + 45 bytes; 
 * for typical case of L=100,C=10,D=60: 905B.
 * So design w/ expectation of ~1kb per structure.
 * 
 * Note that we could reduce that quite a bit, w/ a little more work;
 * there's a lot of redundancy in this information.
 */
#ifndef p7ALIDISPLAY_INCLUDED
#define p7ALIDISPLAY_INCLUDED
#include <p7_config.h>

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_mpi.h"		/* includes MPI_INT64_T, etc. on systems w/ MPI versions < 2.2 */
#include "esl_sq.h"

#include "base/p7_trace.h"
#include "base/p7_profile.h"

typedef struct p7_alidisplay_s {
  char *rfline;                 /* reference coord info; or NULL        */
  char *mmline;                 /* modelmask coord info; or NULL        */
  char *csline;                 /* consensus structure info; or NULL    */
  char *model;                  /* aligned query consensus sequence     */
  char *mline;                  /* "identities", conservation +'s, etc. */
  char *aseq;                   /* aligned target sequence              */
  char *ppline;		        /* posterior prob annotation; or NULL   */
  int   N;                      /* length of strings                    */
  char *hmmname;		/* name of HMM                          */
  char *hmmacc;			/* accession of HMM; or [0]='\0'        */
  char *hmmdesc;		/* description of HMM; or [0]='\0'      */
  int   hmmfrom;		/* start position on HMM (1..M, or -1)  */
  int   hmmto;			/* end position on HMM (1..M, or -1)    */
  int   M;			/* length of model                      */
  uint8_t is_glocal;		/* TRUE if this is a glocal alignment   */

  char *sqname;			/* name of target sequence              */
  char *sqacc;			/* accession of target seq; or [0]='\0' */
  char *sqdesc;			/* description of targ seq; or [0]='\0' */
  int64_t sqfrom;		/* start position on sequence (1..L)    */
  int64_t sqto;		        /* end position on sequence   (1..L)    */
  int64_t L;			/* length of sequence                   */

  int   memsize;                /* size of allocated block of memory    */
  char *mem;			/* memory used for the char data above  */
} P7_ALIDISPLAY;


/* 1. The P7_ALIDISPLAY object */
extern P7_ALIDISPLAY *p7_alidisplay_Create(const P7_TRACE *tr, int which, const P7_PROFILE *gm, const ESL_SQ *sq);
extern P7_ALIDISPLAY *p7_alidisplay_Clone(const P7_ALIDISPLAY *ad);
extern size_t         p7_alidisplay_Sizeof(const P7_ALIDISPLAY *ad);
extern int            p7_alidisplay_Serialize(P7_ALIDISPLAY *ad);
extern int            p7_alidisplay_Deserialize(P7_ALIDISPLAY *ad);
extern void           p7_alidisplay_Destroy(P7_ALIDISPLAY *ad);


/* 2. More P7_ALIDISPLAY API */
extern char           p7_alidisplay_EncodePostProb(float p);
extern float          p7_alidisplay_DecodePostProb(char pc);
extern char           p7_alidisplay_EncodeAliPostProb(float p, float hi, float med, float lo);
extern int            p7_alidisplay_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth, int show_accessions);
extern int            p7_alidisplay_Backconvert(const P7_ALIDISPLAY *ad, const ESL_ALPHABET *abc, ESL_SQ **ret_sq, P7_TRACE **ret_tr);

/* 3. Debugging and development tools */
extern int            p7_alidisplay_TestSample(ESL_RANDOMNESS *rng, int alen, P7_ALIDISPLAY **ret_ad);
extern int            p7_alidisplay_Dump(FILE *fp, const P7_ALIDISPLAY *ad);
extern int            p7_alidisplay_Validate(const P7_ALIDISPLAY *ad, char *errbuf);
extern int            p7_alidisplay_Compare(const P7_ALIDISPLAY *ad1, const P7_ALIDISPLAY *ad2);


#endif /*p7ALIDISPLAY_INCLUDED*/

