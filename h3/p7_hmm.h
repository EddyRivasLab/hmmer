/* The Plan7 HMM structure
 * 
 * SRE, Mon Jan  1 15:38:13 2007 [Casa de Gatos] [The Smiths, Hatful of Hollow]
 * SVN $Id$
 */
#ifndef P7_HMMH_INCLUDED
#define P7_HMMH_INCLUDED

#include "p7_config.h"

#include <esl_alphabet.h>	


/* P7_HMM
 * The core model, in counts or probability form.
 * Some notes:
 *   1. P7_HASPROBS flag is raised when t[][], mat[][], ins[][] contain normalized probabilities.
 *   2. t[0] is special: t[0][TMM] is the begin->M_1 entry probability, and t[0][TMD] 
 *      is the begin->D_1 entry probability. All other t[0] values are set to 0.
 *   3. Note that there is no insert state M.
 *   4. Note that there is no transition distribution t[M]; t[M][TMM] and
 *      t[M][TDM] are implicitly 1.0 to the end state E, and there is no insert
 *      state M.
 */
typedef struct {
  ESL_ALPHABET  *abc;		/* ptr to alphabet info (hmm->abc->K is alphabet size) */

  /*::cexcerpt::plan7_core::begin::*/
  int     M;                    /* length of the model (# nodes)          */
  float **t;                    /* transition prob's. t[(0),1..M-1][0..6] */
  float **mat;                  /* match emissions.  mat[1..M][0..K-1]    */ 
  float **ins;                  /* insert emissions. ins[1..M-1][0..K-1]  */
  /*::cexcerpt::plan7_core::end::*/

  /* Annotation. Everything but <name> is optional. Flags are set when
   * optional values are set.
   */
  char  *name;                  /* name of the model                     (mandatory) */
  char  *acc;			/* accession number of model (Pfam)      (p7_ACC)    */
  char  *desc;                  /* brief (1-line) description of model   (p7_DESC)   */ 
  char  *rf;                    /* reference line from alignment 1..M    (p7_RF)     */
  char  *cs;                    /* consensus structure line      1..M    (p7_CS)     */ 
  char  *ca;			/* consensus accessibility line  1..M    (p7_CA)     */
  char  *comlog;		/* command line(s) that built model      (mandatory) */
  int    nseq;			/* number of training sequences          (mandatory) */
  char  *ctime;			/* creation date                         (mandatory) */
  int   *map;			/* map of alignment cols onto model 1..M (p7_MAP)    */
  int    checksum;              /* checksum of training sequences        (mandatory) */

  /* Pfam-specific score cutoffs.
   * 
   * ga1, ga2 are valid if PLAN7_GA is set in flags.
   * tc1, tc2 are valid if PLAN7_TC is set in flags.
   * nc1, nc2 are valid if PLAN7_NC is set in flags.
   */
  float  ga1, ga2;	/* per-seq/per-domain gathering thresholds (bits) (p7_GA) */
  float  tc1, tc2;	/* per-seq/per-domain trusted cutoff (bits)       (p7_TC) */
  float  nc1, nc2;	/* per-seq/per-domain noise cutoff (bits)         (p7_NC) */

  int flags;
} P7_HMM;


/* Flag codes for hmm->flags.
 * Flags marked with ! may not be changed nor used for other meanings;
 * such flags were stored in old HMM files, and we must preserve their
 * meaning to preserve reverse compatibility.
 */
#define p7_HASBITS (1<<0)    /* obsolete (was: model has log-odds scores)       !*/
#define p7_DESC    (1<<1)    /* description exists                              !*/
#define p7_RF      (1<<2)    /* #RF annotation available                        !*/
#define p7_CS      (1<<3)    /* #CS annotation available                        !*/
#define p7_XRAY    (1<<4)    /* obsolete (was: structural data available)       !*/
#define p7_HASPROB (1<<5)    /* model has probabilities                         !*/
#define p7_HASDNA  (1<<6)    /* obsolete (was: protein HMM->DNA seq params set) !*/
#define p7_STATS   (1<<7)    /* obsolete (was: model has EVD stats calibrated)  !*/
#define p7_MAP     (1<<8)    /* alignment map is available                      !*/
#define p7_ACC     (1<<9)    /* accession number is available                   !*/
#define p7_GA      (1<<10)   /* gathering thresholds available                  !*/
#define p7_TC      (1<<11)   /* trusted cutoffs available                       !*/
#define p7_NC      (1<<12)   /* noise cutoffs available                         !*/
#define p7_CA      (1<<13)   /* surface accessibilities available               !*/


/* Indices for Plan7 main model state transitions.
 * Used for indexing hmm->t[k][]
 * mnemonic: Transition from Match to Match = TMM
 */
#define p7_TMM  0
#define p7_TMI  1
#define p7_TMD  2
#define p7_TIM  3
#define p7_TII  4
#define p7_TDM  5
#define p7_TDD  6 

/* Plan 7 model state types (esp. used in P7_TRACE structure)
 */
#define p7_BOGUS 0
#define p7_STM   1
#define p7_STD   2
#define p7_STI   3
#define p7_STS   4
#define p7_STN   5
#define p7_STB   6
#define p7_STE   7
#define p7_STC   8
#define p7_STT   9
#define p7_STJ   10     
#define p7_STX   11 	/* missing data: used esp. for local entry/exits */


extern P7_HMM *p7_hmm_Create(int M, ESL_ALPHABET *abc);
extern P7_HMM *p7_hmm_CreateShell(void);
extern int     p7_hmm_CreateBody(P7_HMM *hmm, int M, ESL_ALPHABET *abc);
extern void    p7_hmm_Destroy(P7_HMM *hmm);
extern int     p7_hmm_ZeroCounts(P7_HMM *hmm);
extern int     p7_hmm_Dump(FILE *fp, P7_HMM *hmm);
extern char   *p7_hmm_DescribeStatetype(char st);

extern int     p7_hmm_SetName(P7_HMM *hmm, char *name);
extern int     p7_hmm_SetAccession(P7_HMM *hmm, char *acc);
extern int     p7_hmm_SetDescription(P7_HMM *hmm, char *desc);
extern int     p7_hmm_AppendComlog(P7_HMM *hmm, int argc, char **argv);
extern int     p7_hmm_SetCtime(P7_HMM *hmm);

extern int     p7_hmm_Rescale(P7_HMM *hmm, float scale);
extern int     p7_hmm_Renormalize(P7_HMM *hmm);



#endif /* P7_HMM_INCLUDED */
/************************************************************
 * @LICENSE@
 ************************************************************/
