/* evohmmer - funtions to evolve parameters of an HMM
 *
 */
#ifndef EVOHMMER_INCLUDED
#define EVOHMMER_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_scorematrix.h"

#include "hmmer.h"

#include "e2.h"
#include "e1_rate.h"
#include "ratematrix.h"

#define NDTS  30
#define NDTL  1
#define NDTR  1
#define NDT   NDTS+NDTL+2*NDTR
#define DTMIN 0.01
#define DTMAX 5.50
#define DTPRE 0.05
#define DTPOS 5.00


enum emevol_e  { emNONE = 0, emBYRATE = 1, emBYSCMX = 2 };

typedef struct {
  FILE            *statfp;     // stats output stream 
  EMRATE          *emR;        // emissions rate
  ESL_SCOREMATRIX *S;
  EVOM             evomodel;   // the evomodel
  float            fixtime;    // 
  float            betainf;
  double           tol;
} HMMRATE;

typedef struct p7_rate_s {
  int      M;
  EVOM     evomodel;
  int      nrate;
  int      nbern;

  float    betainf;

  char     *name;               /* name of the model                     (mandatory)      */ /* String, \0-terminated   */
  char     *acc;	        /* accession number of model (Pfam)      (optional: NULL) */ /* String, \0-terminated   */
  char     *desc;               /* brief (1-line) description of model   (optional: NULL) */ /* String, \0-terminated   */
 
  /* emission probabilities */
  float   **pzero;              /* match emissions at t=zero     [1..M][0..K-1]                  */
  float   **pstar;              /* match emissions at t=star     [1..M][0..K-1]                  */
  float   **pinfy;              /* match emissions at t=infinity [1..M][0..K-1]                  */
  float   **ins;                /* insert emissions [1..M][0..K-1]                               */
 
  float           fixtime;
  int             ndt;
  float          *dtval;
  ESL_DMATRIX  ***Pdt;          /* conditional matrix at dt [1..M][0..K-1][0..K-1]               */
  float        ***pdt;          /* match emissions at dt [1..M][0..K-1]                          */

 
  /* rates per profile position */
  E1_RATE **e1R;                /* a rate structure for each profile position [(0),1..M]         */

  /* bookkeeping */
  const ESL_ALPHABET *abc_r;	/* reference to the alphabet: includes K, Kp, and sym order */
} P7_RATE;

struct entropy_param_s {
  float      **pref;
  float      **prob;
  int          M;       // prob[1], prob[M] are probability vectors [0..K-1]
  int          K;
  double       etarget;	// information content target, in bits 
};

extern P7_RATE *p7_RateCreateBare(const ESL_ALPHABET *abc, EVOM evomodel, float fixtime, float betainf);
extern P7_RATE *p7_RateCreate(int M, const ESL_ALPHABET *abc, EVOM evomodel, float fixtime, float betainf);
extern int      p7_RateCreateWithEmRate(int M, const ESL_ALPHABET *abc, const P7_BG *bg, HMMRATE *hmmrate, P7_RATE **ret_R, char *errbuf, int verbose);
extern double   p7_RateCompare(P7_RATE *R1, P7_RATE *R2, double tol);
extern int      p7_RateCopy(P7_RATE *R, P7_RATE *Rcopy);
extern P7_RATE *p7_RateClone(P7_RATE *R);
extern void     p7_RateDestroy(P7_RATE *R);
extern void     p7_RateDump(FILE *fp, P7_RATE *R);
extern int      p7_RateCalculate(const P7_HMM *hmm, const P7_BG *bg, HMMRATE *hmmrate, P7_RATE *R,      char *errbuf, int verbose);
extern int      p7_RateConstruct(const P7_HMM *hmm, const P7_BG *bg, HMMRATE *hmmrate, P7_RATE **ret_R, char *errbuf, int verbose);
extern int      p7_RateTransitions(const P7_HMM *hmm, P7_RATE *R, double betainf, double tol, char *errbuf, int verbose);
extern int      p7_RateValidate(P7_RATE *R, double tol, char *errbuf);
extern int      p7_EvolveFromRate(FILE *statfp, P7_HMM *hmm, const P7_RATE *R, const P7_BG *bg, double time, double tol, char *errbuf, int verbose);
extern int      p7_Evolve(P7_HMM *hmm, const P7_BG *bg, double time, HMMRATE *hmmrate, char *errbuf, int verbose);
extern int      p7_RateTest(const P7_HMM *hmm, P7_RATE *R, const P7_BG *bg, double tol, char *errbuf, int verbose);
extern int      p7_RestoreHMMmat(P7_HMM *hmm, P7_RATE *R, char *errbuf, int verbose);
extern int      p7_CalculatePzero(P7_RATE *R, const P7_HMM *hmm, char *errbuf, int verbose);
extern int      p7_CalculatePinfy(P7_RATE *R, const P7_HMM *hmm, const P7_BG *bg, char *errbuf, int verbose);
extern int      er_EntropyWeight(float **prob, int M, int K, float **pref, double etarget, double *ret_cut);
#endif /*EVOHMMER_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
