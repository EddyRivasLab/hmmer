/* 
 * main() for profile HMM construction from a single sequence and a general scoring system (score matrix + gap penalties)
 *
 * This standalone binary creates probabilistic Smith/Waterman models in HMMER3-style
 * It could be a part of hmmbuild with appropiate building options
 * In short, it would have to bypass typical profile building steps
 * and simply use a scoring matrix to parameterize the model (see below)
 *
 * SRE, Fri Mar 23 07:54:02 2007 [Janelia] [Decembrists, Picaresque]
 * SC, Wed Sep 17 17:52:11 EDT 2008 [Janelia, heavily modified] 
 * SVN $Id: seq2hmm.c 2579 2008-09-17 22:18:06Z castellanos $
 *
 * Contents
 *          1- main
 *          2- Public functions  (exposed  API)
 *          3- Private functions (shielded API)
 *
 * Development:
 *
 * gcc -o seq2hmm seq2hmm.c -g -Wall -Wextra -Wstrict-prototypes -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -ansi -std=c99 -pedantic -O3 -L../easel -I../easel -L. -I. -lhmmer -leasel -lm
 *
 * Production:
 *
 * gcc -o seq2hmm seq2hmm.c -g -Wall -O3 -L../easel -I../easel -L. -I. -lhmmer -leasel -lm 
 * 
 * USE: ./seq2hmm -m <matrix> -q gap_open_prob -r gap_extension_prob <hmmfile> <seqfile>
 *
 */


#include "p7_config.h"        /* HMMER hardcoded and system-dependent configuration
                               * Created by the ./configure script from p7_config.h.in
                               * DO NOT EDIT. MUST be first include in file.
                               */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* EASEL library */

#include "easel.h"            /* core easel functionality eg. error handling           */ 
#include "esl_getopts.h"      /* command-line processing                               */
#include "esl_vectorops.h"    /* vector operations eg. double to float vector          */
#include "esl_random.h"       /* random seeds                                          */
#include "esl_alphabet.h"     /* sequence alphabets eg. amino acids                    */ 
#include "esl_sq.h"           /* readseq-like function to format sequences             */
#include "esl_sqio.h"         /* input/output functions                                */
#include "esl_dmatrix.h"      /* algebra operations in matrix of doubles               */
#include "esl_scorematrix.h"  /* score matrix manipulation eg. joint prob. from scores */

/* The global variables, structures and defines in HMMER */

#include "hmmer.h"

/* Array of application defined command-line options */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-m",        eslARG_INFILE,  NULL, NULL, NULL,      NULL,  NULL, NULL, "use substitution score matrix file from <f>",    0 },
  { "-q",        eslARG_REAL,   "0.01", NULL, "0<=x<0.5",NULL,  NULL, NULL, "gap open probability",                          0 },
  { "-r",        eslARG_REAL,   "0.5", NULL, "0<=x<1",  NULL,  NULL, NULL, "gap extend probability",                         0 },
/* Control of E-value calibration */
  { "--Es",      eslARG_INT,    NULL,  NULL,"n>0",       NULL,    NULL,      NULL, "set random number seed to <n>",                     6},
  { "--EvL",     eslARG_INT,    "100", NULL,"n>0",       NULL,    NULL,      NULL, "length of sequences for Viterbi Gumbel mu fit",     6},   
  { "--EvN",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "number of sequences for Viterbi Gumbel mu fit",     6},   
  { "--EfL",     eslARG_INT,    "100", NULL,"n>0",       NULL,    NULL,      NULL, "length of sequences for Forward exp tail mu fit",   6},   
  { "--EfN",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "number of sequences for Forward exp tail mu fit",   6},   
  { "--Eft",     eslARG_REAL,  "0.04", NULL,"0<x<1",     NULL,    NULL,      NULL, "tail mass for Forward exponential tail mu fit",     6},  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes */

struct cfg_s {

  FILE         *ofp;		/* output file (default is stdout) */
  P7_BG	       *bg;		/* null model       */
  ESL_ALPHABET *abc;		/* digital alphabet */
  int          be_verbose; 	/* standard verbose output, as opposed to one-line-per-HMM summary */

};

static char usage[]  = "[-options] <HMM file> <FASTA file>"; /* these files are not options but arguments */
static char banner[] = "Create single sequence HMM for probabilistic Smith-Waterman in HMMER3";

/* PROTOTYPES */

static int
emissionate(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, char *mxfile, ESL_DMATRIX **Q);

int 
p7_Seqmodel(ESL_ALPHABET *abc, ESL_DSQ *dsq, int M, char *name,
		       ESL_DMATRIX *P, P7_BG *bg, const double popen, const double pextend,
		       P7_HMM **ret_hmm);
static int
calibrate(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, P7_HMM *hmm);


/***********
 * 1. MAIN *
 **********/

int main(int argc, char **argv)
{

  /* Parse command-line */

  ESL_GETOPTS    *go    = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);

  const double          popen   = esl_opt_GetReal  (go, "-q"); /* gap open probability  */
  const double          pextend = esl_opt_GetReal  (go, "-r"); /* gap open extenion     */
  char            *mxfile = esl_opt_GetString(go, "-m"); /* ptr to scoring matrix */

  /* Initialize the application config structure */
     
  struct cfg_s     cfg;

  cfg.ofp        = stdout;                                /* info stream                             */ 
  cfg.be_verbose = TRUE;                                  /* be verbose                              */   
  cfg.abc        = esl_alphabet_Create(eslAMINO);         /* set alphabet to amino acids (digitized) */    
  cfg.bg         = p7_bg_Create(cfg.abc);                 /* ptr to the background model (swissprot) */

  /* Files and File ptrs */

  char           *hmmfile = esl_opt_GetArg(go, 1);         /* name of output hmm file                          */
  char           *qfile   = esl_opt_GetArg(go, 2);         /* name of input query sequence file (FASTA format) */ 
  
  ESL_SQFILE     *qfp   = NULL;                            /* file ptr to query sequence file */ 							   
  FILE           *hmmfp = NULL;                            /* file ptr to output hmm file     */  

  /* Digitized sequence */

  ESL_SQ         *qsq   = esl_sq_CreateDigital(cfg.abc);   /* query sequence in digital form */

  /* Others */
 
  ESL_DMATRIX     *Q    = NULL;                             /* structure containing a double **mx matrix,  mx[i][j] is i'th row, j'th col */  
  P7_HMM          *hmm  = NULL;                             /* ptr to our HMM               */

  char            errbuf[eslERRBUFSIZE];                    /* message buffer */
  int             status;                                   /* program status */

  /* Open the query sequence file (FASTA format only) */
  /* Use esl_sqfile_GuessFileFormat() to guess format */

  status = esl_sqfile_Open(qfile, eslSQFILE_FASTA, NULL, &qfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file %s.", qfile);                    /* prints to stderr */
  else if (status == eslEFORMAT)   esl_fatal("Format of %s unrecognized.", qfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        esl_fatal("Open of %s failed, code %d.", qfile, status);

  /* Open the output HMM file */

  if ((hmmfp = fopen(hmmfile, "w")) == NULL) esl_fatal("Failed to open output HMM file %s", hmmfile); /* it will contain as many HMMs as input sequences */

  /* For each sequence, build a model and save it */

  while ((status = esl_sqio_Read(qfp, qsq)) == eslOK)    /* read sequence by sequence */
    {

      /* Compute emission probabilities from score matrix */

      if (emissionate(go, &cfg, errbuf, mxfile, &Q) != eslOK) esl_fatal("Computation of emission probabilities failed: %s\n", errbuf);

      /* Create HMM from single sequence (sets emission and transition probabilities into the model) */
 
      p7_Seqmodel(cfg.abc, qsq->dsq, qsq->n, qsq->name, Q, cfg.bg, popen, pextend, &hmm);

      
      /* Calibrate HMM */

      if (calibrate(go, &cfg, errbuf, hmm) != eslOK) esl_fatal("HMM calibration failed: %s\n", errbuf);

      /* Validate HMM */

      if (p7_hmm_Validate(hmm, errbuf, 1e-5) != eslOK) esl_fatal("HMM validation failed: %s\n", errbuf);

      /* Save HMM */

      if (p7_hmmfile_Write(hmmfp, hmm) != eslOK) esl_fatal("HMM save failed");
      { 
	p7_hmm_Destroy(hmm); /* free hmm */
      }
    }

  /* Check status after reading input sequence file */

  if (status != eslEOF) esl_fatal("Parse failed, line %d, query file %s:\n%s", 
				  qfp->linenumber, qfp->filename, qfp->errbuf);
  
  /* Go clean your room */

  esl_dmatrix_Destroy(Q);          /* free emissions matrix */
  esl_sq_Destroy(qsq);             /* free query sequence */ 
  esl_sqfile_Close(qfp);           /* close input query sequence file*/
  fclose(hmmfp);                   /* close hmm output file */
  esl_alphabet_Destroy(cfg.abc);   /* free alphabet */
  esl_getopts_Destroy(go);         /* free options  */  
  
  return 0;

} /* END main*/


/*************************************
 * 2. Public functions (exposed API) *
 *************************************/

/************************************************************************
 * Usage      : p7_Seqmodel()                              
 * Purpose    : Make a profile HMM from a single sequence  
                for probabilistic Smith/Waterman alignment 
 * Paradigm   : Imperative (procedural)                    
 * Type       : Function                                   
 * Scope      : Public                                     
 * Returns    : <eslOK> in succes                          
                status in error                            
 * Parameters : go    : options    
                abc      : alphabet
                dsq      : digital query sequence
                M        : sequence length (number of match states) 
                qsq->name: model name 
                Q        : emission probabilities
                bg->f    : background probabilities (from swissprot 50.8)
                popen    : tmi = tmd
                pextend  : tii = tdd
                hmm     : pointer of a pointer to an HMM structure
 * Throws     : <eslEMEM> on allocation error, and <*ret_hmm> is <NULL>    
 * Comments   : The query is digital sequence <dsq> of length <M>
                residues in alphabet <abc>, named <name>.             
                The scoring system is given by <Q>, <f>, <popen>, and
                <pextend>. <Q> is a $K \times K$ matrix giving
                conditional residue probabilities $P(a \mid b)}$; these
                are typically obtained by reverse engineering a score
                matrix like BLOSUM62. <f> is a vector of $K$ background
                frequencies $p_a$. <popen> and <pextend> are the
                probabilities assigned to gap-open ($t_{MI}$ and
                $t_{MD}$) and gap-extend ($t_{II}$ and $t_{DD}$)
                transitions.
 * See Also   : N/A                                       
 *********************************************************/

int p7_Seqmodel(ESL_ALPHABET *abc, ESL_DSQ *dsq, int M, char *name,
	    ESL_DMATRIX *Q, P7_BG *bg, const double popen, const double pextend,
	    P7_HMM **ret_hmm)
{
  int     status;

  P7_HMM *hmm    = NULL;                                           /* pointer to a plan7 HMM */

  char   *logmsg = "[HMM created from a single query sequence]";

  int     k;

  /* Allocate Plan 7 HMM structure (from p7_hmm.c)
   * hmm->M   set by p7_hmm_CreateBody when called from p7_hmm_Create
   * hmm->abc set by p7_hmm_CreateBody when called from p7_hmm_Create
   */

  if ((hmm = p7_hmm_Create(M, abc)) == NULL) { status = eslEMEM; goto ERROR; } /* returns error if allocation fails */
  
  /* Set emission and transition probabilities */

  for (k = 0; k <= M; k++)
    {
      /* Use rows of P matrix as source of match emission vectors (hmm->mat[1..M][0..K-1]) */

      if (k > 0) esl_vec_D2F(Q->mx[(int) dsq[k]], abc->K, hmm->mat[k]); /* copy a vector of doubles to a vector of floats */

      /* Set inserts to background for now. 
       * This can be improved using real
       * insertion frequencies. However, for
       * a fair algorithmic comparison with the
       * native S/W, I will stick for now with
       * the background frequencies.
       *
       * Insert emission probabilities derive for the most part 
       * (high alpha) from the slightly polar priors derived from 
       * Pfam 1.0. See p7_prior.c for details.
       *
       * Sets hmm->ins[1..M][0..K-1]
       */

      esl_vec_FCopy(bg->f, abc->K, hmm->ins[k]);                        /* copy background frequencies to insert states */

      /* Set transition probabilities */

      hmm->t[k][p7H_MM] = 1.0 - 2 * popen;
      hmm->t[k][p7H_MI] = popen;
      hmm->t[k][p7H_MD] = popen;
      hmm->t[k][p7H_IM] = 1.0 - pextend;
      hmm->t[k][p7H_II] = pextend;
      hmm->t[k][p7H_DM] = 1.0 - pextend;
      hmm->t[k][p7H_DD] = pextend;
    }

  /* Deal w/ special stuff at node M, overwriting a little of what we just did */

  hmm->t[M][p7H_MM] = 1.0 - popen;
  hmm->t[M][p7H_MD] = 0.;
  hmm->t[M][p7H_DM] = 1.0;
  hmm->t[M][p7H_DD] = 0.;
  
  /* Adds some annotation to the HMM model from p7_hmm.c */

  p7_hmm_SetName(hmm, name);                           /* set HMM name */

  p7_hmm_SetDescription(hmm, "Single-sequence model for probabilistic S/W"); /* one-line description */ 

  p7_hmm_AppendComlog(hmm, 1, &logmsg);                /* command-line used to build the model */

  hmm->nseq     = 0;                                   /* number of training sequences 
                                                        * NOT RELEVANT - emission prob.
                                                        * are obtained from scoring matrix, not counts
                                                        */

  hmm->eff_nseq = 0;                                   /* effective number of sequences
                                                        * NOT RELEVANT - emission prob.
                                                        * are obtained from scoring matrix, not 
							* weighted counts and priors
                                                        */

  p7_hmm_SetCtime(hmm);                                /* construction time */

  hmm->checksum = 0;                                   /* checksum of training sequence */

  /* Pointer to our brand new HMM */

  *ret_hmm = hmm;

  /* Return success */

  return eslOK;
  
  /* Errors */

 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;   /* status is eslEMEM */

} /* p7_seqmodel */

/***************************************
 * 2. Private functions (shielded API) *
 ***************************************/

/****************************************************************************
 * Usage      : emissionate()                                               *
 * Purpose    : Sets the emission probabilities of the model                *
 * Paradigm   : Imperative (procedural)                                     *
 * Type       : Function                                                    *
 * Scope      : Private                                                     *
 * Returns    : eslOK in succes                                             *
                status in error                                             *
 * Parameters : go, ptr to options structure                                *
                cfg_s, ptr to app configuration structure                   *               
                errbuf, ptr to error buffer                                 *
                mxfile, score matrix file                                   *
                Q, probability matrix                                       * 
 * Throws     : No exceptions                                               *
 * Comments   : Reverse engineer a scoring matrix to obtain conditional     *
                prob's that we'll use for the single-seq query HMM.         * 
                Because score mx is symmetric, we can set up                *
                P[a][b] = P(b | a), so we can use the matrix rows as HMM    *
                match emission vectors. This means dividing the joint probs *
                through by f_a.                                             *  
 * See Also   : N/A                                                         *
 ****************************************************************************/

static int
emissionate(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, char *mxfile, ESL_DMATRIX **Q)
{

  ESL_SCOREMATRIX *S    = esl_scorematrix_Create(cfg->abc);  /* matrix object */
   
  double          *fa   = NULL;                              /* vector of query background probabilities p(a) */
  double          *fb   = NULL;                              /* vector of target background probabilities p(b) */

  double          slambda;                                   /* score matrix lambda */ 
  int             a,b;

  /* Open scoring matrix */

  if (mxfile == NULL)  esl_scorematrix_SetBLOSUM62(S); /* BLOSUM62 is default */

  else {

    ESL_FILEPARSER *efp = NULL;

    if ( esl_fileparser_Open(mxfile, &efp) != eslOK) esl_fatal("failed to open score file %s",  mxfile);

    if ( esl_sco_Read(efp, cfg->abc, &S) != eslOK) esl_fatal("failed to read matrix from %s", mxfile);
    {
      esl_fileparser_Close(efp);  
    }
}

  /* Validate matrix simmetry */

  if (! esl_scorematrix_IsSymmetric(S)) esl_fatal("Score matrix isn't symmetric");

  /* Calculate the probabilistic basis of a score matrix
   *
   * ATTENTION: THIS ALGORITHM MAY FAIL FOR MATRICES BUILT USING LARGE LAMBDAS
   *            I MAY WANT TO CHECK esl_sco_ProbifyGivenBG()
   *
   * I need to explore whether I should be reverse engineering the matrix  
   * using the bg distribution in hmmer3 instead of the implicit in the matrix 
   *
   */

  if ( esl_sco_Probify(S, Q, &fa, &fb, &slambda) != eslOK) esl_fatal("failed to determine lambda from %s",  mxfile);

  printf("score matrix lambda: %f\n", slambda);

  /* esl_sco_ProbifyGivenBG(S, cfg->bg, cfg->bg, &slambda, Q);
   *
   * cfg->bg->f[a] to calculate conditionals
   *
   */

  /* Calculate conditional (emission) probabilities
   * from the joint p(i,j) probabilities now in *Q
   */

  for (a = 0; a < cfg->abc->K; a++)

    for (b = 0; b < cfg->abc->K; b++)
      (*Q)->mx[a][b] /= fa[a];	/* Q->mx[a][b] is now P(b | a) */


  /* Go clean your room */

  esl_scorematrix_Destroy(S);
  free(fa);
  free(fb);

  return eslOK;

} /* END emissionate */

/**********************************************************
 * Usage      : calibrate()                               *
 * Purpose    : Sets the E-value parameters of the model  *
 * Paradigm   : Imperative (procedural)                   *
 * Type       : Function                                  *
 * Scope      : Private                                   *
 * Returns    : eslOK in succes                           *
                status in error                           *
 * Parameters : go, ptr to options structure              *
                cfg_s, ptr to app configuration structure *               
                errbuf, ptr to error buffer               *
                hmm, ptr to out hmm in probability form   *
 * Throws     : No exceptions                             *
 * Comments   : Original function in hmmbuild.c           * 
                It is private function there, it          *
                should be moved to the exposed api        *
                Meantime, I include my private copy       *
 * See Also   : N/A                                       *
 **********************************************************/

static int
calibrate(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, P7_HMM *hmm)
{
  ESL_RANDOMNESS *r  = NULL;
  P7_PROFILE     *gm = NULL;
  int             vL = 100;	/* length of random seqs for Viterbi mu sim */
  int             vN = 200;	/* number of random seqs for Viterbi mu sim */
  int             fL = 100;	/* length of random seqs for Forward mu sim */
  int             fN = 200;	/* number of random seqs for Forward mu sim */
  double          ft = 0.04;	/* tail mass for Forward mu sim             */
  double lambda, mu, tau;
  double info, re;
  int    status;

  if (cfg->be_verbose) { fprintf(cfg->ofp, "Calibrating\n");    fflush(cfg->ofp); }

  if (esl_opt_IsDefault(go, "--Es"))  r = esl_randomness_CreateTimeseeded();                      /* default is null (no option in command-line) */
  else                                r = esl_randomness_Create(esl_opt_GetInteger(go, "--Es"));
  if (r == NULL) ESL_XFAIL(eslEMEM, errbuf, "failed to create random number generator");

  /* These functions are defined in evalues.c */

  if ((gm     = p7_profile_Create(hmm->M, cfg->abc))                  == NULL) ESL_XFAIL(eslEMEM, errbuf, "failed to allocate profile");
  if ((status = p7_ProfileConfig(hmm, cfg->bg, gm, vL, p7_LOCAL))    != eslOK) ESL_XFAIL(status,  errbuf, "failed to configure profile");
  if ((status = p7_Lambda(hmm, cfg->bg, &lambda))                    != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine lambda");
  if ((status = p7_Mu    (r, gm, cfg->bg, vL, vN, lambda, &mu))      != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine mu");
  if ((status = p7_Tau   (r, gm, cfg->bg, fL, fN, lambda, ft, &tau)) != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine tau");

  hmm->evparam[p7_LAMBDA] = lambda;
  hmm->evparam[p7_MU]     = mu;
  hmm->evparam[p7_TAU]    = tau;
  hmm->flags             |= p7H_STATS;   /* STATS is a  (1<<7) left shift bit operator: 1 * 2^7 */

  /* Relative entropy of the model
   * Mean Information content per match state emission distribution
   * Mean relative entropy per match state
   *
   */

  info = p7_MeanMatchInfo(hmm,cfg->bg);

  re = p7_MeanMatchRelativeEntropy(hmm,cfg->bg);

  printf("Information content: %f ", info);

  printf("Relative entropy: %f ", re);

  /* Print model/calibration info */

  printf("Lambda: %f ", lambda);
  printf("Mu: %f ", mu);
  printf("Tau: %f\n", tau);

  /* Go clean your room */

  p7_profile_Destroy(gm);
  esl_randomness_Destroy(r);
  if (cfg->be_verbose) fprintf(cfg->ofp, "done.\n\n");
  return eslOK;

 ERROR:
  esl_randomness_Destroy(r);
  p7_profile_Destroy(gm);
  if (cfg->be_verbose) fprintf(cfg->ofp, "FAILED.\n");
  return status;

} /* END calibrate */

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
