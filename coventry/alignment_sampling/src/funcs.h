/* funcs.h 
 * Declarations of external functions in HMMER.
 *
 * SVN $Id$
 */            
#ifndef FUNCSH_INCLUDED
#define FUNCSH_INCLUDED

#include "config.h"
#include "squid.h"
#include "msa.h"

#include "plan7.h"
#include "structs.h"

#include <easel.h>
#include <esl_random.h>		/* ESL_RANDOMNESS */

/* alphabet.c
 * Configuration of global alphabet information
 */
extern void           DetermineAlphabet(char **rseqs, int  nseq);
extern void           SetAlphabet(int type);
extern unsigned char  SymbolIndex(char sym);
extern unsigned char *DigitizeSequence(char *seq, int L);
extern char          *DedigitizeSequence(unsigned char *dsq, int L);
extern void           DigitizeAlignment(MSA *msa, unsigned char ***ret_dsqs);
extern void           P7CountSymbol(float *counters, unsigned char sym, float wt);
extern void           DefaultGeneticCode(int *aacode);
extern void           DefaultCodonBias(float *codebias);

/* from calibration.c
 * Determination of parameters needed for E-values.
 */
extern void P7CalibrateV(ESL_RANDOMNESS *r, struct plan7_s *hmm, double *fq, int N, int L,
			 float *ret_mu, float *ret_kappa);

/* from core_algorithms.c
 * Clean research/demonstration versions of basic algorithms.
 */
extern struct dpmatrix_s *CreatePlan7Matrix(int N, int M, int padN, int padM);
extern void   ResizePlan7Matrix(struct dpmatrix_s *mx, int N, int M, 
				int ***xmx, int ***mmx, int ***imx, int ***dmx);
extern struct dpmatrix_s *AllocPlan7Matrix(int rows, int M, 
					   int ***xmx, int ***mmx, int ***imx, int ***dmx);
extern struct dpshadow_s *AllocShadowMatrix(int rows, int M, char ***xtb, 
					    char ***mtb, char ***itb, char ***dtb);
extern void  FreePlan7Matrix(struct dpmatrix_s *mx);
extern void  FreeShadowMatrix(struct dpshadow_s *tb);
extern int   P7ViterbiSpaceOK(int L, int M, struct dpmatrix_s *mx);
extern int   P7ViterbiSize(int L, int M);
extern int   P7SmallViterbiSize(int L, int M);
extern int   P7WeeViterbiSize(int L, int M);
extern float P7Forward(unsigned char *dsq, int L, struct plan7_s *hmm, 
			  struct dpmatrix_s **ret_mx);
extern float P7Viterbi(unsigned char *dsq, int L, struct plan7_s *hmm, struct dpmatrix_s *mx,
			  struct p7trace_s **ret_tr);
extern float P7ViterbiNoTrace(unsigned char *dsq, int L, struct plan7_s *hmm,
			      struct dpmatrix_s *mx);
extern void  P7ViterbiTrace(struct plan7_s *hmm, unsigned char *dsq, int L,
			   struct dpmatrix_s *mx, struct p7trace_s **ret_tr);
extern float P7SmallViterbi(unsigned char *dsq, int L, struct plan7_s *hmm, 
			    struct dpmatrix_s *mx, struct p7trace_s **ret_tr);
extern float P7ParsingViterbi(unsigned char *dsq, int L, struct plan7_s *hmm, 
			      struct p7trace_s **ret_tr);
extern float P7WeeViterbi(unsigned char *dsq, int L, struct plan7_s *hmm, 
			  struct p7trace_s **ret_tr);
extern float Plan7ESTViterbi(unsigned char *dsq, int L, struct plan7_s *hmm, 
			     struct dpmatrix_s **ret_mx);
extern struct p7trace_s *P7ViterbiAlignAlignment(MSA *msa, struct plan7_s *hmm);
extern struct p7trace_s *ShadowTrace(struct dpshadow_s *tb, struct plan7_s *hmm, int L);
extern float  PostprocessSignificantHit(struct tophit_s *ghit, struct tophit_s *dhit,
					struct p7trace_s   *tr, struct plan7_s *hmm, unsigned char *dsq, 
					int L, char *seqname, char *seqacc, char *seqdesc, 
					int do_forward, float sc_override, int do_null2,
					struct threshold_s *thresh, int hmmpfam_mode);



/* from debug.c
 * Debugging output of various sorts.
 */
extern char *Statetype(char st);
extern char *AlphabetType2String(int type);
extern void P7PrintTrace(FILE *fp, struct p7trace_s *tr, 
			 struct plan7_s *hmm, unsigned char *dsq);
extern void P7PrintPrior(FILE *fp, struct p7prior_s *pri);
extern int  TraceCompare(struct p7trace_s *t1, struct p7trace_s *t2);
extern int  TraceVerify(struct p7trace_s *tr, int M, int N);

/* 
 * from display.c
 * Ian Holmes' functions for displaying HMMER2 data structures, especially
 * for posterior probabilities in alignments.
 */
extern void DisplayPlan7Matrix(unsigned char *dsq, int L, struct plan7_s *hmm,
			       struct dpmatrix_s *mx);
extern void DisplayPlan7Posteriors(int L, struct plan7_s *hmm,
				   struct dpmatrix_s *forward, struct dpmatrix_s *backward,
				   struct p7trace_s *viterbi, struct p7trace_s *optacc);
extern void DisplayPlan7PostAlign(int L, struct plan7_s *hmm,
				  struct dpmatrix_s *forward, struct dpmatrix_s *backward,
				  struct p7trace_s **alignment, int A);


/* from emit.c
 * Generation of sequences/traces from an HMM
 */
extern void EmitSequence(struct plan7_s *hmm, unsigned char **ret_dsq, int *ret_L, 
			 struct p7trace_s **ret_tr, ESL_RANDOMNESS *randomness);
extern void EmitConsensusSequence(struct plan7_s *hmm, char **ret_seq, unsigned char **ret_dsq,
				  int *ret_L, struct p7trace_s **ret_tr);
extern void StateOccupancy(struct plan7_s *hmm, float **ret_mp, float **ret_ip, float **ret_dp);


/* from emulation.c
 * Interfaces between HMMER and other software packages
 */
extern void WriteProfile(FILE *fp, struct plan7_s *hmm, int do_xsw);


/* from evolution.c
 * Phylogenetic extrapolation of profile HMMs.
 */
extern int EvolveOneTransitionVector(float *qs, float ts, int n, float *q0, float *qz, 
				     float t, float *q);


/* from histogram.c
 * accumulation of scores
 */
extern struct histogram_s *AllocHistogram(int min, int max, int lumpsize);
extern void FreeHistogram(struct histogram_s *h);
extern void UnfitHistogram(struct histogram_s *h);
extern void AddToHistogram(struct histogram_s *h, double sc);
extern void PrintASCIIHistogram(FILE *fp, struct histogram_s *h); 
extern void PrintXMGRHistogram(FILE *fp, struct histogram_s *h);
extern void PrintXMGRDistribution(FILE *fp, struct histogram_s *h);
extern void PrintXMGRRegressionLine(FILE *fp, struct histogram_s *h);
extern void EVDBasicFit(struct histogram_s *h);
extern int  ExtremeValueFitHistogram(struct histogram_s *h, int censor,
				     double high_hint);
extern void ExtremeValueSetHistogram(struct histogram_s *h, double mu, double lambda, 
				     double low, double high, int ndegrees);
extern int  GaussianFitHistogram(struct histogram_s *h, double high_hint);
extern void GaussianSetHistogram(struct histogram_s *h, double mean, double sd);

extern double EVDDensity(double x, double mu, double lambda);
extern double EVDDistribution(double x, double mu, double lambda);
extern double ExtremeValueP (double x, double mu, double lambda);
extern double ExtremeValueP2(double x, double mu, double lambda, int N);
extern double ExtremeValueE (double x, double mu, double lambda, int N);
extern double EVDrandom(double mu, double lambda);
extern int    EVDMaxLikelyFit(double *x, int *y, int n, 
			      double *ret_mu, double *ret_lambda);
extern int    EVDCensoredFit(double *x, int *y, int n, int z, double c, 
			     double *ret_mu, double *ret_lambda);
extern void   Lawless416(double *x, int *y, int n, double lambda, 
			 double *ret_f, double *ret_df);
extern void   Lawless422(double *x, int *y, int n, int z, double c,
			 double lambda, double *ret_f, double *ret_df);

/* from hmmio.c
 * Input/output (saving/reading) of models
 */
extern HMMFILE *HMMFileOpen(char *hmmfile, char *env);
extern int      HMMFileRead(HMMFILE *hmmfp, struct plan7_s **ret_hmm);
extern void     HMMFileClose(HMMFILE *hmmfp);
extern int      HMMFileFormat(HMMFILE *hmmfp);
extern void     HMMFileRewind(HMMFILE *hmmfp);
extern int      HMMFilePositionByName(HMMFILE *hmmfp, char *name);
extern int      HMMFilePositionByIndex(HMMFILE *hmmfp, int idx);
extern void     WriteAscHMM(FILE *fp, struct plan7_s *hmm);
extern void     WriteBinHMM(FILE *fp, struct plan7_s *hmm);

/* from infocontent.c
 * Evolving to specified information content
 */
extern void  AdjustAveInfoContent (struct plan7_s *hmm, float desired,
		char *matrixfile);
extern void  EvolveEmits (double *temp_emits, double *P, int L);
extern float CalculateBackgroundEntropy ();
extern float CalculateEmitsEntropy (double *emits, int L);
extern void  NormalizeEmits (double *temp_emits, int L);
extern void  PrintAveInfoContent (struct plan7_s *hmm);

/* island.c
 * Hwa's island method for fitting EVD.
 */
extern float IslandViterbi(unsigned char *dsq, int L, struct plan7_s *hmm,
			   int **ret_isle_sc, int **ret_isle_len, int *ret_inum);


/* masks.c
 * Repetitive sequence masking.
 */
extern int   XNU(unsigned char *dsq, int len);
extern float TraceScoreCorrection(struct plan7_s *hmm, struct p7trace_s *tr, unsigned char *dsq);
extern float ForwardScoreCorrection(struct plan7_s *hmm, unsigned char *dsq, int length,
				    struct dpmatrix_s *mx, ESL_RANDOMNESS *randomness);


/* mathsupport.c
 * Much of this code deals with Dirichlet prior mathematics.
 */
extern int   Prob2Score(float p, float null);
extern int   LL2Score(float ll, float null);
extern float Score2Prob(int sc, float null);
extern float Scorify(int sc);
extern double PValue(struct plan7_s *hmm, double sc);
extern float LogSum(float p1, float p2);
extern int   ILogsum(int p1, int p2);
extern void  LogNorm(float *vec, int n);
extern float Logp_cvec(float *cvec, int n, float *alpha);
extern float P_PvecGivenDirichlet(float *p, int n, float *alpha);

/* from matrices.c
 * for matrix manipulation
 */
extern void 	ReadAAMatrices(double **ret_Sij, double **ret_pi,
		char *matrixfile, int environ, int L);
extern void	AssignWagMatrix(double **Sij, double **pi);
extern void 	ReadMatrices (double **ret_Sij, double **ret_pi,
		char *matrixfile, int environ, int L);
extern void 	PrintMatrices(double *prob, int L, int environ);
extern void 	UnlogAndPrintMatrices(double *prob, int L, int environ);
extern void 	SymToRateMatrices(double *Qij, double *Sij, double *pi, int L,
		int environ);
extern void 	NormRateMatrices(double *Qij, double *pi, int L, int environ);
extern void  	AssignMatrixNotLog (double *Qij, int matrix_n, double time,
		double *Pij);
extern double 	*Cal_Id(int L);
extern double  *Cal_M_Exp(double *M, int L, double power);
extern void    Comp_M_Exp(double *M, int L, double power);
extern void    Comp_M_N_Prod(double *M, double *N, int L);
extern void    CheckSingleProb(double *psingle, int size);
extern void    CopyMatrix (double *copyQ, double *Q, int N);
extern int     Check_Accuracy(double *vec, int L);
extern void	LogifyMatrix(double *M, int L);

/* from misc.c
 * Miscellaneous functions with no home
 */
extern void  HMMERBanner(FILE *fp, char *banner);
extern char *Getword(FILE *fp, int type);
extern char *Getline(char *s, int n, FILE *fp);
extern int   SetAutocuts(struct threshold_s *thresh, struct plan7_s *hmm);

/* from modelconfig.c
 * Model configuration, from core model probabilities
 * to the full Plan7 score model
 */
extern void  P7Config(struct plan7_s *hmm, enum p7_algmode mode);
extern void  P7ReconfigLength (struct plan7_s *hmm, int L);
extern float P7EdgeCorrection (struct plan7_s *hmm, int L);
extern float P7ScoreCorrection(struct plan7_s *hmm, int L);


/* from modelmakers.c
 * Model construction algorithms
 */
extern void P7Handmodelmaker(MSA *msa, unsigned char **dsq, char *isfrag,
			     struct plan7_s **ret_hmm,
			     struct p7trace_s ***ret_tr);
extern void P7Fastmodelmaker(MSA *msa, unsigned char **dsq, char *isfrag,
			     float symfrac, struct plan7_s **ret_hmm,
			     struct p7trace_s ***ret_tr);

/* from plan7.c
 * Plan7 HMM structure support
 */
extern struct plan7_s *AllocPlan7(int M);
extern struct plan7_s *AllocPlan7Shell(void);
extern void AllocPlan7Body(struct plan7_s *hmm, int M);
extern void FreePlan7(struct plan7_s *hmm);
extern void ZeroPlan7(struct plan7_s *hmm);
extern void Plan7SetName(struct plan7_s *hmm, char *name);
extern void Plan7SetAccession(struct plan7_s *hmm, char *acc);
extern void Plan7SetDescription(struct plan7_s *hmm, char *desc);
extern void Plan7ComlogAppend(struct plan7_s *hmm, int argc, char **argv);
extern void Plan7SetCtime(struct plan7_s *hmm);
extern void Plan7SetNullModel(struct plan7_s *hmm, float null[MAXABET], float p1);
extern void Plan7Rescale(struct plan7_s *hmm, float scale);
extern void Plan7Renormalize(struct plan7_s *hmm);
extern void PrintPlan7Stats(FILE *fp, struct plan7_s *hmm, unsigned char **dsq, 
			    int nseq, struct p7trace_s **tr);
extern int  DegenerateSymbolScore(float *p, float *null, int ambig);
extern void Plan7_DumpScores(FILE *fp, struct plan7_s *hmm);
extern void Plan7_DumpCounts(FILE *fp, struct plan7_s *hmm);

/* 
 * from plan9.c
 * Backwards compatibility for the Plan 9 data structures of HMMER 1.x
 */
extern struct plan9_s *P9AllocHMM(int M);
extern void P9ZeroHMM(struct plan9_s *hmm);
extern int  P9FreeHMM(struct plan9_s *hmm);
extern void P9Renormalize(struct plan9_s *hmm);
extern void P9DefaultNullModel(float *null);
extern void Plan9toPlan7(struct plan9_s *hmm, struct plan7_s **ret_plan7);

/* 
 * from postprob.c
 * Functions for working with posterior probabilities within alignments
 */
extern float P7OptimalAccuracy(unsigned char *dsq, int L, struct plan7_s *hmm, struct p7trace_s **ret_tr);
extern float P7Backward(unsigned char *dsq, int L, struct plan7_s *hmm,	struct dpmatrix_s **ret_mx);
extern void  P7EmitterPosterior(int L, struct plan7_s *hmm, struct dpmatrix_s *forward,
				struct dpmatrix_s *backward, struct dpmatrix_s *mx);
extern float P7FillOptimalAccuracy(int L, int M, struct dpmatrix_s *posterior,
				   struct dpmatrix_s *mx, struct p7trace_s **ret_tr);
extern void  P7OptimalAccuracyTrace(int L, int M, struct dpmatrix_s *posterior,
				    struct dpmatrix_s *mx, struct p7trace_s **ret_tr);
extern char *PostalCode(int L, struct dpmatrix_s *mx, struct p7trace_s *tr);

/* from prior.c
 * Dirichlet priors
 */
extern struct p7prior_s *P7AllocPrior(void);
extern struct p7prior_s *P7LaplacePrior(void);
extern struct p7prior_s *P7DefaultPrior(void);
extern struct p7prior_s *P7ReadPrior(char *prifile);
extern void P7FreePrior(struct p7prior_s *pri);
extern void PAMPrior(char *pamfile, struct p7prior_s *pri, float pamwgt);
extern void P7DefaultNullModel(float *null, float *ret_p1);
extern void P7ReadNullModel(char *rndfile, float *null, float *ret_p1);
extern void P7PriorifyHMM(struct plan7_s *hmm, struct p7prior_s *pri);
extern void P7PriorifyTransitionVector(float *t, struct p7prior_s *prior, 
				       float tq[MAXDCHLET]);
extern void P7PriorifyEmissionVector(float *vec, struct p7prior_s *pri, 
				     int num, float eq[MAXDCHLET], 
				     float e[MAXDCHLET][MAXABET],
				     float *ret_mix);


#ifdef HMMER_PVM
/* from pvm.c
 * PVM Parallel Virtual Machine implementation
 */
extern void              PVMSpawnSlaves(char *slave, int **ret_tid, int *ret_nslaves);
extern void              PVMConfirmSlaves(int *slave_tid, int nslaves);
extern void              PVMCheckSlaves(int *slave_tid, int nslaves);
extern void              PVMKillSlaves(int *slave_tid, int nslaves);
extern int               PVMPackString(char *s);
extern char *            PVMUnpackString(void);
extern int               PVMPackTrace(struct p7trace_s *tr);
extern struct p7trace_s *PVMUnpackTrace(void);
extern int               PVMPackHMM(struct plan7_s *hmm);
extern struct plan7_s *  PVMUnpackHMM(void);
#endif /*HMMER_PVM*/

/* from threads.c
 * util for POSIX threads implementation
 */
extern int   ThreadNumber(void);

/* from tophits.c
 * Support for keeping/sorting top scoring hit/alignment lists
 */
extern struct tophit_s *AllocTophits(int lumpsize);
extern void   GrowTophits(struct tophit_s *h);
extern void   FreeTophits(struct tophit_s *h);
extern struct fancyali_s *AllocFancyAli(void);
extern void   FreeFancyAli(struct fancyali_s *ali);
extern void   RegisterHit(struct tophit_s *h, double sortkey, 
			  double pvalue, float score, 
			  double motherp, float mothersc,
			  char *name, char *acc, char *desc, 
			  int sqfrom, int sqto, int sqlen, 
			  int hmmfrom, int hmmto, int hmmlen, 
			  int domidx, int ndom, 
			  struct fancyali_s *ali);
extern void GetRankedHit(struct tophit_s *h, int rank, 
			 double *r_pvalue, float *r_score, 
			 double *r_motherp, float *r_mothersc,
			 char **r_name, char **r_acc, char **r_desc,
			 int *r_sqfrom, int *r_sqto, int *r_sqlen,
			 int *r_hmmfrom, int *r_hmmto, int *r_hmmlen,
			 int *r_domidx, int *r_ndom,
			 struct fancyali_s **r_ali);
extern int    TophitsMaxName(struct tophit_s *h);
extern void   FullSortTophits(struct tophit_s *h);
extern void   TophitsReport(struct tophit_s *h, double E, int nseq);

/* from trace.c
 * Support for traceback (state path) structure
 */
extern void  P7AllocTrace(int tlen, struct p7trace_s **ret_tr);
extern void  P7ReallocTrace(struct p7trace_s *tr, int tlen);
extern void  P7FreeTrace(struct p7trace_s *tr);
extern void  TraceSet(struct p7trace_s *tr, int tpos, char type, int idx, int pos);
extern struct p7trace_s **MergeTraceArrays(struct p7trace_s **t1, int n1, struct p7trace_s **t2, int n2);
extern void  P7ReverseTrace(struct p7trace_s *tr);
extern void  P7TraceCount(struct plan7_s *hmm, unsigned char *dsq, float wt, 
			  struct p7trace_s *tr);
extern float P7TraceScore(struct plan7_s *hmm, unsigned char *dsq, struct p7trace_s *tr);
extern MSA  *P7Traces2Alignment(unsigned char **dsq, SQINFO *sqinfo, float *wgt, 
				int nseq, int M, 
				struct p7trace_s **tr, int matchonly);
extern int  TransitionScoreLookup(struct plan7_s *hmm, char st1, 
				  int k1, char st2, int k2);
extern struct fancyali_s *CreateFancyAli(struct p7trace_s *tr, struct plan7_s *hmm,
					 unsigned char *dsq, char *name);
extern void PrintFancyAli(FILE *fp, struct fancyali_s *ali);
extern void TraceDecompose(struct p7trace_s *otr, struct p7trace_s ***ret_tr,
			   int *ret_ntr);
extern int  TraceDomainNumber(struct p7trace_s *tr);
extern int  Trace_GetAlignmentBounds(struct p7trace_s *tr, int which,
				     int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2,
				     int *ret_avlen);
extern int  P7TraceCountAnnotated(struct p7trace_s *tr);
extern struct p7trace_s *MasterTraceFromMap(int *map, int M, int alen);
extern void ImposeMasterTrace(char **aseq, int nseq, struct p7trace_s *mtr, 
			      struct p7trace_s ***ret_tr);

/* from sample.c
 */
extern void P7SampleAlignment(struct plan7_s *hmm, unsigned char *dsq, int N,
			      struct dpmatrix_s *mx, struct p7trace_s **ret_tr,
			      ESL_RANDOMNESS *randomness);

#endif /*FUNCSH_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/

