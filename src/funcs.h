/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the 
 *   GNU General Public License. See the files COPYING and 
 *   GNULICENSE for details.
 *    
 ************************************************************/

/* funcs.h 
 * RCS $Id$
 *
 * Declarations of external functions in HMMER.
 */            

#ifndef FUNCSH_INCLUDED
#define FUNCSH_INCLUDED

#include "config.h"
#include "structs.h"
#include "squid.h"

/* alphabet.c
 * Configuration of global alphabet information
 */
extern void  DetermineAlphabet(char **rseqs, int  nseq);
extern void  SetAlphabet(int type);
extern int   SymbolIndex(char sym);
extern char *DigitizeSequence(char *seq, int L);
extern char *DedigitizeSequence(char *dsq, int L);
extern void  DigitizeAlignment(char **aseqs, AINFO *ainfo, char ***ret_dsqs);
extern void  P7CountSymbol(float *counters, char sym, float wt);
extern void  DefaultGeneticCode(int *aacode);
extern void  DefaultCodonBias(float *codebias);

/* from core_algorithms.c
 * Clean research/demonstration versions of basic algorithms.
 */
extern struct dpmatrix_s *AllocPlan7Matrix(int rows, int M, int ***xmx, 
					   int ***mmx, int ***imx, int ***dmx);
extern struct dpshadow_s *AllocShadowMatrix(int rows, int M, char ***xtb, 
					    char ***mtb, char ***itb, char ***dtb);
extern void  FreePlan7Matrix(struct dpmatrix_s *mx);
extern void  FreeShadowMatrix(struct dpshadow_s *tb);
extern int   P7ViterbiSize(int L, int M);
extern int   P7SmallViterbiSize(int L, int M);
extern int   P7WeeViterbiSize(int L, int M);
extern float P7Forward(char *dsq, int L, struct plan7_s *hmm, 
			  struct dpmatrix_s **ret_mx);
extern float P7Viterbi(char *dsq, int L, struct plan7_s *hmm, 
			  struct p7trace_s **ret_tr);
extern void  P7ViterbiTrace(struct plan7_s *hmm, char *dsq, int L,
			   struct dpmatrix_s *mx, struct p7trace_s **ret_tr);
extern float P7SmallViterbi(char *dsq, int L, struct plan7_s *hmm, struct p7trace_s **ret_tr);
extern float P7ParsingViterbi(char *dsq, int L, struct plan7_s *hmm, 
			      struct p7trace_s **ret_tr);
extern float P7WeeViterbi(char *dsq, int L, struct plan7_s *hmm, 
			  struct p7trace_s **ret_tr);
extern float Plan7ESTViterbi(char *dsq, int L, struct plan7_s *hmm, 
			     struct dpmatrix_s **ret_mx);
extern struct p7trace_s *P7ViterbiAlignAlignment(char **aseq, AINFO *ainfo, struct plan7_s *hmm);
extern struct p7trace_s *ShadowTrace(struct dpshadow_s *tb, struct plan7_s *hmm, int L);



/* from debug.c
 * Debugging output of various sorts.
 */
extern char *Statetype(char st);
extern void P7PrintTrace(FILE *fp, struct p7trace_s *tr, 
			 struct plan7_s *hmm, char *dsq);
extern void P7PrintPrior(FILE *fp, struct p7prior_s *pri);
extern int  TraceCompare(struct p7trace_s *t1, struct p7trace_s *t2);
extern int  TraceVerify(struct p7trace_s *tr, int M, int N);

/* from emit.c
 * Generation of sequences/traces from an HMM
 */
extern void EmitSequence(struct plan7_s *hmm, char **ret_dsq, int *ret_L, struct p7trace_s **ret_tr);


/* from emulation.c
 * Interfaces between HMMER and other software packages
 */
extern void WriteProfile(FILE *fp, struct plan7_s *hmm, int do_xsw);


/* from histogram.c
 * accumulation of scores
 */
extern struct histogram_s *AllocHistogram(int min, int max, int lumpsize);
extern void FreeHistogram(struct histogram_s *h);
extern void UnfitHistogram(struct histogram_s *h);
extern void AddToHistogram(struct histogram_s *h, float sc);
extern void PrintASCIIHistogram(FILE *fp, struct histogram_s *h); 
extern void PrintXMGRHistogram(FILE *fp, struct histogram_s *h);
extern void PrintXMGRDistribution(FILE *fp, struct histogram_s *h);
extern void PrintXMGRRegressionLine(FILE *fp, struct histogram_s *h);
extern void EVDBasicFit(struct histogram_s *h);
extern int  ExtremeValueFitHistogram(struct histogram_s *h, int censor,
				     float high_hint);
extern void ExtremeValueSetHistogram(struct histogram_s *h, float mu, float lambda, 
				     float low, float high, int ndegrees);
extern int  GaussianFitHistogram(struct histogram_s *h, float high_hint);
extern void GaussianSetHistogram(struct histogram_s *h, float mean, float sd);
extern double EVDDensity(float x, float mu, float lambda);
extern double EVDDistribution(float x, float mu, float lambda);
extern double ExtremeValueP (float x, float mu, float lambda);
extern double ExtremeValueP2(float x, float mu, float lambda, int N);
extern double ExtremeValueE (float x, float mu, float lambda, int N);
extern float  EVDrandom(float mu, float lambda);
extern int    EVDMaxLikelyFit(float *x, int *y, int n, 
			      float *ret_mu, float *ret_lambda);
extern int    EVDCensoredFit(float *x, int *y, int n, int z, float c, 
			     float *ret_mu, float *ret_lambda);
extern void   Lawless416(float *x, int *y, int n, float lambda, 
			 float *ret_f, float *ret_df);
extern void   Lawless422(float *x, int *y, int n, int z, float c,
			 float lambda, float *ret_f, float *ret_df);

/* from hmmio.c
 * Input/output (saving/reading) of models
 */
extern HMMFILE *HMMFileOpen(char *hmmfile, char *env);
extern int      HMMFileRead(HMMFILE *hmmfp, struct plan7_s **ret_hmm);
extern void     HMMFileClose(HMMFILE *hmmfp);
extern int      HMMFileFormat(HMMFILE *hmmfp);
extern void     HMMFileRewind(HMMFILE *hmmfp);
extern void     HMMFilePositionByOffset(HMMFILE *hmmfp, long offset);
extern int      HMMFilePositionByName(HMMFILE *hmmfp, char *name);
extern int      HMMFilePositionByIndex(HMMFILE *hmmfp, int idx);
extern void     WriteAscHMM(FILE *fp, struct plan7_s *hmm);
extern void     WriteBinHMM(FILE *fp, struct plan7_s *hmm);

/* masks.c
 * Repetitive sequence masking.
 */
extern int   XNU(char *dsq, int len);
extern float TraceScoreCorrection(struct plan7_s *hmm, struct p7trace_s *tr, char *dsq);

/* mathsupport.c
 * Much of this code deals with Dirichlet prior mathematics. 
 */
extern int   Prob2Score(float p, float null);
extern float Score2Prob(int sc, float null);
extern float Scorify(int sc);
extern double PValue(struct plan7_s *hmm, float sc);
extern float LogSum(float p1, float p2);
extern int   ILogsum(int p1, int p2);
extern void  LogNorm(float *vec, int n);
extern float Logp_cvec(float *cvec, int n, float *alpha);
extern void  SampleDirichlet(float *alpha, int n, float *p);
extern float SampleGamma(float alpha);
extern void  SampleCountvector(float *p, int n, int c, float *cvec);
extern float P_PvecGivenDirichlet(float *p, int n, float *alpha);

/* from misc.c
 * Miscellaneous functions with no home
 */
extern void  Banner(FILE *fp, char *banner);
extern char *Getword(FILE *fp, int type); 
extern char *Getline(char *s, int n, FILE *fp);


/* from modelmakers.c
 * Model construction algorithms
 */
extern void P7Handmodelmaker(char **aseq, char **dsq, AINFO *ainfo,
			     struct plan7_s **ret_hmm,
			     struct p7trace_s ***ret_tr);
extern void P7Fastmodelmaker(char **aseq, char **dsq, AINFO *ainfo,
			     float maxgap, struct plan7_s **ret_hmm, 
			     struct p7trace_s ***ret_tr);
extern void P7Maxmodelmaker(char **aseqs, char **dsq, AINFO *ainfo,
			    float maxgap, struct p7prior_s *prior, 
			    float *null, float null_p1, float mpri, 
			    struct plan7_s **ret_hmm,
			    struct p7trace_s  ***ret_tr);

/* from plan7.c
 * Plan7 HMM structure support
 */
extern struct plan7_s *AllocPlan7(int M);
extern struct plan7_s *AllocPlan7Shell(void);
extern void AllocPlan7Body(struct plan7_s *hmm, int M);
extern void FreePlan7(struct plan7_s *hmm);
extern void ZeroPlan7(struct plan7_s *hmm);
extern void Plan7SetName(struct plan7_s *hmm, char *name);
extern void Plan7SetDescription(struct plan7_s *hmm, char *desc);
extern void Plan7ComlogAppend(struct plan7_s *hmm, int argc, char **argv);
extern void Plan7SetCtime(struct plan7_s *hmm);
extern void Plan7SetNullModel(struct plan7_s *hmm, float null[MAXABET], float p1);
extern void P7Logoddsify(struct plan7_s *hmm, int viterbi_mode);
extern void Plan7Renormalize(struct plan7_s *hmm);
extern void Plan7RenormalizeExits(struct plan7_s *hmm);
extern void Plan7NakedConfig(struct plan7_s *hmm);
extern void Plan7GlobalConfig(struct plan7_s *hmm);
extern void Plan7LSConfig(struct plan7_s *hmm);
extern void Plan7SWConfig(struct plan7_s *hmm, float pentry, float pexit);
extern void Plan7FSConfig(struct plan7_s *hmm, float pentry, float pexit); 
extern void PrintPlan7Stats(FILE *fp, struct plan7_s *hmm, char **dsq, 
			    int nseq, struct p7trace_s **tr);
extern int  DegenerateSymbolScore(float *p, float *null, int ambig);
extern void Plan9toPlan7(struct plan9_s *hmm, struct plan7_s **ret_plan7);

/* 
 * from plan9.c
 * Backwards compatibility for the Plan 9 data structures of HMMER 1.x
 */
extern struct plan9_s *P9AllocHMM(int M);
extern void P9ZeroHMM(struct plan9_s *hmm);
extern int  P9FreeHMM(struct plan9_s *hmm);
extern void P9Renormalize(struct plan9_s *hmm);
extern void P9DefaultNullModel(float *null);

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
extern void P7PriorifyTransitionVector(float *t, struct p7prior_s *prior);
extern void P7PriorifyEmissionVector(float *vec, struct p7prior_s *pri, 
				     int num, float eq[MAXDCHLET], 
				     float e[MAXDCHLET][MAXABET],
				     float *ret_mix);


#ifdef HMMER_PVM
/* from pvm.c
 * PVM Parallel Virtual Machine implementation
 */
extern void              PVMSpawnSlaves(char *slave, int **ret_tid, int *ret_nslaves);
extern void              PVMCheckSlaves(int *slave_tid, int nslaves);
extern void              PVMKillSlaves(int *slave_tid, int nslaves);
extern int               PVMPackString(char *s);
extern char *            PVMUnpackString(void);
extern int               PVMPackTrace(struct p7trace_s *tr);
extern struct p7trace_s *PVMUnpackTrace(void);
extern int               PVMPackHMM(struct plan7_s *hmm);
extern struct plan7_s *  PVMUnpackHMM(void);
#endif /*HMMER_PVM*/

#ifdef HMMER_THREADS
/* from threads.c
 * POSIX threads implementation
 */
extern int   ThreadNumber(void);
#endif /*HMMER_THREADS*/


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
			  char *name, char *desc, 
			  int sqfrom, int sqto, int sqlen, 
			  int hmmfrom, int hmmto, int hmmlen, 
			  int domidx, int ndom, 
			  struct fancyali_s *ali);
extern void GetRankedHit(struct tophit_s *h, int rank, 
			 double *r_pvalue, float *r_score, 
			 double *r_motherp, float *r_mothersc,
			 char **r_name, char **r_desc,
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
extern void  P7TraceCount(struct plan7_s *hmm, char *dsq, float wt, 
			  struct p7trace_s *tr);
extern float P7TraceScore(struct plan7_s *hmm, char *dsq, struct p7trace_s *tr);
extern void  P7Traces2Alignment(char **dsq, SQINFO *sqinfo, float *wgt, int nseq, int M, 
				struct p7trace_s **tr, int matchonly, 
				char ***ret_aseqs, AINFO *ainfo);
extern int  TransitionScoreLookup(struct plan7_s *hmm, char st1, 
				  int k1, char st2, int k2);
extern struct fancyali_s *CreateFancyAli(struct p7trace_s *tr, struct plan7_s *hmm,
					 char *dsq, char *name);
extern void PrintFancyAli(FILE *fp, struct fancyali_s *ali);
extern void TraceDecompose(struct p7trace_s *otr, struct p7trace_s ***ret_tr,
			   int *ret_ntr);
extern int  TraceDomainNumber(struct p7trace_s *tr);
extern void TraceSimpleBounds(struct p7trace_s *tr, int *ret_i1, int *ret_i2, 
			      int *ret_k1,  int *ret_k2);
extern struct p7trace_s *MasterTraceFromMap(int *map, int M, int alen);
extern void ImposeMasterTrace(char **aseq, int nseq, struct p7trace_s *mtr, 
			      struct p7trace_s ***ret_tr);


#endif /*FUNCSH_INCLUDED*/
