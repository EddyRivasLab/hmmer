/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1997 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the 
 *   GNU General Public License. See the files COPYING and 
 *   GNULICENSE for details.
 *    
 ************************************************************/

/* 
 * funcs.h 
 *
 * Declarations of external functions in HMMER.
 * 
 * RCS $Id$
 */            

#ifndef FUNCSH_INCLUDED
#define FUNCSH_INCLUDED

#ifdef BLAMM
#include <ncbi.h>
#include <dfa.h>
#endif /* BLAMM */

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
extern void  DigitizeAlignment(char **aseqs, AINFO *ainfo, char ***ret_dsqs);
extern void  P7CountSymbol(float *counters, char sym, float wt);
extern void  DefaultGeneticCode(int *aacode);
extern void  DefaultCodonBias(float *codebias);

/* 
 * from blast.c
 * BLAMM: Experimental BLAST/HMM hybrid code.
 */
#ifdef BLAMM
extern void NeighborhoodWords(struct hmm_struc *phmm, DFAPtr dfa, 
			      struct wordpool_s *pool, int target_len, float T,
			      FILE *ofp, int *ret_nwords);
extern struct wordpool_s *Init_HMMWordPool(void);
extern void               Free_HMMWordPool(struct wordpool_s *pool);
extern int  ForwardExtension(struct shmm_s *shmm, char *seq1, int L, 
			     int starti, int startk, int startsc, int blast_X,
			     float *ret_sc, int *ret_j);
#endif /*BLAMM*/

/* 
 * from camJul97.c
 * Island HMM experimental code
 */
extern void MakeStarHMM(struct plan7_s *hmm, float **mx, float *pq,
			int nq, struct p7prior_s *pri);
extern void ReadGJMMatrices(FILE *fp, float ***ret_mx, float **ret_pq, int *ret_nq);
extern void Joint2SubstitutionMatrix(float **jmx, float **smx, int N);
extern struct plan7_s *MakeIslandHMM(char **aseqs, AINFO *ainfo, int idx, 
				     float null[MAXABET], float ***mx,
				     float **bprior, float *qprior, int nmx);

/* 
 * from core_algorithms.c
 * Clean research/demonstration versions of basic algorithms.
 */
extern struct dpmatrix_s *AllocPlan7Matrix(int rows, int M, int ***xmx, 
					   int ***mmx, int ***imx, int ***dmx);
extern void  FreePlan7Matrix(struct dpmatrix_s *mx);
extern float Plan7Forward(char *dsq, int L, struct plan7_s *hmm, 
			  struct dpmatrix_s **ret_mx);
extern float Plan7Viterbi(char *dsq, int L, struct plan7_s *hmm, 
			  struct dpmatrix_s **ret_mx);
extern void  P7ViterbiTrace(struct plan7_s *hmm, char *dsq, int L,
			   struct dpmatrix_s *mx, struct p7trace_s **ret_tr);
extern float Plan7ESTViterbi(char *dsq, int L, struct plan7_s *hmm, 
			     struct dpmatrix_s **ret_mx);

/* 
 * from debug.c
 * Debugging output of various sorts.
 */
extern void VerboseWorry(int level, char *file, int line, char *fmt, ...);
extern char *Statetype(enum p7stype st);
extern void P7PrintTrace(FILE *fp, struct p7trace_s *tr, 
			 struct plan7_s *hmm, char *dsq);
extern void P7PrintPrior(FILE *fp, struct p7prior_s *pri);

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
extern void ExtremeValueSetHistogram(struct histogram_s *h, float mu, 
				     float lambda, int ndegrees);
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
extern void     WriteAscHMM(FILE *fp, struct plan7_s *hmm);
extern void     WriteBinHMM(FILE *fp, struct plan7_s *hmm);

/* masks.c
 * Repetitive sequence masking.
 */
extern void XNU(char *dsq, int len);

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

/* 
 * from misc.c
 * Miscellaneous functions with no home
 */
extern void  Banner(FILE *fp, char *banner);
extern void  BlockRaggedEdgedAlignment(char **aseqs, int nseq, int alen);
extern int   AlignmentTooBig(int L, int M);
extern void  SuppressChatter(struct hmm_struc *hmm);
extern char *Getword(FILE *fp, int type); 
extern char *Getline(char *s, int n, FILE *fp);


/* 
 * from modelmakers.c
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
 * Experimental: Plan7 HMM structure support
 */
extern struct plan7_s *AllocPlan7(int M);
extern void FreePlan7(struct plan7_s *hmm);
extern void ZeroPlan7(struct plan7_s *hmm);
extern void Plan7SetName(struct plan7_s *hmm, char *name);
extern void Plan7SetDescription(struct plan7_s *hmm, char *desc);
extern void Plan7SetComline(struct plan7_s *hmm, int argc, char **argv);
extern void Plan7SetCtime(struct plan7_s *hmm);
extern void Plan7SetNullModel(struct plan7_s *hmm, float null[MAXABET], float p1);
extern void Plan7Logoddsify(struct plan7_s *hmm);
extern void Plan7Probify(struct plan7_s *hmm);
extern void Plan7Renormalize(struct plan7_s *hmm);
extern void Plan7GlobalConfig(struct plan7_s *hmm);
extern void Plan7LSConfig(struct plan7_s *hmm);
extern void Plan7SWConfig(struct plan7_s *hmm, float pentry, float pexit);
extern void PrintPlan7Stats(FILE *fp, struct plan7_s *hmm, char **dsq, 
			    int nseq, struct p7trace_s **tr);
extern int  DegenerateSymbolScore(float *p, float *null, int ambig);
extern void Plan9toPlan7(struct hmm_struc *hmm, struct plan7_s **ret_plan7);
extern void Plan7toPlan9Search(struct plan7_s *hmm, struct shmm_s **ret_shmm);

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

/* 
 * from states.c
 * Support for the basic data structures
 */
extern struct hmm_struc *AllocHMM(int M);
extern void ZeroHMM(struct hmm_struc *hmm);
extern void LogifyHMM(struct hmm_struc *hmm);
extern void LogoddsifyHMM(struct hmm_struc *hmm);
extern int  WriteFlatPriorHMM(struct hmm_struc *hmm, struct prior_s *prior);
extern struct hmm_struc *HMMDup(struct hmm_struc *hmm);
extern void HMMCopy(struct hmm_struc *hmm1, struct hmm_struc *hmm2);
extern int  FreeHMM(struct hmm_struc *hmm);
extern struct shmm_s *AllocSearchHMM(int M);
extern void  FreeSearchHMM(struct shmm_s *shmm);
extern int  CountSymbol(char sym, float wt, float *counters);
extern float HMMDistance(struct hmm_struc *newhmm, struct hmm_struc *oldhmm);
extern void Renormalize(struct hmm_struc *hmm);

/* from tophits.c
 * Support for keeping/sorting top scoring hit/alignment lists
 */
extern struct tophit_s *AllocTophits(int H, int A);
extern void   FreeTophits(struct tophit_s *hitlist);
extern struct fancyali_s *AllocFancyAli(void);
extern void   FreeFancyAli(struct fancyali_s *ali);
extern void   RegisterHit(struct tophit_s *hitlist, 
			  double sortkey, double evalue, float score, 
			  char *name, char *desc, int sqfrom, int sqto, int sqlen, 
			  int hmmfrom, int hmmto, int hmmlen, 
			  int domidx, int ndom, 
			  struct fancyali_s *ali);
extern void GetRankedHit(struct tophit_s *h, int rank, 
			 double *r_evalue, float *r_score, 
			 char **r_name, char **r_desc,
			 int *r_sqfrom, int *r_sqto, int *r_sqlen,
			 int *r_hmmfrom, int *r_hmmto, int *r_hmmlen,
			 int *r_domidx, int *r_ndom,
			 struct fancyali_s **r_ali);
extern int    TophitsMaxName(struct tophit_s *h);
extern void   FastSortTophits(struct tophit_s *h);
extern void   FullSortTophits(struct tophit_s *h);

/* from trace.c
 * Support for traceback (state path) structure
 */
extern void  P7AllocTrace(int tlen, struct p7trace_s **ret_tr);
extern void  P7ReallocTrace(struct p7trace_s *tr, int tlen);
extern void  P7FreeTrace(struct p7trace_s *tr);
extern void  P7ReverseTrace(struct p7trace_s *tr);
extern void  P7TraceCount(struct plan7_s *hmm, char *dsq, float wt, 
			  struct p7trace_s *tr);
extern float P7TraceScore(struct plan7_s *hmm, char *dsq, struct p7trace_s *tr);
extern void  P7Traces2Alignment(char **dsq, SQINFO *sqinfo, float *wgt, int nseq, int M, 
				struct p7trace_s **tr, int matchonly, 
				char ***ret_aseqs, AINFO *ainfo);
extern int  TransitionScoreLookup(struct plan7_s *hmm, enum p7stype st1, 
				  int k1, enum p7stype st2, int k2);
extern struct fancyali_s *CreateFancyAli(struct p7trace_s *tr, struct plan7_s *hmm,
					 char *dsq, char *name);
extern void PrintFancyAli(FILE *fp, struct fancyali_s *ali);
extern void TraceDecompose(struct p7trace_s *otr, struct p7trace_s ***ret_tr,
			   int *ret_ntr);
extern int  TraceDomainNumber(struct p7trace_s *tr);
extern void TraceSimpleBounds(struct p7trace_s *tr, int *ret_i1, int *ret_i2, 
			      int *ret_k1,  int *ret_k2);


#endif /*FUNCSH_INCLUDED*/
