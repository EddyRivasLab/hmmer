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

/* 
 * from align.c
 * Make multiple alignments from individual tracebacks
 */
extern int Traces2Alignment(char **rseqs, SQINFO *sqinfo, int nseq, int M, 
			    struct trace_s **tr, int matchonly, 
			    char ***ret_aseqs, AINFO *ret_ainfo);

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
 * from p7debug.c
 * Debugging output of various sorts.
 */
extern char *Statetype(enum p7stype st);
extern void P7PrintTrace(FILE *fp, struct p7trace_s *tr, 
			 struct plan7_s *hmm, char *dsq);
extern void P7PrintPrior(FILE *fp, struct p7prior_s *pri);

/* 
 * from emit.c
 * Emit sequences from a model.
 */
extern int EmitSequence(struct hmm_struc *hmm, char **ret_seq, 
			struct trace_s **ret_tr);
extern int EmitBestSequence(struct hmm_struc *hmm, char **ret_seq, 
			    struct trace_s **ret_tr);
extern int RD_EmitBestSequence(struct hmm_struc *hmm, char **ret_seq, 
			       struct trace_s **ret_tr);
extern float  HMMExpectedScore(struct hmm_struc *hmm, float **ret_pernode);
extern float  HMMExpectRandomScore(struct hmm_struc *hmm);
extern float  HMMSimpleExpectRandom(struct hmm_struc *hmm);


/* from emulation.c
 * Emulation of GCG PROFILE package by Michael Gribskov
 * Emulation of classical pairwise non-probabilistic alignment
 */
extern int WriteProfile(FILE *fp, struct hmm_struc *hmm, struct shmm_s *shmm); 
extern struct shmm_s *PairwiseEmulation(char *seqfile, char *pamfile, 
					float gop, float gex);
extern struct hmm_struc *PairHMMEmulation(char *seqfile, char *pamfile, 
					  float gop, float gex);
extern float **PAM2SubstitutionMatrix(int **pam, float scale, float *pri);
extern float **PAM2Joint(int **pam, float scale, float *pri);

/* from forback.c
 * Forward-backward algorithm
 */
extern void  Forward(struct hmm_struc *hmm, char *s, struct forback_s ***ret_fwd, float **ret_scl);
extern void  Backward(struct hmm_struc *hmm, char *s, float *scl, struct forback_s ***ret_bck);
extern void  AlignmentConfidence(struct hmm_struc *hmm, int L, 
				 struct forback_s **fwd, struct forback_s **bck, float *scl,
				 struct forback_s ***ret_conf);
extern void  TraceConfidence(struct forback_s **cmx, int L, struct trace_s *tr, float **ret_conf);
extern void  ForbackCount(struct hmm_struc *hmm, char *seq, int L, float weight,
			  struct forback_s **fwd, struct forback_s **bck,
			  float *scl, struct hmm_struc *count);
extern float ForwardScore(struct forback_s **fwd, int L, int M, float *scl);
extern float BackwardScore(struct forback_s **bck, int L, float *scl);
extern float RandomScore(float *randomseq, char *seq);
extern float ForbackSymscore(char x, float *scores, int hyperbayes);
extern void  DumpForbackMatrix(struct forback_s **mx, int L, int M, float *scalefactors);


/* from fragviterbi.c
 * alignment algorithm: multiple-hit Smith/Waterman scanning version (hmmfs)
 */
extern float FragViterbi(struct shmm_s *shmm, char *seq, int L, int singlehit, int hmmls,
			 float P1, float P2, float P3, float thresh,
			 int (*gotone_f)(struct shmm_s *, char *, int,int,int,float));

/* from histogram.c
 * accumulation of scores
 */
extern struct histogram_s *AllocHistogram(int min, int max, int lumpsize);
extern void FreeHistogram(struct histogram_s *h);
extern void UnfitHistogram(struct histogram_s *h);
extern void AddToHistogram(struct histogram_s *h, float sc);
extern void PrintASCIIHistogram(FILE *fp, struct histogram_s *h); 
extern int  ExtremeValueFitHistogram(struct histogram_s *h, float high_hint);
extern void ExtremeValueSetHistogram(struct histogram_s *h, float mu, float lambda);
extern int  GaussianFitHistogram(struct histogram_s *h, float high_hint);
extern void GaussianSetHistogram(struct histogram_s *h, float mean, float sd);
extern double ExtremeValueP (float x, float mu, float lambda);
extern double ExtremeValueP2(float x, float mu, float lambda, int N);
extern double ExtremeValueE (float x, float mu, float lambda, int N);
extern float  EVDrandom(float mu, float lambda);

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
 * from maxmodelmaker.c
 * Construction of models from multiple alignments
 */
extern void Maxmodelmaker(char **aseqs, AINFO *ainfo, int nseq,
			  struct prior_s *prior, float *randomseq, float mpri, 
			  struct hmm_struc **ret_hmm, struct trace_s ***tr);
extern void Handmodelmaker(char **aseqs, AINFO *ainfo, int nseq,
			   struct prior_s *prior, struct hmm_struc **ret_hmm, 
			   struct trace_s ***tr);

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

/* from output.c
 * "User-friendly" output
 */
extern void PrintFancyTrace(FILE *ofp, struct shmm_s *shmm, 
			    struct trace_s *tr,
			    char *seq, char *seqname, int from_pos);

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

/* from p7prior.c
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
/* from p7trace.c
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

/* OBSOLETE FUNCS: */
extern void P7PrintFancyTrace(FILE *fp, struct p7trace_s *tr, struct plan7_s *hmm,
			      char *dsq, char *name);


/* 
 * from prior.c
 * Dirichlet priors
 */
extern void DefaultSimplePrior(struct prior_s **ret_prior);
extern void ReadPrior(char *pfile, struct prior_s **ret_prior);
extern void ToPAMPrior(int **pam, float scale, float wt, struct prior_s *prior);
extern void PrintPAMPrior(struct prior_s *pri);
extern struct prior_s *AllocPrior(void);
extern void FreePrior(struct prior_s *prior);
extern void DefaultNullModel(float *null);
extern void ReadNullModel(char *fname, float *null);
extern void PriorifyMatchVector(float *vec, struct prior_s *prior,
				float *ret_mix);
extern void PriorifyInsertVector(float *vec, struct prior_s *prior);
extern void PriorifyTransitionVectors(float *tm, float *ti, float *td, 
				      struct prior_s *prior);
extern void PriorifyHMM(struct hmm_struc *hmm, struct prior_s *prior);
extern void SpecialPriorifyHMM(struct hmm_struc *hmm, struct prior_s *prior);
extern void StructurePerceptron(struct prior_s *prior, float *xray);
extern void AnnotateAlignment(char **aseq, int nseq, AINFO *ainfo, 
			      float **ret_inputs);
extern float P_ModelGivenDirichlet(struct hmm_struc *hmm,struct prior_s *pri);



/* 
 * from saviterbi.c
 * Alignment algorithm: simulated annealing version
 */
extern int    SaFill(struct sa_hmm_s *sahmm, char *s, struct sa_s ***ret_mx);
extern int    SaTrace(struct sa_s **mx, int L, struct sa_hmm_s *sahmm, 
		      struct trace_s **ret_tr);
extern void   DumpSaMatrix(struct sa_s **mx, int L, int M, 
			   double *scalefactors);
extern struct sa_hmm_s *CreateSahmm(struct hmm_struc *hmm, float kT);
extern void   DestroySahmm(struct sa_hmm_s *sahmm);
extern double SaSymscore(char x, double *scores, int hyperbayes);


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

/* 
 * from swviterbi.c
 * Alignment algorithm: single-hit Smith/Waterman version
 */
extern void SWViterbi(struct shmm_s *shmm, char *seq, int L, 
		      float P1, float P2, float P3, int fullseq,
		      int *ret_i, int *ret_j, int *ret_kstart, int *ret_kend, 
		      float *ret_score, struct trace_s **ret_tr);


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
			  int hmmfrom, int hmmto, int hmmlen, struct fancyali_s *ali);
extern void GetRankedHit(struct tophit_s *h, int rank, 
			 double *r_evalue, float *r_score, 
			 char **r_name, char **r_desc,
			 int *r_sqfrom, int *r_sqto, int *r_sqlen,
			 int *r_hmmfrom, int *r_hmmto, int *r_hmmlen,
			 struct fancyali_s **r_ali);
extern int    TophitsMaxName(struct tophit_s *h, int top_howmany);
extern void   PrintTopHits(FILE *fp, struct tophit_s *h, int n, int evd, 
			   float mu, float lambda, int dbseqs);
extern void   FastSortTophits(struct tophit_s *h);
extern void   FullSortTophits(struct tophit_s *h);


/* from trace.c
 * Support for the alignment traceback data structure
 */
extern void AllocTrace(int tlen, struct trace_s **ret_tr);
extern void ReallocTrace(struct trace_s *tr, int tlen);
extern void FreeTrace(struct trace_s *tr);
extern void ReverseTrace(struct trace_s *tr, int tlen);
extern void PrintTrace(struct trace_s *tr);
extern void TraceCount(struct hmm_struc *hmm, char *seq, float wt, struct trace_s *tr);
extern int  TraceScore(struct shmm_s *shmm, char *seq, struct trace_s *tr,
		       float *ret_score);
extern void DealignTrace(struct trace_s *tr, char *aseq, int alen);
extern void PositionAverageScore(struct shmm_s *shmm, char **seq, float *wgt,
				 int nseq, struct trace_s **tr, 
				 float **ret_scarr, float *ret_avg);

/* from viterbi.c
 * Alignment algorithm: global alignment (hmms, hmma, hmmt)
 */
extern int   Symscore(char x, float *scores, float *priors);
extern void  MakeSearchHMM(struct hmm_struc *hmm, struct shmm_s *shmm);
extern struct shmm_s *MakePAMSearchHMM(char *seq, SQINFO *sqinfo, int **pam, 
				       float scale, float gop, float gex);
extern void  PrintSearchHMM(FILE *fp, struct shmm_s *shmm);
extern void  PrepareSequence(char *s, char **ret_seq, int *ret_L);
extern int   ViterbiFill(struct shmm_s *shmm, char *seq, int L, 
			 struct vit_s ***ret_mx, float *ret_score);
extern int   ViterbiTrace(struct vit_s **mx, struct shmm_s *shmm, char *seq, int window,
			  int end_i, int end_k, 
			  struct trace_s **ret_tr, int *ret_i, int *ret_k);
extern void  FreeViterbiMatrix(struct vit_s **mx, int L);	
extern void  PrintViterbiMatrix(struct vit_s **mx, char *seq1, int L, int M);
extern void  ViterbiAlignAlignment(struct shmm_s *shmm, char **aseqs, int alen, int nseq,
				   struct trace_s ***ret_tr, float *ret_sc);
extern void PrintFragViterbiMatrix(struct fvit_s **mx, int L, int M);



/* from weeviterbi.c
 * Linear memory alignment for Needleman/Wunsch and full seq Smith Waterman
 */
extern void WeeViterbi(struct shmm_s *shmm, char *seq, int L, 
		       int smithwaterman, float P2, float P3,
		       struct trace_s **ret_tr, float *ret_sc);



#endif /*FUNCSH_INCLUDED*/
