/* P7_PIPELINE is the standardized pipeline for one profile/sequence
 * comparison, from the fast filters down through domain postprocessing,
 * alignment, and scoring.
 */
#ifndef p7PIPELINE_INCLUDED
#define p7PIPELINE_INCLUDED

#include <p7_config.h>

#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"

#include "base/p7_bg.h"
#include "base/p7_hmmfile.h"
#include "base/p7_hmmwindow.h"
#include "base/p7_masstrace.h"
#include "base/p7_scoredata.h"
#include "base/p7_tophits.h"

#include "dp_sparse/p7_sparsemx.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/p7_filtermx.h"


enum p7_pipemodes_e { p7_SEARCH_SEQS = 0, p7_SCAN_MODELS = 1 };
enum p7_zsetby_e    { p7_ZSETBY_NTARGETS = 0, p7_ZSETBY_OPTION = 1, p7_ZSETBY_FILEINFO = 2 };
enum p7_complementarity_e { p7_NOCOMPLEMENT    = 0, p7_COMPLEMENT   = 1 };

/* We break out a separate chunk of structure, P7_PIPELINE_STATS, to
 * contain accounting count accumulators; in threaded/mpi
 * implementations, workers need to send these data back to be folded
 * into the master, but the rest of each worker's pipeline does not
 * need to be communicated. See p7_pipeline_stats_Merge(), and MPI
 * comm routines for P7_PIPELINE_STATS in p7_pipeline_mpi.c.
 */
typedef struct p7_pipeline_stats_s {
  uint64_t      nmodels;        /* # of HMMs searched                       */
  uint64_t      nseqs;	        /* # of sequences searched                  */
  uint64_t      nres;	        /* # of residues searched                   */
  uint64_t      nnodes;	        /* # of model nodes searched                */
  uint64_t      n_past_ssv;	/* # comparisons that pass MSVFilter()      */
  uint64_t      n_past_bias;	/* # comparisons that pass bias filter      */
  uint64_t      n_past_vit;	/* # comparisons that pass ViterbiFilter()  */
  uint64_t      n_past_fwd;	/* # comparisons that pass ForwardFilter()  */

  /* Additional accounting in nhmmer */
  uint64_t      n_output;	/* # alis that make it to final output      */
  uint64_t      pos_past_msv;	/* # positions that pass MSVFilter()        */
  uint64_t      pos_past_bias;	/* # positions that pass bias filter        */
  uint64_t      pos_past_vit;	/* # positions that pass ViterbiFilter()    */
  uint64_t      pos_past_fwd;	/* # positions that pass ForwardFilter()    */
  uint64_t      pos_output;	/* # positions that make it to final output */
} P7_PIPELINE_STATS;



/* P7_PIPELINE
 */
typedef struct p7_pipeline_s {
  /* Dynamic programming matrices                                           */
  P7_FILTERMX   *fx;	        /* one-row vector DP: MSV, Vit filter       */
  P7_CHECKPTMX  *cx;		/* checkpointed vector Fwd/Bck decoding     */
  P7_SPARSEMASK *sm;		/* sparse mask created by F/B decoding      */
  P7_SPARSEMX   *sxf;		/* sparse Forward (glocal/local)            */
  P7_SPARSEMX   *sxb;		/* sparse Backward                          */
  P7_SPARSEMX   *sxd;		/* sparse Decoding                          */
  P7_SPARSEMX   *sxx;		/* sparse Viterbi, plus other uses          */
  P7_TRACE      *tr;		/* Viterbi trace of seq/profile comparison  */
  P7_MASSTRACE  *mt;		/* envelope endpoint determination          */
  float         *n2sc;		/* null2 scores per position, 1..L          */
  float         *wrk;		/* null2 calculation needs M+1 wrkspace     */

  /* Reporting threshold settings                                           */
  int     by_E;		        /* TRUE to cut per-target report off by E   */
  double  E;	                /* per-target E-value threshold             */
  double  T;	                /* per-target bit score threshold           */
  int     dom_by_E;             /* TRUE to cut domain reporting off by E    */
  double  domE;	                /* domain E-value threshold                 */
  double  domT;	                /* domain bit score threshold               */
  int     use_bit_cutoffs;      /* (FALSE | p7H_GA | p7H_TC | p7H_NC)       */

  /* Inclusion threshold settings                                           */
  int     inc_by_E;		/* TRUE to threshold inclusion by E-values  */
  double  incE;			/* per-target inclusion E-value threshold   */
  double  incT;			/* per-target inclusion score threshold     */
  int     incdom_by_E;		/* TRUE to threshold domain inclusion by E  */
  double  incdomE;		/* per-domain inclusion E-value threshold   */
  double  incdomT;		/* per-domain inclusion E-value threshold   */

  /* Tracking search space sizes for E value calculations                   */
  double  Z;			/* eff # targs searched (per-target E-val)  */
  double  domZ;			/* eff # signific targs (per-domain E-val)  */
  enum p7_zsetby_e Z_setby;   	/* how Z was set                            */
  enum p7_zsetby_e domZ_setby;	/* how domZ was set                         */
  
  /* Threshold settings for pipeline                                        */
  int     do_max;	        /* TRUE to run in slow/max mode             */
  double  F1;		        /* MSV filter threshold                     */
  double  F2;		        /* Viterbi filter threshold                 */
  double  F3;		        /* uncorrected Forward filter threshold     */
  int     B1;                   /* window len, biased-comp modifier - MSV   */
  int     B2;                   /* window len, biased-comp modifier - Vit   */
  int     B3;                   /* window len, biased-comp modifier - Fwd   */
  int     do_biasfilter;	/* TRUE to use biased comp HMM filter       */
  int     do_null2;		/* TRUE to use null2 score corrections      */

  /* State config, flags */
  enum p7_pipemodes_e mode;    	/* p7_SCAN_MODELS | p7_SEARCH_SEQS          */
  int           show_accessions;/* TRUE to output accessions not names      */
  int           show_alignments;/* TRUE to output alignments (default)      */
  P7_HMMFILE   *hfp;		/* COPY of open HMM db (if hmmscan mode)    */

  /* Additional state config for nhmmer */
  int           long_targets;   /* TRUE if targ seqs expected to be v. long                     */
  enum p7_strands_e strand;     /* p7_STRAND_TOPONLY | p7_STRAND_BOTTOMONLY | p7_STRAND_BOTH    */
  int 		W;              /* window len for nhmmer; ~max len expected                     */
  int           block_length;   /* overlapping block len, threaded; p7_NHMMER_MAX_RESIDUE_COUNT */

  /* Accounting. (reduceable in threaded/MPI parallel version)              */
  P7_PIPELINE_STATS stats;	/* # of things searched, things past filters*/

  /* Diagostic info */
  char          errbuf[eslERRBUFSIZE];
} P7_PIPELINE;



extern P7_PIPELINE *p7_pipeline_Create(ESL_GETOPTS *go, int M_hint, int L_hint, int do_longtargets, enum p7_pipemodes_e mode);
extern int          p7_pipeline_Reuse  (P7_PIPELINE *pli);
extern void         p7_pipeline_Destroy(P7_PIPELINE *pli);

extern int p7_pipeline_NewModel          (P7_PIPELINE *pli, const P7_OPROFILE *om, P7_BG *bg);
extern int p7_pipeline_NewModelThresholds(P7_PIPELINE *pli, const P7_OPROFILE *om);
extern int p7_pipeline_NewSeq            (P7_PIPELINE *pli, const ESL_SQ *sq);

extern int p7_pipeline_TargetReportable  (P7_PIPELINE *pli, float score,     double lnP);
extern int p7_pipeline_DomainReportable  (P7_PIPELINE *pli, float dom_score, double lnP);
extern int p7_pipeline_TargetIncludable  (P7_PIPELINE *pli, float score,     double lnP);
extern int p7_pipeline_DomainIncludable  (P7_PIPELINE *pli, float dom_score, double lnP);

extern int p7_pipeline_stats_Init (P7_PIPELINE_STATS *stats);
extern int p7_pipeline_stats_Merge(P7_PIPELINE *p1, const P7_PIPELINE_STATS *stats);
extern int p7_pipeline_WriteStats(FILE *ofp, P7_PIPELINE *pli, ESL_STOPWATCH *w);

extern int p7_Pipeline              (P7_PIPELINE *pli, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, P7_TOPHITS *th);
extern int p7_Pipeline_LongTarget   (P7_PIPELINE *pli, P7_PROFILE *gm, P7_OPROFILE *om, P7_SCOREDATA *msvdata, P7_BG *bg, const ESL_SQ *sq, P7_TOPHITS *hitlist, int64_t seqidx);

extern int p7_pipeline_AccelerationFilter(ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_BG *bg,
					  P7_FILTERMX *fx, P7_CHECKPTMX *cx, P7_SPARSEMASK *sm);

#endif /*p7PIPELINE_INCLUDED*/

