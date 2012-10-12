enum p7_pipemodes_e { p7_SEARCH_SEQS = 0, p7_SCAN_MODELS = 1 };
enum p7_zsetby_e    { p7_ZSETBY_NTARGETS = 0, p7_ZSETBY_OPTION = 1, p7_ZSETBY_FILEINFO = 2 };

typedef struct p7_pipeline_s {
  /* Dynamic programming matrices                                           */
  P7_OMX     *oxf;		/* one-row Forward matrix, accel pipe       */
  P7_OMX     *oxb;		/* one-row Backward matrix, accel pipe      */
  P7_OMX     *fwd;		/* full Fwd matrix for domain envelopes     */
  P7_OMX     *bck;		/* full Bck matrix for domain envelopes     */

  /* Domain postprocessing                                                  */
  ESL_RANDOMNESS *r;		/* random number generator                  */
  int             do_reseeding; /* TRUE: reseed for reproducible results    */
  P7_DOMAINDEF   *ddef;		/* domain definition workflow               */

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
  int     B1;               /* window length for biased-composition modifier - MSV*/
  int     B2;               /* window length for biased-composition modifier - Viterbi*/
  int     B3;               /* window length for biased-composition modifier - Forward*/
  int     do_biasfilter;	/* TRUE to use biased comp HMM filter       */
  int     do_null2;		/* TRUE to use null2 score corrections      */

  /* Accounting. (reduceable in threaded/MPI parallel version)              */
  uint64_t      nmodels;        /* # of HMMs searched                       */
  uint64_t      nseqs;	        /* # of sequences searched                  */
  uint64_t      nres;	        /* # of residues searched                   */
  uint64_t      nnodes;	        /* # of model nodes searched                */
  uint64_t      n_past_msv;	/* # comparisons that pass MSVFilter()      */
  uint64_t      n_past_bias;	/* # comparisons that pass bias filter      */
  uint64_t      n_past_vit;	/* # comparisons that pass ViterbiFilter()  */
  uint64_t      n_past_fwd;	/* # comparisons that pass ForwardFilter()  */
  uint64_t      n_output;	    /* # alignments that make it to the final output (used for nhmmer) */
  uint64_t      pos_past_msv;	/* # positions that pass MSVFilter()  (used for nhmmer) */
  uint64_t      pos_past_bias;	/* # positions that pass bias filter  (used for nhmmer) */
  uint64_t      pos_past_vit;	/* # positions that pass ViterbiFilter()  (used for nhmmer) */
  uint64_t      pos_past_fwd;	/* # positions that pass ForwardFilter()  (used for nhmmer) */
  uint64_t      pos_output;	    /* # positions that make it to the final output (used for nhmmer) */

  enum p7_pipemodes_e mode;    	/* p7_SCAN_MODELS | p7_SEARCH_SEQS          */
  int           long_targets;   /* TRUE if the target sequences are expected to be very long (e.g. dna chromosome search in nhmmer) */
  int           strand;         /* TRUE if the search should ignore the revcomp (used for nhmmer only) */
  int 		    	W;              /* window length for nhmmer scan - essentially maximum length of model that we expect to find*/
  int           block_length;   /* length of overlapping blocks read in the multi-threaded variant (default MAX_RESIDUE_COUNT) */

  int           show_accessions;/* TRUE to output accessions not names      */
  int           show_alignments;/* TRUE to output alignments (default)      */

  P7_HMMFILE   *hfp;		/* COPY of open HMM database (if scan mode) */
  char          errbuf[eslERRBUFSIZE];
} P7_PIPELINE;



extern P7_PIPELINE *p7_pipeline_Create(ESL_GETOPTS *go, int M_hint, int L_hint, int long_targets, enum p7_pipemodes_e mode);
extern int          p7_pipeline_Reuse  (P7_PIPELINE *pli);
extern void         p7_pipeline_Destroy(P7_PIPELINE *pli);
extern int          p7_pipeline_Merge  (P7_PIPELINE *p1, P7_PIPELINE *p2);

extern int p7_pli_ExtendAndMergeWindows (P7_OPROFILE *om, const P7_SCOREDATA *msvdata, P7_HMM_WINDOWLIST *windowlist, int L, float pct_overlap);
extern int p7_pli_TargetReportable  (P7_PIPELINE *pli, float score,     double lnP);
extern int p7_pli_DomainReportable  (P7_PIPELINE *pli, float dom_score, double lnP);

extern int p7_pli_TargetIncludable  (P7_PIPELINE *pli, float score,     double lnP);
extern int p7_pli_DomainIncludable  (P7_PIPELINE *pli, float dom_score, double lnP);
extern int p7_pli_NewModel          (P7_PIPELINE *pli, const P7_OPROFILE *om, P7_BG *bg);
extern int p7_pli_NewModelThresholds(P7_PIPELINE *pli, const P7_OPROFILE *om);
extern int p7_pli_NewSeq            (P7_PIPELINE *pli, const ESL_SQ *sq);
extern int p7_Pipeline              (P7_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, P7_TOPHITS *th);
extern int p7_Pipeline_LongTarget   (P7_PIPELINE *pli, P7_OPROFILE *om, P7_SCOREDATA *msvdata, P7_BG *bg, const ESL_SQ *sq, P7_TOPHITS *hitlist, int64_t seqidx);
extern int p7_Pipeline_FM           (P7_PIPELINE *pli, P7_OPROFILE *om, P7_SCOREDATA *msvdata, P7_BG *bg, P7_TOPHITS *hitlist, int64_t seqidx,
                                     const FM_DATA *fmf, const FM_DATA *fmb, FM_CFG *fm_cfg);

extern int p7_pli_Statistics(FILE *ofp, P7_PIPELINE *pli, ESL_STOPWATCH *w);

