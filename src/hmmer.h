
/* Master include file; this includes the entire HMMER3 header tree.
 * Programs can be built against a HMMER3 library (libhmmer.a) by including this
 * single header file. 
 *
 * Also serves as an overview of the package's source code.
 */
#include "p7_config.h"

/* 'base' subdir: contains code for many of HMMER3's data structures                                                             */
#include "base/general.h"	     /* widely used definitions, error handling, initialization                                  */
#include "base/p7_hmm.h"	     /* P7_HMM        : basic model, probability parameters or counts                            */
#include "base/p7_hmm_mpi.h"	     /*               :    ... add-on: MPI communication                                         */
#include "base/p7_profile.h"	     /* P7_PROFILE    : search model, glocal/local, with additional states for nonhomologous seq */
#include "base/p7_profile_mpi.h"     /*               :    ... add-on: MPI communication                                         */
#include "base/p7_hmmfile.h"	     /* P7_HMMFILE    : reading models from files                                                */
#include "base/p7_trace.h"	     /* P7_TRACE      : alignment of a model to a sequence: an HMM state path                    */
#include "base/p7_bg.h"		     /* P7_BG         : null model, of an entirely nonhomologous target seq                      */
#include "base/p7_prior.h"	     /* P7_PRIOR      : Dirichlet mixture prior on model parameters                              */
#include "base/p7_masstrace.h"	     /* P7_MASSTRACE  : used in calculating envelope bounds for a domain                         */
#include "base/p7_domain.h"	     /* P7_DOMAIN     : information about a match to a model in a target seq                     */
#include "base/p7_alidisplay.h"	     /* P7_ALIDISPLAY : an alignment formatted for output                                        */
#include "base/p7_tophits.h"	     /* P7_HIT, P7_TOPHITS : accumulated information about hits (scores, alis) during search     */
#include "base/p7_hmmwindow.h"	     /* P7_HMM_WINDOW, P7_HMM_WINDOWLIST : {nhmmer}                                              */
#include "base/p7_scoredata.h"	     /* P7_SCOREDATA  : {nhmmer}                                                                 */

/* 'misc' subdir: various other support functions                                                                         */
#include "misc/emit.h"		           /* emitting (sampling) sequences from HMM or profile                           */
#include "misc/mpisupport.h"	           /* MPI (Message Passing Interface) support                                     */
#include "misc/h2_io.h"		           /* Legacy support for HMMER2 model file formats                                */
#include "misc/p7_trace_metrics.h"         /* Benchmarking utilities for alignment accuracy                               */
#include "misc/logsum.h"	           /* Fast lookup-table-driven log-sum-exp-2 function                             */
#include "misc/tracealign.h"	           /* Conversion of trace structures to multiple alignment                        */

/* 'build' subdir: building and calibrating a new model */
#include "build/build.h"	           /* Constructing a model from a multiple seq alignment                          */
#include "build/seqmodel.h"	           /* Constructing a model from a single query seq                                */
#include "build/evalues.h"	           /* Model calibration, setting parameters for E-value determination             */
#include "build/eweight.h"	           /* Entropy-weighting: ad hoc absolute sequence weights                         */
#include "build/modelsample.h"	           /* Sampling models randomly; used extensively in unit testing                  */
#include "build/modelstats.h"	           /* Various summary statistics for a model                                      */
#include "build/p7_builder.h"	           /* P7_BUILDER: an aggregrated pipeline for building new models                 */

/* 'search' subdir: compare a model against a sequence, to calculate score and alignment(s)                               */
#include "search/modelconfig.h"	           /* Building a search profile (P7_PROFILE) from a basic model (P7_HMM)          */
#include "search/tophits_output.h"	   /* Formatted output of search results, "human-readable" verbose                */
#include "search/tophits_output_tabular.h" /* Formatted output of search results, tabular                                 */
#include "search/null3.h"		   /* {nhmmer}: 'null3' model for correcting biased composition                   */
#include "search/p7_pipeline.h"	           /* P7_PIPELINE: an aggregated pipeline for comparing one model against one seq */

/* 'dp_vector' subdir: SIMD-vector-parallel accelerated dynamic programming filters; H3's speed heuristics                */
#include "dp_vector/simdvec.h"	           /* general definitions and utilities for SIMD vector code                      */
#include "dp_vector/p7_oprofile.h"	   /* P7_OPROFILE  : striped/vectorized profile                                   */
#include "dp_vector/p7_filtermx.h"	   /* P7_FILTERMX  : a one-row O(M) memory DP matrix for SSV, MSV, VF             */
#include "dp_vector/p7_checkptmx.h"	   /* P7_CHECKPTMX : checkpointed O(M sqrt(L)) memory DP for F/B/Decode           */
#include "dp_vector/ssvfilter.h"           /* "single segment ungapped Viterbi" filter, from Bjarne Knudsen               */
#include "dp_vector/msvfilter.h"	   /* "multiple seqment ungapped Viterbi" filter                                  */
#include "dp_vector/vitfilter.h"	   /* Viterbi filter                                                              */
#include "dp_vector/fwdfilter.h"	   /* Checkpointed Forward/Backward/posterior decoding filter                     */
#include "dp_vector/io.h"		   /* Reading/writing P7_OPROFILEs to files                                       */
#include "dp_vector/p7_oprofile_mpi.h"     /* MPI support functions for communicating P7_OPROFILE                         */

/* 'dp_sparse' subdir: glocal/local alignment in sparse dynamic programming */
#include "dp_sparse/p7_sparsemx.h"         /* P7_SPARSEMASK, P7_SPARSEMX : sparse mask and sparse DP matrix  */
#include "dp_sparse/sparse_fwdback.h"	   /* Viterbi glocal/local alignment in sparse DP                    */
#include "dp_sparse/sparse_viterbi.h"	   /* Forward, Backward glocal/local in sparse DP                    */
#include "dp_sparse/sparse_decoding.h"	   /* Posterior decoding in sparse DP                                */
#include "dp_sparse/sparse_trace.h"	   /* glocal/local traceback (Viterbi, stochastic)                   */
#include "dp_sparse/sparse_masstrace.h"	   /* 'Mass trace' algorithm for determining envelope bounds         */
#include "dp_sparse/sparse_envscore.h"	   /* Modified Forward algorithm for scoring a single envelope       */
#include "dp_sparse/sparse_null2.h"	   /* 'null2' algorithm for compensating for biased composition      */

/* 'dp_reference' subdir: "reference" implementations of key algorithms, used for regression and other testing */
#include "dp_reference/p7_refmx.h"            /* P7_REFMX: a O(ML) DP matrix for reference DP algorithms       */
#include "dp_reference/reference_fwdback.h"   /* reference implementation of glocal/local Forward, Backward    */
#include "dp_reference/reference_viterbi.h"   /* reference implementation of glocal/local Viterbi              */
#include "dp_reference/reference_decoding.h"  /* reference implementation of glocal/local posterior Decoding   */
#include "dp_reference/reference_trace.h"     /* reference implementation of glocal/local tracebacks           */

/* 'daemon' subdir: Michael Farrar's daemon for HMMER3 web services                                       */
#include "daemon/p7_hmmcache.h"	              /* P7_HMMCACHE : cached profile DB (such as Pfam)           */
#include "daemon/cachedb.h"		      /* P7_SEQCACHE : cached sequence database (such as UniProt) */
#include "daemon/hmmdutils.h"		      /* support utilities for the daemon and client              */
#include "daemon/hmmdwrkr.h"		      /* worker process on a node                                 */
#include "daemon/hmmdmstr.h"		      /* master process on the daemon's head node                 */
#include "daemon/hmmpgmd2msa.h"		      /* utility for efficient MSA creation (Rob Finn)            */



/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
