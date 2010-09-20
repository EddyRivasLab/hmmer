/* The all-encompassing include file for HMMER.
 * All-encompassing because there's a lot of crossdependency.
 * There's some opportunity for modularity, but not a lot.
 *
 *    1. P7_HMM:         a core model.
 *    2. Copyright and license information.
 *   
 * SRE, Wed Jan  3 13:46:42 2007 [Janelia]
 * SVN $Id: hmmer.h 3350 2010-08-19 17:36:13Z farrarm $
 */
#ifndef P7_HMMPGMD_INCLUDED
#define P7_HMMPGMD_INCLUDED

/*****************************************************************
 * 1.
 *****************************************************************/

typedef struct {
  uint32_t   status;            /* error status                             */
  uint32_t   err_len;           /* if status not zero, the stream will next */
                                /* contain the error message.               */
} HMMER_SEARCH_STATUS;

typedef struct {
  double     elapsed;         	/* elapsed time, seconds                    */
  double     user;            	/* CPU time, seconds                        */
  double     sys;             	/* system time, seconds                     */

  double  Z;			/* eff # targs searched (per-target E-val)  */
  double  domZ;			/* eff # signific targs (per-domain E-val)  */
  enum p7_zsetby_e Z_setby;   	/* how Z was set                            */
  enum p7_zsetby_e domZ_setby;	/* how domZ was set                         */

  uint64_t   nmodels;         	/* # of HMMs searched                       */
  uint64_t   nseqs;           	/* # of sequences searched                  */
  uint64_t   nres;            	/* # of residues searched                   */
  uint64_t   nnodes;          	/* # of model nodes searched                */
  uint64_t   n_past_msv;      	/* # comparisons that pass MSVFilter()      */
  uint64_t   n_past_bias;     	/* # comparisons that pass bias filter      */
  uint64_t   n_past_vit;      	/* # comparisons that pass ViterbiFilter()  */
  uint64_t   n_past_fwd;      	/* # comparisons that pass ForwardFilter()  */

  uint64_t   nhits;           	/* number of hits in list now               */
  uint64_t   nreported;       	/* number of hits that are reportable       */
  uint64_t   nincluded;       	/* number of hits that are includable       */
} HMMER_SEARCH_STATS;

#endif /*P7_HMMPGMD_INCLUDED*/

size_t writen(int fd, const void *vptr, size_t n);
size_t readn(int fd, void *vptr, size_t n);

/************************************************************
 * @LICENSE@
 ************************************************************/
