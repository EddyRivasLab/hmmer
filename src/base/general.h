#ifndef P7_GENERAL_INCLUDED
#define P7_GENERAL_INCLUDED

#include "p7_config.h"
#include "esl_getopts.h"

/* Search modes. */
#define p7_NO_MODE   0
#define p7_LOCAL     1		/* multihit local:  "fs" mode   */
#define p7_GLOCAL    2		/* multihit glocal: "ls" mode   */
#define p7_UNILOCAL  3		/* unihit local: "sw" mode      */
#define p7_UNIGLOCAL 4		/* unihit glocal: "s" mode      */
#define p7_MIXED     5          /* mixed mode, local+glocal     */
#define p7_UNIMIXED  

#define p7_IsLocal(mode)  (mode == p7_LOCAL || mode == p7_UNILOCAL)
#define p7_IsMulti(mode)  (mode == p7_LOCAL || mode == p7_GLOCAL)

#define p7_NEVPARAM 6	/* number of statistical parameters stored in models                      */
#define p7_NCUTOFFS 6	/* number of Pfam score cutoffs stored in models                          */
#define p7_NOFFSETS 3	/* number of disk offsets stored in models for hmmscan's fast model input */
enum p7_evparams_e {    p7_MMU  = 0, p7_MLAMBDA = 1,     p7_VMU = 2,  p7_VLAMBDA = 3, p7_FTAU = 4, p7_FLAMBDA = 5 };
enum p7_cutoffs_e  {     p7_GA1 = 0,     p7_GA2 = 1,     p7_TC1 = 2,      p7_TC2 = 3,  p7_NC1 = 4,     p7_NC2 = 5 };
enum p7_offsets_e  { p7_MOFFSET = 0, p7_FOFFSET = 1, p7_POFFSET = 2 };

#define p7_EVPARAM_UNSET -99999.0f  /* if evparam[0] is unset, then all unset                         */
#define p7_CUTOFF_UNSET  -99999.0f  /* if cutoff[XX1] is unset, then cutoff[XX2] unset, XX={GA,TC,NC} */
#define p7_COMPO_UNSET   -1.0f      /* if compo[0] is unset, then all unset                           */

/* Option flags when creating multiple alignments with p7_tracealign_*() */
#define p7_DEFAULT             0
#define p7_DIGITIZE            (1<<0)
#define p7_ALL_CONSENSUS_COLS  (1<<1)
#define p7_TRIM                (1<<2)

/* Option flags when creating faux traces with p7_trace_FauxFromMSA() */
#define p7_MSA_COORDS	       (1<<0) /* default: i = unaligned seq residue coords     */

/* MEA traceback routines have to check consistency, all path transitions must be nonzero */
#define P7_DELTAT(val, tsc) ( ((tsc) == -eslINFINITY) ? -eslINFINITY : (val))

#define p7O_EXTRA_SB 17    /* see ssvfilter.c for explanation */

/* Which strand(s) should be searched */
enum p7_strands_e {    p7_STRAND_TOPONLY  = 0, p7_STRAND_BOTTOMONLY = 1,  p7_STRAND_BOTH = 2};

/* Macros below implement indexing idioms for generic DP routines.
 * They require the following setup, for profile <gm> and matrix <gx>:
 *   float const *tsc = gm->tsc;
 *   float      **dp  = gx->dp;
 *   float       *xmx = gx->xmx;
 * and for row i, i+1 (target residues x_i, x_i+1 in digital seq <dsq>):
 *   float const *rsc = gm->rsc[dsq[i]];
 *   float const *rsn = gm->rsc[dsq[i+1]];
 */
#define MMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_M])
#define IMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_I])
#define DMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_D])
#define XMX(i,s) (xmx[(i) * p7G_NXCELLS + (s)])

#define TSC(s,k) (tsc[(k) * p7P_NTRANS + (s)])
#define MSC(k)   (rsc[(k) * p7P_NR     + p7P_M])
#define ISC(k)   (rsc[(k) * p7P_NR     + p7P_I])
#define MSN(k)   (rsn[(k) * p7P_NR     + p7P_M])
#define ISN(k)   (rsn[(k) * p7P_NR     + p7P_I])

/* Flags that control P7_GMX debugging dumps */
#define p7_HIDE_SPECIALS (1<<0)
#define p7_SHOW_LOG      (1<<1)


extern int          p7_Init(void);
extern void         p7_banner(FILE *fp, char *progname, char *banner);
extern void         p7_Die (char *format, ...);
extern void         p7_Fail(char *format, ...);
extern ESL_GETOPTS *p7_CreateDefaultApp(ESL_OPTIONS *options, int nargs, int argc, char **argv, char *banner, char *usage);
extern int          p7_AminoFrequencies(float *f);

#endif /*P7_GENERAL_INCLUDED*/

/************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 ************************************************************/
