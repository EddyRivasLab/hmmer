/* src/p7_config.h.  Generated from p7_config.h.in by configure.  */
/* @configure_input@
 * p7config.h.in -> p7config.h
 * 
 * p7config.h is generated from p7config.h.in by the ./configure script.
 * DO NOT EDIT p7config.h; only edit p7config.h.in.
 *
 * Configuration of HMMER, including both system-dependent configuration
 * (done by ./configure) and hardcoded configuration that someone might
 * want to alter someday.
 *
 * Because this header may configure the behavior of system headers
 * (for example, LFS support), it must be included before any other
 * header file.
 * 
 * SRE, Mon Jan  1 16:07:28 2007 [Casa de Gatos] [Nirvana, Nevermind]
 * SVN $Id: p7_config.h.in 3152 2010-02-07 22:55:22Z eddys $
 */
#ifndef P7_CONFIGH_INCLUDED
#define P7_CONFIGH_INCLUDED


/*****************************************************************
 * 1. Compile-time constants that control HMMER's computational 
 *    behavior (memory and processor use), and output formatting.
 *    It can be edited and configured manually before compilation.
 *****************************************************************/

/* p7_RAMLIMIT controls the switch from fast full DP to slow
 * linear-memory divide and conquer. Default devotes 32 MB/thread.
 */
#ifndef p7_RAMLIMIT
#define p7_RAMLIMIT   32
#endif

/* p7_ALILENGTH controls length of displayed alignment lines.
 */
#ifndef p7_ALILENGTH
#define p7_ALILENGTH       50
#endif

/*****************************************************************
 * 2. Compile-time constants that control empirically tuned HMMER
 *    default parameters. You can edit it, but you ought not to, 
 *    unless you're trying to improve on our empirical data.
 *****************************************************************/

/* Relative entropy target defaults:
 * For proteins, hmmbuild's effective sequence number calculation
 * aims to achieve a certain relative entropy per match emission.
 * (= average score per match emission).
 * These are empirically tuned constants, from the work of Steve Johnson.
 */
#define p7_ETARGET_AMINO  0.59 /* bits */
#define p7_ETARGET_DNA    0.59 /* bits */
#define p7_ETARGET_OTHER  1.0  /* bits */ /* if you define your own alphabet, set this */

#define p7_SEQDBENV          "BLASTDB"
#define p7_HMMDBENV          "PFAMDB"

/*****************************************************************
 * 3. The next section probably shouldn't be edited at all, unless
 *    you really know what you're doing. It controls some fundamental
 *    parameters in HMMER that occasionally get reconfigured in
 *    experimental versions, or for variants of HMMER that work on
 *    non-biological alphabets.
 *****************************************************************/

/* The symbol alphabet is handled by ESL_ALPHABET objects, which
 * dynamically allocate; but sometimes HMMER uses statically-allocated
 * space, and it's useful to know a reasonable maximum for
 * symbol alphabet size.
 */
#define p7_MAXABET    20      /* maximum size of alphabet (4 or 20)              */
#define p7_MAXCODE    29      /* maximum degenerate alphabet size (18 or 29)     */

/* p7_MAX_SC_TXTLEN has to be large enough to represent a score as a
 * string, including \0 and a sign.
 */
#define p7_MAX_SC_TXTLEN   11	      

/* In Forward algorithm implementations, we use a table lookup in
 * p7_FLogsum() to calculate summed probabilities in log
 * space. p7_INTSCALE defines the precision of the calculation; the
 * default of 1000.0 means rounding differences to the nearest 0.001
 * nat. p7_LOGSUM_TBL defines the size of the lookup table; the
 * default of 16000 means entries are calculated for differences of 0
 * to 16.000 nats (when p7_INTSCALE is 1000.0).  e^{-p7_LOGSUM_TBL /
 * p7_INTSCALE} should be on the order of the machine FLT_EPSILON,
 * typically 1.2e-7.
 */
#define p7_INTSCALE     1000.0f
#define p7_LOGSUM_TBL   16000

/* Other stuff.
 */
#define p7_MAXDCHLET  20      /* maximum # Dirichlet components in mixture prior */


/*****************************************************************
 * 4. The final section isn't meant to be human editable at all.
 *    It is configured automatically by the ./configure script. 
 *****************************************************************/

/* Version info - set once for whole package in configure.ac
 */
#define HMMER_VERSION "3.0rc1"
#define HMMER_DATE "February 2010"
#define HMMER_COPYRIGHT "Copyright (C) 2010 Howard Hughes Medical Institute."
#define HMMER_LICENSE "Freely distributed under the GNU General Public License (GPLv3)."
#define HMMER_URL "http://hmmer.org/"

/* Large file support (must precede any header file inclusion.)
 */
/* #undef _FILE_OFFSET_BITS */
/* #undef _LARGE_FILES */
/* #undef _LARGEFILE_SOURCE */

/* Choice of optimized implementation (one and only one must be set)
 */
#define p7_IMPL_SSE 1
/* #undef p7_IMPL_VMX */
/* #undef p7_IMPL_DUMMY */

/* Optional parallel implementations
 */
#define HAVE_SSE2 1
/* #undef HAVE_MPI */
/* #undef HMMER_PVM */
#define HMMER_THREADS 1
/* #undef HAVE_PTHREAD_ATTR_SETSCOPE */
/* #undef HAVE_PTHREAD_SETCONCURRENCY */

/* Optional processor specific support
 */
#define HAVE_FLUSH_ZERO_MODE 1

/* Debugging hooks
 */
/* #undef p7_DEBUGGING */

#endif /*P7_CONFIGH_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
