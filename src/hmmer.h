/* hmmer.h
 * Most of the compile-time configuration constants for HMMER.
 * 
 * Also see p7config.h, for system dependencies that are automatically
 * determined by the configure script.
 * 
 * SRE, Thu Apr  6 13:18:13 2006 [AA890 enroute to Boston]
 * SVN $Id$
 */


/* The RAMLIMIT determines the point at which we switch from fast,
 * full dynamic programming to slow, linear-memory divide and conquer
 * dynamic programming algorithms. It is the minimum amount of available
 * RAM on the systems the package will run on. It can be overridden
 * from the Makefile (hence the #ifdef).
 * By default, we assume we have 32 Mb RAM available (per thread).
 */
#ifndef p7_RAMLIMIT
#define p7_RAMLIMIT        32
#endif





/* Scores are kept internally as scaled integer log-odds scores:
 *     sc_x = INTSCALE * log_2 p_x / f_x, rounded off.
 * 
 * We use 32-bit signed integers to hold scores, giving a dynamic
 * range of -2.147e9 to 2.147e9. With INTSCALE at its default 
 * 1000, this means we hold scores of -2.1e6..2.1e6 bits, with
 * three digits of fractional precision.
 * 
 * We therefore need at least 11 characters width to guarantee
 * space for a raw, printed integer score; thus TXTLEN. See
 * p7_Score2ASCII() for use.
 * 
 * To deal with numerical underflow (in DP calculations, for example),
 * we define a threshold P7_IMPOSSIBLE for -infinity. We have to be
 * able to add two -infinities together without underflow, so we can
 * check the threshold and reset a result of an addition to our
 * internal -infinity.
 * 
 * In principle, one might attempt to store scores in less space
 * (particularly handy for taking advantage of fast graphics hardware,
 * which tend to provide vector instructions for 8-bit and 16-bit
 * quantities, but not always for 32-bit quantities).  16 bits only
 * gives us a range of +/-32767, though; HMMER scores range within
 * about +/-10000, and the three digits of fractional precision are
 * thought to be important, so we don't want to decrease INTSCALE.
 * It's just barely conceivable that you could use INTSCALE=100
 * and try to capture scores between -327.67 and 327.67; you'd have to deal 
 * with scores that underflow or overflow, though.
 * 
 */
#define p7_INTSCALE        1000.0     /* scaling for integer scores      */
#define p7_MAX_SC_TXTLEN   11	      /* maxlen of sc as string; inc \0  */
#define p7_IMPOSSIBLE     -987654321  /* an internal score of -infinity  */


/* HMMER_NCPU: 
 * If defined, this determines the number of threads for
 * multithreading.  It can be defined here, or by -DHMMER_NCPU=x in the
 * Makefile, or by a setenv HMMER_NCPU x in the environment, or by a
 * command line option to hmmpfam or hmmsearch. Default behaviour
 * (when NCPU is unset) is to autodetect the number of processors
 * dynamically, and use them all. On systems (FreeBSD and older Linuxen, notably)
 * where autodetection is 
 * can't autodetect the available # of cpus. On these systems we
 * assume 2 processors by default - dual processor Intel servers
 * are common. That assumption can be overridden
 * here if HMMER_NCPU is uncommented.  
 */
/* #define p7_NCPU 4 */ 
