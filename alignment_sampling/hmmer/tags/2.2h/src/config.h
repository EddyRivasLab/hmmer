/************************************************************
 * @LICENSE@
 ************************************************************/

/* config.h
 * 
 * Configurable compile-time parameters in HMMER.
 */

#ifndef CONFIGH_INCLUDED
#define CONFIGH_INCLUDED

/* RAMLIMIT determines the point at which we switch from fast,
 * full dynamic programming to slow, linear-memory divide and conquer
 * dynamic programming algorithms. It is the minimum amount of available
 * RAM on the systems the package will run on. It can be overridden
 * from the Makefile.
 * By default, we assume we have 32 Mb RAM available (per thread).
 */
#ifndef RAMLIMIT
#define RAMLIMIT 32
#endif

/* HMMER_NCPU determines the number of threads/processors that
 * a threads version will parallelize across. This can be overridden
 * by -DHMMER_NCPU=x in the Makefile, and by a setenv HMMER_NCPU x
 * in the environment, and usually by a command line option.
 * Usually we detect the number of processors dynamically, but
 * on some systems (FreeBSD and Linux, notably), we can't. On
 * these systems we assume 2 processors by default. That assumption
 * can be overridden here if HMMER_NCPU is uncommented.
 */
/* #define HMMER_NCPU 4 */ 

#define INTSCALE    1000.0      /* scaling constant for floats to integer scores   */
#define MAXABET     20	        /* maximum size of alphabet (4 or 20)              */
#define MAXCODE     23	        /* maximum degenerate alphabet size (17 or 23)     */
#define MAXDCHLET   200	        /* maximum # Dirichlet components in mixture prior */
#define NINPUTS     4	        /* number of inputs into structural prior          */
#define INFTY       987654321   /* infinity for purposes of integer DP cells       */
#define NXRAY       4           /* number of structural inputs                */
#define LOGSUM_TBL  20000       /* controls precision of ILogsum()            */
#define ALILENGTH   50		/* length of displayed alignment lines        */

#endif /*CONFIGH_INCLUDED*/

