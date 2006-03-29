/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1997 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the 
 *   GNU General Public License. See the files COPYING and 
 *   GNULICENSE for details.
 *    
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
 * from the Makefile. By default, we assume we have 32 Mb RAM available.
 */
#ifndef RAMLIMIT
#define RAMLIMIT 32
#endif

#define INTSCALE    1000.0      /* scaling constant for floats to integer scores   */
#define MAXABET     20	        /* maximum size of alphabet (4 or 20)              */
#define MAXCODE     23	        /* maximum degenerate alphabet size (17 or 23)     */
#define MAXDCHLET   200	        /* maximum # Dirichlet components in mixture prior */
#define NINPUTS     4	        /* number of inputs into structural prior          */
#define INFTY       987654321   /* infinity for purposes of integer DP cells       */
#define NXRAY       4           /* number of structural inputs                */
#define LOGSUM_TBL  20000       /* controls precision of Logsum()             */
#define ALILENGTH   50		/* length of displayed alignment lines        */

/* Debugging levels in HMMER are controlled by conditional preprocessing
 */
#ifndef DEBUGLEVEL
#define DEBUGLEVEL 0
#endif
#define DEBUG_NONE        0	/* no debugging output                          */
#define DEBUG_LIGHT       1	/* turn on most assertions                      */
#define DEBUG_SOME        2	/* turn on some output; low verbosity           */
#define DEBUG_AGGRESSIVE  3	/* turn on all assertions; mild verbosity       */
#define DEBUG_LOTS        4	/* turn on most output; moderate/high verbosity */
#define DEBUG_EXTREME     5	/* turn on all output; intolerable verbosity    */


#endif /*CONFIGH_INCLUDED*/

