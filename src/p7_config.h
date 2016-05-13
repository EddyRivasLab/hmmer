/* @configure_input@  
 * p7_config.h is generated from p7_config.h.in by the ./configure script.
 * DO NOT EDIT p7_config.h; only edit p7_config.h.in.
 *
 *
 * Configuration of HMMER, including both system-dependent
 * configuration (done by ./configure) and hardcoded configuration
 * that someone might want to alter manually.
 *
 * Because this header may configure the behavior of system headers
 * (for example, LFS support), it must be included before any other
 * header file.
 */
#ifndef P7_CONFIGH_INCLUDED
#define P7_CONFIGH_INCLUDED


/*****************************************************************
 * 1. Default parameters.
 *****************************************************************/

/* Model parameterization, relative entropy target defaults, empirically tuned, in bits: */
#define p7_ETARGET_AMINO         0.59  //  .. for protein; from work of Steve Johnson. 
#define p7_ETARGET_DNA           0.45  //  .. for nucleic; from work of Travis Wheeler. 
#define p7_ETARGET_OTHER          1.0  //  .. for other (custom) alphabets.

/* Sparsification in checkpointed/vectorized local decoding: fwdfilter.c */
#define p7_SPARSIFY_RAMLIMIT      128  // Memory "redline" cap on the O(M sqrt L) checkpoint mx
#define p7_SPARSIFY_THRESH       0.01  // per-cell posterior probability inclusion threshold 

/* MPAS algorithm: {reference,sparse}_anchors.c */
#define p7_MPAS_LOSS_THRESHOLD  0.001  // Controls main convergence criterion
#define p7_MPAS_MAX_ITERATIONS   1000  // Limits to maximum # iterations when unconverged
#define p7_MPAS_NMAX_SAMPLING   FALSE  // TRUE turns off conv. criteria, always samples to max #
#define p7_MPAS_BE_VERBOSE      FALSE  // TRUE turns on internal printf info dumping

/* Comparison engine */
#define p7_ENGINE_FIXED_SEED       42  // if 0, RNG is seeded randomly
#define p7_ENGINE_REPRODUCIBLE   TRUE  // TRUE reseeds RNG for every comparison, making results order-independent
#define p7_ENGINE_DO_BIASFILTER  TRUE  // Use ad hoc "bias filter" after MSV/SSV step

#define p7_SEQDBENV          "BLASTDB"
#define p7_HMMDBENV          "PFAMDB"

#define p7_ALILENGTH   50   // Length of alignment output lines, in characters

/* Controls block_size when reading seq data in nhmmer */
#define p7_NHMMER_MAX_RESIDUE_COUNT (1024 * 256) /* 0.25 MB */



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

/* Other stuff.
 */
#define p7_MAXDCHLET  20      /* maximum # Dirichlet components in mixture prior */


/*****************************************************************
 * 4. The following constants define our SIMD vector layout and memory
 *    alignment.  Although SSE, Altivec/VMX are 128b/16B vectors, we
 *    must anticipate different vector sizes.  For example, Intel AVX
 *    is already roadmapped out to 1024b/128B vectors.  See note [1]
 *    in src/impl_{sse/vmx}.h on memory alignment, SIMD vectors, and
 *    malloc(). We need these constants and macros not only in the
 *    vector implementation, but also in the P7_SPARSEMASK code that
 *    interfaces with the vector f/b filter, which is why we need
 *    these constants in the general config.h.
 *****************************************************************/

#define p7_VALIGN   16		/* Vector memory must be aligned on 16-byte boundaries   */
#define p7_VNF      4		/* Number of floats per SIMD vector (Forward, Backward)  */
#define p7_VNW      8		/* Number of shorts (words) per SIMD vector (Viterbi)    */
#define p7_VNB      16		/* Number of bytes per SIMD vector (SSV, MSV)            */
#define p7_VALIMASK (~0xf)      /* Ptrs are aligned using & p7_VALIMASK                  */




/*****************************************************************
 * 5. The final section isn't meant to be human editable at all.
 *    It is configured automatically by the ./configure script. 
 *****************************************************************/

/* Version info - set once for whole package in configure.ac
 */
#define HMMER_VERSION "4.0dev"
#define HMMER_DATE "Apr 2015"
#define HMMER_COPYRIGHT "Copyright (C) 2015 Howard Hughes Medical Institute."
#define HMMER_LICENSE "Freely distributed under the open source BSD license."
#define HMMER_URL "hmmer.org"

/* Large file support (must precede any header file inclusion.)
 */
/* #undef _FILE_OFFSET_BITS */
/* #undef _LARGE_FILES */
/* #undef _LARGEFILE_SOURCE */

/* System headers
 */
#define HAVE_STRINGS_H 1

/* #undef HAVE_ENDIAN_H */
#define HAVE_INTTYPES_H 1
#define HAVE_STDINT_H 1
/* #undef HAVE_UNISTD_HS */
#define HAVE_SYS_TYPES_H 1
#define HAVE_NETINET_IN_H 1

#define HAVE_SYS_PARAM_H 1
#define HAVE_SYS_SYSCTL_H 1

#define HAVE_XMMINTRIN_H 1
#define HAVE_EMMINTRIN_H 1
#define HAVE_PMMINTRIN_H 1

/* Optional parallel implementations
 */
#define HAVE_SSE2 1
#define HAVE_AVX2 1
/* #undef HAVE_MPI */
/* #undef HMMER_PVM */
#define HMMER_THREADS 1
/* #undef HAVE_PTHREAD_ATTR_SETSCOPE */
/* #undef HAVE_PTHREAD_SETCONCURRENCY */

/* Which SIMD code should we build?
HAVE_x macros define what the compiler supports.
p7_build_x defines what should get built.
*/
#define p7_build_SSE 1
#define p7_build_AVX2 1
/* #undef p7_build_AVX512 */
#define p7_build_check_AVX2 1
/* #undef p7_build_check_AVX512 */

/* Optional processor specific support
 */
#define HAVE_FLUSH_ZERO_MODE 1

/* Debugging and development hooks
 */
#define p7_DEVELOPMENT  1            // TRUE when code is in development. Allow things like stochastic test failures that would frighten civilians.
/* #undef p7_DEBUGGING */
/* #undef p7_LOGSUM_SLOWEXACT */
#define HAVE_FUNC_ATTRIBUTE_NORETURN 1

#endif /*P7_CONFIGH_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
