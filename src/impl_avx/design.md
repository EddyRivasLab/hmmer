#Design notes for the AVX and AVX-512 implementation of HMMER#

## Overview ##
Our support for the newer SIMD instruction sets on Intel processors differs somewhat from the other vector
implementations in HMMER 3.  The src/impl_sse directory contains SSE implementations of all vectorized
functions, and is provided for compatibility with older processors and compilers.  This directory
(src/impl_avx) contains SSE, AVX, and AVX-512 implementations of our vectorized functions, as well as
code that selects the implementation (or multiple implementations, for testing) to use at runtime based on
the capabilities of the hardware HMMER is running on.

Selecting an x86 vector implementation is a two-phase process.  At compile time, the configure script tests
to see which of the SSE, AVX, and AVX-512  SIMD ISAs the compiler supports, and defines one or more of
eslENABLE_SSE, esl_ENABLE_AVX, and eslENABLE_AVX512 in src/p7_config.h and easel/esl_config.h,  These
defines control the compilation of all SIMD code, such that HMMER gets built with support for all of the
vector ISAs the compiler can handle.

Our vectorized functions are declared as function pointers that initially point to a "dispatcher" routine
instead of the actual function.  When called, the dispatcher routine queries the hardware to find out
what SIMD ISAs it supports, using functions defined in easel/esl_cpu.c, and sets the function pointer that
originally pointed to it to a function that uses the highest-performance ISA that the hardware supports.
For testing purposes, we also provide versions of most functions that call multiple SIMD implementations
and compare the results.

The file src/impl_avx/impl_avx.h contains declarations for our vector API, as well as definitions for
vectorized data structures.

Within src/impl_avx, most functions are distributed across four files, one (e.g. foo.c) that contains the
function pointer, the dispatch routines, any multi-ISA test functions, and our utest and benchmark
functions; a second (foo_sse.c) that contains the SSE implementation; and two that contain the AVX and
AVX-512 implementations (foo_avx.c and foo_avx512.c).  This makes it easier to ensure that no ISA-specific
code is compiled into a function for a different ISA and simplifies handling of static functions that are
only visible within a file.