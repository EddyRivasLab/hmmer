MODULE = dp_vector

LIBOBJS =  \
           fwdfilter.o\
           io.o\
	   msvfilter.o\
	   p7_checkptmx.o\
	   p7_filtermx.o\
	   p7_oprofile.o\
	   p7_oprofile_mpi.o\
	   simdvec.o\
	   ssvfilter.o\
	   vitfilter.o

SSELIBOBJS =  \
        io_sse.o \
        msvfilter_sse.o\
        p7_checkptmx_sse.o\
        p7_filtermx_sse.o\
        p7_oprofile_sse.o\
        ssvfilter_sse.o\
        vitfilter_sse.o\
        fwdfilter_sse.o 

AVXLIBOBJS = \
      io_avx.o \
       msvfilter_avx.o\
       p7_checkptmx_avx.o\
       p7_filtermx_avx.o\
       p7_oprofile_avx.o\
        ssvfilter_avx.o \
        vitfilter_avx.o \
        fwdfilter_avx.o

AVX512LIBOBJS = \
      io_avx512.o\
       msvfilter_avx512.o\
       p7_checkptmx_avx512.o\
       p7_filtermx_avx512.o\
       p7_oprofile_avx512.o\
        ssvfilter_avx512.o \
        vitfilter_avx512.o \
        fwdfilter_avx512.o

NEONLIBOBJS=  \
              io_neon.o \
              msvfilter_neon.o \
              p7_checkptmx_neon.o \
    
LIBHDRS =  

OBJS =

HDRS =

PROGS = 

UTESTS = fwdfilter_utest \
         io_utest \
         msvfilter_utest \
	 p7_oprofile_mpi_utest \
         vitfilter_utest

STATS  = fwdfilter_stats

BENCHMARKS = fwdfilter_benchmark \
             io_benchmark \
             msvfilter_benchmark  \
             p7_oprofile_benchmark \
             vitfilter_benchmark

EXAMPLES = fwdfilter_example \
           io_example \
           msvfilter_example \
           p7_oprofile_example \
           vitfilter_example

LIBTMP =

top_srcdir   = ../..
top_builddir = ../..
srcdir       = .

SHELL        = /bin/sh

include ../Makefile-subdirs.mk

################################################################
# @LICENSE@
#
# SVN $URL$
# SVN $Id$
################################################################

