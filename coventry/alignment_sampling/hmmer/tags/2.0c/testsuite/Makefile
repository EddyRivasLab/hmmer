# Makefile for Shiva testsuite in HMMER
# RCS $Id$
##########
# HMMER - Biological sequence analysis with HMMs
# Copyright (C) 1992-1998 Sean R. Eddy
#
#   This source code is distributed under the terms of the 
#   GNU General Public License. See the files COPYING and 
#   LICENSE for details.
#    
###########

## your compiler and flags
#
CC     = cc
CFLAGS = -g -Wall

## archiving command, and ranlib command if you need one
#  if you need ranlib and it's in /bin or /usr/bin, you
#  don't need to change RANLIB. If you don't need ranlib,
#  you also don't need to change anything.
#
AR     = ar rcv 
RANLIB = ranlib

#######
## should not need to modify below this line
#######
SHELL  = /bin/sh
LIBS   = -lhmmer -lm

SHIVA = evd_test masks_test tophits_test trace_test

SHIVASRC  = evd_test.c masks_test.c tophits_test.c trace_test.c

SHIVADATA = trace_test.hmm trace_test.seq trace_test.slx 

.c.o:
	$(CC) $(CFLAGS) -I../src/ -c $<		

#######
## Targets defining how to make Shiva executables.
#######

all: 	$(SHIVA)

evd_test: $(OBJ) evd_test.o
	$(CC) $(CFLAGS) -o $@ -L../src/ evd_test.o $(LIBS)

masks_test: $(OBJ) masks_test.o
	$(CC) $(CFLAGS) -o $@ -L../src/ masks_test.o $(LIBS)

tophits_test: $(OBJ) tophits_test.o
	$(CC) $(CFLAGS) -o $@ -L../src/ tophits_test.o $(LIBS)

trace_test: $(OBJ) trace_test.o
	$(CC) $(CFLAGS) -o $@ -L../src/ trace_test.o $(LIBS)

#######
## `make shiva` runs the tests.
#######

shiva:   $(SHIVA)
	@for shivatest in $(SHIVA); do\
	   $$shivatest;\
	   echo $$shivatest: ok;\
	done

#######
## Miscellaneous
#######

clean:
	-rm -f *.o *~ Makefile.bak core $(SHIVA) TAGS gmon.out

tags:
	etags -t *.c

