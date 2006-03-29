#! /bin/sh

# Various exercises to test the package.
# SRE, Fri Oct 23 10:38:44 1998
# RCS $Id$

# Test binary formats and interconversion.
#    (tests for bug detected in 2.1, fixed in 2.1.1a.)
#
../binaries/hmmconvert -F fn3-bin      ex1.tmp > /dev/null
../binaries/hmmconvert -F fn3-bin-swap ex2.tmp > /dev/null
diff ex1.tmp ex2.tmp > /dev/null
if (test $? != 0) then
  echo FAILED: hmmconvert byteswap test
fi
rm ex1.tmp ex2.tmp

