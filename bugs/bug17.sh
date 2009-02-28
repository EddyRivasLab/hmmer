#! /bin/sh

# Exercises and tests bug #17. 
# CVS $Id$
../binaries/shuffle -i -n 100 -t 200 --amino > bug17.fa
../binaries/hmmsearch ../testsuite/fn3.hmm bug17.fa > bug17.out1
../binaries/hmmsearch --cpu 0 ../testsuite/fn3.hmm bug17.fa > bug17.out2
diff bug17.out1 bug17.out2 
if (test $? != 0) then
  echo "FAILED: threaded/nonthreaded hmmsearch identity test (bug 17)"
fi
rm bug17.fa bug17.out1 bug17.out2