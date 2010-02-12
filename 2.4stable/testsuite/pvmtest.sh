#!/bin/sh

PVM_ARCH=LINUX
PVM_ROOT=/usr/seshare/Linux/pvm3
PVM_RSH=/usr/bin/ssh
PVM_EXPORT=HMMERDB
HMMERDB=/nfs/wol2/people/eddy/src/hmmer/testsuite/
export PVM_ARCH PVM_ROOT PVM_RSH PVM_EXPORT HMMERDB
echo halt | pvm 
echo quit | pvm pvm.conf
../binaries/hmmindex fn3.hmm
../binaries/hmmpfam --pvm fn3.hmm titin.fa
rm fn3.hmm.ssi