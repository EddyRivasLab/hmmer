#! /bin/sh

# Usage:
#   ./i4-zerolength-seqs.sh <hmmscan binary> <hmmsearch> <phmmer> <jackhmmer> <HMM database> <tmpfile prefix>
# 
# Example:
#   ../src/hmmbuild minifam.hmm minifam
#   ../src/hmmpress minifam.hmm
#   ./i4-zerolength-seqs.sh ../src/hmmscan ../src/hmmsearch ../src/phmmer ../src/jackhmmer minifam.hmm baz
#   rm minifam.hmm*

# Verifies that the search programs can take zero length sequences as
# input without crashing.
# 
# The <HMM database> must be press'ed, for hmmscan to work on it.

if test ! $# -eq 6; then 
  echo "Usage: $0 <hmmscan binary> <hmmsearch> <phmmer> <jackhmmer> <HMM database> <tmpfile prefix>"
  exit 1
fi

hmmscan=$1;   if test ! -x $hmmscan;   then echo "FAIL: $hmmbuild not executable";  exit 1; fi
hmmsearch=$2; if test ! -x $hmmsearch; then echo "FAIL: $hmmsearch not executable"; exit 1; fi
phmmer=$3;    if test ! -x $phmmer;    then echo "FAIL: $phmmer not executable";    exit 1; fi
jackhmmer=$4; if test ! -x $jackhmmer; then echo "FAIL: $jackhmmer not executable"; exit 1; fi
hmmfile=$5;   if test ! -r $hmmfile;   then echo "FAIL: $hmmfile not readable";     exit 1; fi

fafile=$6.fa

cat > $fafile <<EOF
>foo
>bar
YYYYY
EOF

$hmmscan   $hmmfile $fafile > /dev/null 2>&1; if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi
$hmmsearch $hmmfile $fafile > /dev/null 2>&1; if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi
$phmmer    $fafile $fafile  > /dev/null 2>&1; if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi
$jackhmmer $fafile $fafile  > /dev/null 2>&1; if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi

echo "ok"

rm $fafile
exit 0

