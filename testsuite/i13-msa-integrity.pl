#! /usr/bin/perl

# Look for any problems in hmmalign that corrupt the input sequences.
#
# Usage:   ./i13-msa-integrity.pl <builddir> <srcdir> <tmpfile prefix>
# Example: ./i13-msa-integrity.pl ..         ..       tmpfoo
#
# SRE, Tue Mar  9 09:19:22 2010 [Janelia]
# SVN $Id$

BEGIN {
    $builddir  = shift;
    $srcdir    = shift;
    $tmppfx    = shift;
}

# Verify that we have all the executables we need for the test.
if (! -x "$builddir/src/hmmalign")                  { die "FAIL: didn't find hmmalign binary in $builddir/src\n";  }
if (! -x "$builddir/src/hmmemit")                   { die "FAIL: didn't find hmmemit binary in $builddir/src\n";  }
if (! -x "$builddir/easel/miniapps/esl-reformat")   { die "FAIL: didn't find esl-reformat binary in $builddir/easel/miniapps\n";  }
if (! -x "$builddir/easel/miniapps/esl-shuffle")    { die "FAIL: didn't find esl-reformat binary in $builddir/easel/miniapps\n";  }

# Verify that we have all the datafiles we need.
if (! -e "$srcdir/testsuite/RRM_1.hmm")  { die "FAIL: didn't find RRM_1.hmm in $srcdir/testsuite\n";  }
$profile = "$srcdir/testsuite/RRM_1.hmm";

foreach $trial (1..10)
{
    foreach $n (1, 2, 10, 1000)
    {
	# homologous sequence fragments: generated from local profile
	`$builddir/src/hmmemit -o $tmppfx.fa -N $n -L 0 -p --local $profile`;
	if ($? != 0) { die "FAIL: hmmemit failed\n"; }

	&msa_integrity_check("$tmppfx.fa", $profile);

	# random sequences
	`$builddir/easel/miniapps/esl-shuffle -G -N $n -L 50 --amino -o $tmppfx.fa`;
	if ($? != 0) { die "FAIL: esl-shuffle failed\n"; }

	&msa_integrity_check("$tmppfx.fa", $profile);
    }
}


print "ok\n";
unlink "$tmppfx.sto";
unlink "$tmppfx.fa";
unlink "$tmppfx.fa1";
unlink "$tmppfx.fa2";
exit 0;



sub msa_integrity_check
{
    my ($fafile, $hmmfile) = @_;

    `$builddir/src/hmmalign -o $tmppfx.sto $hmmfile $fafile 2>&1 > /dev/null`;
    if ($? != 0) { die "FAIL: hmmalign failed\n"; }
    
    `$builddir/easel/miniapps/esl-reformat -u fasta $tmppfx.sto 2>&1 > $tmppfx.fa1`;
    if ($? != 0) { die "FAIL: first esl-reformat failed\n"; }
    
    `$builddir/easel/miniapps/esl-reformat -u fasta $fafile     2>&1 > $tmppfx.fa2`;
    if ($? != 0) { die "FAIL: second esl-reformat failed\n"; }

    `diff -b $tmppfx.fa1 $tmppfx.fa2 2>&1 > /dev/null`;
    if ($? != 0) { die "FAIL: alignment corrupted\n"; }
}




