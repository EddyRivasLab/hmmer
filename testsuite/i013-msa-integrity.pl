#! /usr/bin/perl

# Look for any problems in hmmalign that corrupt the input sequences.
#
# Usage:   ./i013-msa-integrity.pl <builddir> <srcdir> <tmpfile prefix>
# Example: ./i013-msa-integrity.pl ..         ..       tmpfoo

$builddir  = shift;
$srcdir    = shift;
$tmppfx    = shift;

# Verify that we have all the executables we need for the test.
if (! -x "$builddir/src/programs/hmmalign")             { die "FAIL: didn't find hmmalign binary in $builddir/src/programs";  }
if (! -x "$builddir/src/programs/hmmemit")              { die "FAIL: didn't find hmmemit binary in $builddir/src/programs";  }
if (! -x "$builddir/src/programs/hmmsearch")            { die "FAIL: didn't find hmmsearch binary in $builddir/src/programs";  }
if (! -x "$builddir/lib/easel/miniapps/esl-reformat")   { die "FAIL: didn't find esl-reformat binary in $builddir/lib/easel/miniapps";  }
if (! -x "$builddir/lib/easel/miniapps/esl-shuffle")    { die "FAIL: didn't find esl-shuffle binary in $builddir/lib/easel/miniapps";  }

# Verify that we have all the datafiles we need.
if (! -e "$srcdir/testsuite/RRM_1.hmm")  { die "FAIL: didn't find RRM_1.hmm in $srcdir/testsuite";  }
$profile = "$srcdir/testsuite/RRM_1.hmm";

foreach $trial (1..5)
{
    foreach $n (1, 10, 100)
    {

	# homologous sequence fragments: generated from local profile
	`$builddir/src/programs/hmmemit -o $tmppfx.fa -N $n -L 0 -p --unilocal $profile`;
	if ($? != 0) { die "FAIL: hmmemit"; }

	&hmmalign_msa_integrity_check ("$tmppfx.fa", $profile);
	&hmmsearch_msa_integrity_check("$tmppfx.fa", $profile);

	# random sequences
	`$builddir/lib/easel/miniapps/esl-shuffle -G -N $n -L 50 --amino -o $tmppfx.fa`;
	if ($? != 0) { die "FAIL: esl-shuffle"; }

	&hmmalign_msa_integrity_check ("$tmppfx.fa", $profile);
	&hmmsearch_msa_integrity_check("$tmppfx.fa", $profile);
    }
}


print "ok\n";
unlink "$tmppfx.sto";
unlink <$tmppfx.fa*>;
unlink "$tmppfx.dtbl";
exit 0;



sub hmmalign_msa_integrity_check
{
    my ($fafile, $hmmfile) = @_;

    `$builddir/src/programs/hmmalign -o $tmppfx.sto $hmmfile $fafile > /dev/null 2>&1`;
    if ($? != 0) { die "FAIL: hmmalign failed"; }
    
    `$builddir/lib/easel/miniapps/esl-reformat -u fasta $tmppfx.sto > $tmppfx.fa1 2>/dev/null`;
    if ($? != 0) { die "FAIL: first esl-reformat failed"; }
    
    `$builddir/lib/easel/miniapps/esl-reformat -u fasta $fafile    > $tmppfx.fa2 2>/dev/null`;
    if ($? != 0) { die "FAIL: second esl-reformat failed"; }

    `diff -b $tmppfx.fa1 $tmppfx.fa2 > /dev/null 2>&1`;
    if ($? != 0) { die "FAIL: alignment corrupted\n"; }
    0;
}


sub hmmsearch_msa_integrity_check
{
    my ($fafile, $hmmfile) = @_;
    my $i;

    # need report = include threshold to have same hits in .sto, .domtbl
    `$builddir/src/programs/hmmsearch -E 0.01 --domE 0.01 -A $tmppfx.sto --domtbl $tmppfx.dtbl $hmmfile $fafile > /dev/null 2>&1`;
    if ($? != 0) { die "FAIL: hmmsearch failed"; }

    $nlines = `cat $tmppfx.dtbl | grep -v "^#" | wc -l`;
    if ($nlines == 0) { return 0; }

    `$builddir/lib/easel/miniapps/esl-sfetch --index $fafile > /dev/null 2>&1`;
    if ($? != 0) { die "FAIL: esl-sfetch --index failed"; }
    
    `cat $tmppfx.dtbl | grep -v "^#" | awk '{print \$1, \$18, \$19, \$1}' > $tmppfx.gdf`;
    if ($? != 0) { die "FAIL: gdf table"; }

    `$builddir/lib/easel/miniapps/esl-sfetch -Cf $fafile $tmppfx.gdf > $tmppfx.fa3`;
    if ($? != 0) { die "FAIL: esl-sfetch failed"; }

    # Grep out name/desc lines because they'll differ: esl-sfetch adds descline.
    `$builddir/lib/easel/miniapps/esl-reformat -u fasta $tmppfx.sto | grep -v "^>" > $tmppfx.fa1 2>/dev/null`;
    if ($? != 0) { die "FAIL: first esl-reformat failed"; }
    
    `$builddir/lib/easel/miniapps/esl-reformat -u fasta $tmppfx.fa3 | grep -v "^>" > $tmppfx.fa2 2>/dev/null`;
    if ($? != 0) { die "FAIL: second esl-reformat failed"; }

    `diff -b $tmppfx.fa1 $tmppfx.fa2 > /dev/null 2>&1`;
    if ($? != 0) { die "FAIL: alignment corrupted"; }

    0;
}




