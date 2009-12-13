#! /usr/bin/perl

# Check that we can deal with profiles and sequences that contain
# duplicate names, both as queries and targets. 
#
# Usage:    ./i10-duplicate-names.pl <bindir> <tmpfile prefix>
# Example:  ./i10-duplicate-names.pl ../src   tmpfoo
#
# SRE, Sun Dec 13 14:41:31 2009 [Yokohama, Japan]
# SVN $Id$

use h3;

$bindir = shift;
$tmppfx = shift;

# Verify that we have all the executables we need for the test.
if (! -x "$bindir/hmmbuild")   { die "FAIL: didn't find hmmbuild binary in $bindir\n";  }
if (! -x "$bindir/hmmpress")   { die "FAIL: didn't find hmmpress binary in $bindir\n";  }
if (! -x "$bindir/hmmsearch")  { die "FAIL: didn't find hmmsearch binary in $bindir\n"; }
if (! -x "$bindir/hmmscan")    { die "FAIL: didn't find hmmscan binary in $bindir\n";   }
if (! -x "$bindir/phmmer")     { die "FAIL: didn't find phmmer binary in $bindir\n"; }
if (! -x "$bindir/jackhmmer")  { die "FAIL: didn't find jackhmmer binary in $bindir\n";   }

# Create our test files
if (! open(ALI1, ">$tmppfx.sto")) { print "FAIL: couldn't open $tmppfx.sto for write";  exit 1; }
if (! open(SEQ1, ">$tmppfx.fa"))  { print "FAIL: couldn't open $tmppfx.fa for write";   exit 1; }

print ALI1 <<"EOF";
# STOCKHOLM 1.0
#=GF ID profile
#=GF AC XX01234.5
#=GF DE A test description
seq1 ACDEFGHIKLMNPQRSTVWY
seq2 ACDEFGHIKLMNPQRSTVWY
seq3 ACDEFGHIKLMNPQRSTVWY
//
# STOCKHOLM 1.0
#=GF ID profile
#=GF AC XX01234.5
#=GF DE A test description
seq1 ACDEFGHIKLLMNPQRSTVWY
seq2 ACDEFGHIKLLMNPQRSTVWY
seq3 ACDEFGHIKLLMNPQRSTVWY
//
EOF

print SEQ1 << "EOF";
>seq
ACDEFGHIKLMNPQRSTVWY
>seq
ACDEFGHIKLLMNPQRSTVWY
EOF

close ALI1;
close SEQ1;

# Build profiles from the test alignments
@output = `$bindir/hmmbuild $tmppfx.hmm $tmppfx.sto 2>&1`;
if ($? != 0) { die "FAIL: hmmbuild failed\n"; }
@output = `$bindir/hmmpress $tmppfx.hmm             2>&1`;
if ($? != 0) { die "FAIL: hmmpress failed\n"; }



# phmmer should show four results
$output = `$bindir/phmmer --tblout $tmppfx.tbl $tmppfx.fa $tmppfx.fa 2>&1`;
if ($? != 0) { die "FAIL: phmmer failed\n"; }

&h3::ParseTbl("$tmppfx.tbl");
if ($h3::ntbl != 4) { die "FAIL: on expected number of hits, phmmer\n"; } 


# jackhmmer should show four results
$output = `$bindir/jackhmmer --tblout $tmppfx.tbl $tmppfx.fa $tmppfx.fa 2>&1`;
if ($? != 0) { die "FAIL: jackhmmer failed\n"; }

&h3::ParseTbl("$tmppfx.tbl");
if ($h3::ntbl != 4) { die "FAIL: on expected number of hits, jackhmmer\n"; } 


# hmmsearch should show four results
$output = `$bindir/hmmsearch --tblout $tmppfx.tbl $tmppfx.hmm $tmppfx.fa 2>&1`;
if ($? != 0) { die "FAIL: hmmsearch failed\n"; }

&h3::ParseTbl("$tmppfx.tbl");
if ($h3::ntbl != 4) { die "FAIL: on expected number of hits, hmmsearch\n"; } 



# hmmscan should show four results
$output = `$bindir/hmmscan --tblout $tmppfx.tbl $tmppfx.hmm $tmppfx.fa 2>&1`;
if ($? != 0) { die "FAIL: hmmscan failed\n"; }

&h3::ParseTbl("$tmppfx.tbl");
if ($h3::ntbl != 4) { die "FAIL: on expected number of hits, hmmscan\n"; } 

print "ok\n";
unlink "$tmppfx.sto";
unlink "$tmppfx.fa";
unlink "$tmppfx.tbl";
unlink <$tmppfx.hmm*>;
exit 0;





