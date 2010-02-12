#! /usr/bin/perl

# Test hmmconvert.
#
# Adapted from Exercises.sh, xref SRE, Fri Oct 23 10:38:44 1998
# conversion bug detected in 2.1 and fixed in 2.1.1a 
#
# CVS $Id$

$usage = "test1-conversion.pl <hmmconvert>\n";
if ($#ARGV != 0) { die "Wrong argument number.\n$usage"; }

$hmmconvert = shift;
$ok      = 1;

if ($ok) {
    system("$hmmconvert -F fn3-bin test1.tmp1 > /dev/null 2> /dev/null");
    if ($? != 0) { $ok = 0; }
}
if ($ok) {
    system("$hmmconvert -F fn3-bin-swap test1.tmp2 > /dev/null 2> /dev/null");
    if ($? != 0) { $ok = 0; }
}
if ($ok) {
    system("diff test1.tmp1 test1.tmp2 > /dev/null");
    if ($? != 0) { $ok = 0; }
}

foreach $tmpfile ("test1.tmp1", "test1.tmp2") {
    unlink $tmpfile if -e $tmpfile;
}

if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }



