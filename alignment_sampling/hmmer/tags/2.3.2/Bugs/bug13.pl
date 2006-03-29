#! /usr/local/bin/perl

# Any program that reads multiple HMMs should refuse to 
# deal with "mixed" databases of DNA and protein HMMs,
# because behavior will be undefined. 
#
# First detected as an hmmconvert problem (by GCG's Christiane
# VanSchlun), and it's easiest to test there too. The other
# executables are not tested because they share the
# same mechanism (mixed HMMs are detected in hmmio.c).
#

$usage = "Usage: ./bug13.pl [BINARYPATH]\n";

$binpath = "";
if    ($#ARGV > 0)  { print $usage; exit 1; }
elsif ($#ARGV == 0) { $binpath = shift;     }

$hmmconvert = "$binpath/hmmconvert";
print "Bug 13 ...\t";

# These two commands should fail with non-zero exit status.
#
$output = `$hmmconvert -F bug13-mixed1.hmm bug13.tmp 2>/dev/null`;
if ($? == 0) { print "DETECTED! (1)\n"; unlink "bug13.tmp"; exit 1; }
$output = `$hmmconvert -F bug13-mixed2.hmm bug13.tmp 2>/dev/null`;
if ($? == 0) { print "DETECTED! (2)\n"; unlink "bug13.tmp"; exit 1; }

print "not detected.\n";
unlink "bug13.tmp";
exit 0;
