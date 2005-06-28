#! /usr/local/bin/perl

# hmmconvert should refuse to append ASCII HMMs to binary files,
# and vice versa.
#    bug14-ascii.hmm  - a small ascii HMM
#    bug14-binary.hmm - a small binary HMM

$usage = "Usage: ./bug14.pl [BINARYPATH]\n";

$binpath = "";
if    ($#ARGV > 0)  { print $usage; exit 1; }
elsif ($#ARGV == 0) { $binpath = shift;     }

$hmmconvert = "$binpath/hmmconvert";
print "Bug 14 ...\t";

# These two commands should fail with non-zero exit status.
#
# 1. An attempt to append binary to an existing ascii file.
system("cp bug14-ascii.hmm bug14.tmp");
$output = `$hmmconvert -Ab bug14-ascii.hmm bug14.tmp 2>/dev/null`;
if ($? == 0) { print "DETECTED! (1)\n"; unlink "bug14.tmp"; exit 1; }

# 2. An attempt to append ASCII (converted) to an existing binary file.
system("cp bug14-binary.hmm bug14.tmp");
$output = `$hmmconvert -A bug14-binary.hmm bug14.tmp 2>/dev/null`;
if ($? == 0) { print "DETECTED! (2)\n"; unlink "bug14.tmp"; exit 1; }

# These two commands should succeed with zero exit status.
#
# 3. An attempt to append ASCII (converted) to an existing ascii file
system("cp bug14-ascii.hmm bug14.tmp");
$output = `$hmmconvert -A bug14-binary.hmm bug14.tmp 2>/dev/null`;
if ($? != 0) { print "DETECTED! (3)\n"; unlink "bug14.tmp"; exit 1; }

# 4. An attempt to append binary (converted) to an existing binary file
system("cp bug14-binary.hmm bug14.tmp");
$output = `$hmmconvert -Ab bug14-ascii.hmm bug14.tmp 2>/dev/null`;
if ($? != 0) { print "DETECTED! (4)\n"; unlink "bug14.tmp"; exit 1; }

print "not detected.\n";
unlink "bug14.tmp";
exit 0;
