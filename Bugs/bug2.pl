#! /usr/local/bin/perl

# The sum of domain hits should equal the total sequence score.
# Else, it's possible to have a sequence that actually scores
# *less* than some of its individual domains... such as
# the example in the bug 2 report, Pfam's spectrin HMM vs.
# S30431.
#
# xref STL3 notebook p. 127.
# SRE, Tue Dec 19 14:47:30 2000

use hmmer;

$usage = "Usage: ./bug2.pl [BINARYPATH]\n";

$binpath = "";
if    ($#ARGV > 0)  { print $usage; exit 1; }
elsif ($#ARGV == 0) { $binpath = shift;     }

$hmmpfam   = "$binpath/hmmpfam";
$hmmsearch = "$binpath/hmmsearch";
print "Bug 2 ...\t";

$output = `$hmmpfam -E100 bug2-spectrin.hmm bug2-S30431.fa`;
&hmmer::ParseHMMER($output);
$seqsc          = $hmmer::seqscore{"spectrin"};
$sum_of_domains = 0.;
foreach $sc (@hmmer::domscore) { $sum_of_domains += $sc; }
if (abs($sum_of_domains-$seqsc) >=0.1) { print "DETECTED! (1)\n";  exit 1;}

$output = `$hmmsearch -E100 bug2-spectrin.hmm bug2-S30431.fa`;
&hmmer::ParseHMMER($output);
$seqsc          = $hmmer::seqscore{"S30431"};
$sum_of_domains = 0.;
foreach $sc (@hmmer::domscore) { $sum_of_domains += $sc; }
if (abs($sum_of_domains-$seqsc) >= 0.1) { print "DETECTED! (2)\n"; exit 1;}

print "not detected.\n";
exit 0;
