#!/usr/bin/perl -w
#pdb_parse.pl

use strict;
use Class::Struct;

# find directory where the script is installed

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  compare.pl [options] <pdbfile> <stofile> <rscapebin> <gnuplotdir> \n\n";
        print "options:\n";
        exit;
}

my $file1 = shift;
my $file2 = shift;

my $nhit1 = 0;
my @name1 = 0;
parse_hmm($file1, \$nhit1, \@name1);
my $nhit2 = 0;
my @name2 = 0;
parse_hmm($file2, \$nhit2, \@name2);

my $h1_not_h2 = 0;
for (my $h1 = 0; $h1 < $nhit1; $h1++) {
    my $found = 0;
    my $h2;
    for ($h2 = 0; $h2 < $nhit2; $h2++) {
	if ($name1[$h1] =~ /$name2[$h2]/) { last; }
    }
    if ($h2 == $nhit2) { $h1_not_h2 ++; printf "%d>$name1[$h1] not found in $file2\n", $h1_not_h2; }
}
print "$h1_not_h2/$nhit1 not found in $file2\n\n";

my $h2_not_h1 = 0;
for (my $h2 = 0; $h2 < $nhit2; $h2++) {
    my $found = 0;
    my $h1;
    for ($h1 = 0; $h1 < $nhit1; $h1++) {
	if ($name2[$h2] =~ /$name1[$h1]/) { last; }
    }
    if ($h1 == $nhit1) { $h2_not_h1 ++; printf "%d>$name2[$h2] not found in $file1\n", $h2_not_h1; }
}
print "$h2_not_h1/$nhit2 not found in $file1\n";

sub parse_hmm {
    my ($file, $ret_nhit, $name_ref) = @_;

    my $name = $file;
    if ($name =~ /^(\S+)\./)    { $name = $1; }
    if ($name =~ /^(\S+)\.max/) { $name = $1; }
    print "name $name\n";
    
    open (F, "$file") || die;
    while(<F>) {
	if (/^\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\d\s+($name\S+)\s+domain/) {
	    $name_ref->[$$ret_nhit] = $1;
	    $$ret_nhit ++;
	}
	elsif (/inclusion threshold/) { last; }
    }
    close (F);
    
    printf("$file nhit %d\n", $$ret_nhit);
    
}
