#! /usr/local/bin/perl

# licenseadd.pl <licensefile> <source file>
#
# Replaces occurrences of @LICENSE@ in the source file
# with the text found in <licensefile>. Since it expands
# to a multiline license within a comment, it needs to
# be a little careful about the context of the comment
# A (C code) line like
#     * @LICENSE@ 
# is replaced with 
#     * An example license.
#     * Copyright (C) ...
#
# A (shell or Perl script) line like
#     # @LICENSE@
# is replaced with 
#     # An example license
#     # Copyright (C) ...
# 
# An HTML section like
#   <!--  
#     -- @LICENSE@ 
#     -->
# is replaced with
#     <!-- 
#       -- An example license 
#       -- Copyright (C) ...
#       -->
#
$licensefile = shift;
$sourcefile  = shift;

if (! -e $sourcefile) { die "no such file $sourcefile"; }
($dev,$ino,$mode) = stat($sourcefile);

open(LICENSE,$licensefile) || die;
$nlines = 0;
while (<LICENSE>)
{
    chomp;
    $licenseline[$nlines] = $_;
    $nlines++;
}
close(LICENSE);

open(TMPFILE,">/tmp/tmplicense") || die "Fatal: can't open /tmp/tmplicense : $!\n";
open(SOURCE,$sourcefile) || die;
while (<SOURCE>) 
{
    if (/^(.*)\@LICENSE\@(.*)$/) 
    {
	$start = $1;
	$end   = $2;
	foreach $line (@licenseline) 
	{
	    print TMPFILE "$start$line$end\n";
	}
    } else { print TMPFILE $_;}
}
close SOURCE;
close TMPFILE;

# Replace the original file with the new one, and restore the original
# file's mode.
#
unlink $sourcefile;
system("mv /tmp/tmplicense $sourcefile");
chmod $mode, $sourcefile;
