#! /usr/local/bin/perl

# rmanprocess.pl <rman LaTeX2e output>
# 
# Example:
#    rman -f LaTeX2e foo.man | rmanprocess.pl > foo.tex
#
# Converts a man page to a HMMER User's Guide section.
# Written to operate with PolyglotMan v3.0.5, by Thomas Phelps
# Obtain from ftp.cs.berkeley.edu:/ucb/people/phelps/tcltk/rman.tar.Z
#
# - removes document declarations
# - removes See Also and Author sections
# - converts sections to subsections
# - adds a section declaration for program name
# 
# 
# SRE, Mon May 25 11:06:58 1998


while (<>)
{
    if (/--/) { s/--/{-}{-}/g; }

    if (/^\\documentclass/) { 
	print "\\setlength{\\sresavei}{\\parindent}\n";
	print "\\setlength{\\sresaves}{\\parskip}\n";
	next;
    }
    if (/^\\begin\{document\}/) { next; }
    
    if (/^\\section\{See Also/) {
	print "\\setlength{\\parindent}{\\sresavei}\n";
	print "\\setlength{\\parskip}{\\sresaves}\n";
	print "\\newpage";
	last;
    }

    if (/\\begin\{itemize\}/ || /\\end\{itemize\}/) {
	s/itemize/wideitem/;
	print;
	next;
    }

    if (/^\\section\{Name/) {
	$line = <>;			# get begin itemize
	$line = <>;			# get item
	if ($line =~ /^\\item\s*\[(\S+)\s*-\s*(.+)\]/) {
	    print "\\section{\\texttt{$1} - $2}\n";
	}
	<>;			# get end itemize
	next;
    }

    if (/^\\section/) {
	s/section/subsection/;
	print;
	next;
    }

    print;

}
