#! /usr/local/bin/perl

# checkoptions.pl <program>
#
# - Runs binaries/<program> -h to extract a list of options.
# - Checks testsuite/Optiontests.pl to be sure each one is run at least
#   once in an option test.
#     Special case: don't check --pvm or --cpu, because we don't know
#     whether support was compiled in.
# - Checks documentation/man/<program>.man to check that each option is documented.
#
$program = shift;

die("no program binaries/$program") unless -e "binaries/$program";

# Get options from binaries/<program>
#
$output = `binaries/$program -h`;
@lines  = split(/^/, $output);
foreach $line (@lines) {
    if ($line =~ /^\s+(\S+)\s*(\S*)\s*:/) { $progoption{$1} = 1; $progarg{$1} = $2;}
}

# Get tested options from testsuite/Optiontests.pl
#
open(FILE, "testsuite/Optiontests.pl") || die "oops";
while (<FILE>) {
    if (/^\s+\"(\S+) (\S+?)[\" ]/) {
	$prg = $1;
	$opt = $2;

	if ($prg eq $program && $opt =~ /^-/) {
	    $testedoption{$opt} = 1;
	}
    }
}
close FILE;
	

# Get documented options from Man/<program>.man
#
die "documentation/man/$program.man doesn't exist!\n" unless -e "documentation/man/$program.man";
open(MAN, "documentation/man/$program.man") || die "oops";
while (<MAN>) {
    if (/^\.SH OPTIONS/ || /^\.SH EXPERT OPTIONS/) {$inoptions = 1;}
    elsif (/^\.SH/)      {$inoptions = 0;}

    if ($inoptions) {
	if ($next_is_option) {
	    if (/^\.B (\S+)/) {
		$docoption{$1} = 1;
		$docarg{$1}    = "";
	    } elsif (/^\.BI (\S+) \" (\S+)\"/) {
		$docoption{$1} = 1;
		$docarg{$1}    = $2;
	    } elsif (/^\.BI (\S+)/) {
		$docoption{$1} = 1;
	    }
	} 

	if (/^\.TP/)  { $next_is_option = 1; } else {$next_is_option = 0;} 
    }
}
close MAN;

# Check that all prog options are tested, and no others
#
foreach $opt (keys(%progoption)) {
    next if ($opt == "--cpu" || $opt == "--pvm"); 
    if (! $testedoption{$opt}) { print("$program $opt is not tested in Shiva\n"); }
}
foreach $opt (keys(%testedoption)) {
    if (! $progoption{$opt}) { print("$program $opt is tested in Shiva, but is not in $program -h\n"); }
}

# Check that all prog options are documented, and no others
#
foreach $opt (keys(%progoption)) {
    if (! $docoption{$opt}) { print("$program $opt is not documented in Man\n"); }
    else {
	if ($docarg{$opt} ne $progarg{$opt}) {
	    printf("%s %s documented for arg of %s but -h shows arg of %s\n",
		   $program, $opt, $docarg{$opt}, $progarg{$opt});
	}
    }
}
