# @LICENSE@

# hmmer.pm
# Perl routines for parsing HMMER output
#
# SRE, Wed Oct 28 11:27:17 1998
# CVS $Id$ 

package hmmer;

#------------ ParseHMMER ------------
#
# Parse hmmsearch or hmmpfam output into
# arrays that we can use in perl scripts.
#
# Nothing fancy, we just put the info into
# module-specific variables, expecting that we'll just
# be writing small scripts with this module.
#
# Illustrative example:
#    use hmmer;
#    $output = `hmmsearch foo.hmm swissprot36`;
#    &hmmer::ParseHMMER($output);
#    printf "The top scoring sequence is %s\n", $hmmer:targname[0];
#    printf "The total number of domains hit is %d\n", $hmmer::ndom;
#
# Data made available:
#    $query        - name of query HMM or sequence, e.g. "rrm"
#    $querydesc    - description of query, e.g. "RNA recognition motif"
#
#    $ntarget      - total number of targets hit in per-seq output
#    @targname     - array of target names (e.g. RU1A_HUMAN)
#    %targdesc     - target descriptions (indexed by target name)
#    %seqscore     - per-seq score (indexed by target name)
#    %seqeval      - per-seq E-value (indexed by target name)
#    %seqndom      - number of domains hit (indexed by target name)
#
#    $ndom         - total number of hits in domain output
#    @domname      - target names hit e.g. $domname[0] = "RU1A_HUMAN"
#    @domnum       - e.g. "1/2"
#    @domsqfrom    - sequence from coords (start positions)
#    @domsqto      - sequence to coords (end positions)
#    @domsqbounds  - e.g. "[]" or ".." for seq
#    @domhmmfrom   - array of hmm-from coords
#    @domhmmto     - array of hmm-to coords
#    @domhmmbounds - e.g. "[]" or ".." for HMM
#    @domscore     - domain scores
#    @domevalue    - domain E-values
#
#    $aligndata    - the raw alignment text (currently not parsed further)
#
sub ParseHMMER {
    my($output) = @_;
    my($indom, $inseq, $inali);
    my(@lines, $line);

    $query       = "";
    $querydesc   = "";

    $ntarget     = 0;
    @targname    = ();
    %targdesc    = ();
    %seqscore    = ();
    %seqeval     = ();
    %seqndom     = ();

    $ndom        = 0;
    @domname     = ();
    @domnum      = ();
    @domsqfrom   = ();
    @domsqto     = ();
    @domsqbounds = ();
    @domhmmfrom  = ();
    @domhmmto    = ();
    @domhmmbounds= ();
    @domscore    = ();
    @domevalue   = ();
    $aligndata   = "";

    @lines = split(/^/, $output);
    $ndom=0;
    $ntarget=0;
    foreach $line (@lines) 
    {
	if ($line =~ /^Query:\s+(\S+)\s+(.+)$/)        {$query = $1; $querydesc = $2;}
	if ($line =~ /^Scores for/)                    {$indom = 0;  $inseq = 1;}
	if ($line =~ /^Parsed for domains/)            {$indom = 1;  $inseq = 0;}
	if ($line =~ /^Histogram of all scores/)       {$indom = 0;  $inseq = 0; $inali = 0; }
	if ($line =~ /^Alignments of top-scoring/)     {$inali = 1;  $indom = 0; $inseq = 0; }

	if ($inseq && $line =~ /^\s*(\S+)\s+(.+?)\s+(\S+)\s+(\S+)\s+(\d+)\s*$/)
	{
	    $targname[$ntarget]= $1;
	    $targdesc{$1}      = $2;
	    $seqscore{$1}      = $3;
	    $seqeval{$1}       = $4;
	    $seqndom{$1}       = $5;
	    $ntarget++;
	}
	if ($indom && $line =~ /^\s*\S+\s+\d+\/\d+/)
	{
	    ($domname[$ndom],
	     $domnum[$ndom],
	     $domsqfrom[$ndom],
	     $domsqto[$ndom],
	     $domsqbounds[$ndom],
	     $domhmmfrom[$ndom],
	     $domhmmto[$ndom],
	     $domhmmbounds[$ndom],
	     $domscore[$ndom],
	     $domevalue[$ndom]) = split ' ', $line;

	     $ndom++;
	}

	if ($inali) { $aligndata .= $line; }
    }
    1;
}



1;
