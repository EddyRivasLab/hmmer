# Demotic parser for HMMER3 output
#
# SVN $Id$
# SRE, Tue Apr 22 13:13:51 2008

package demotic_hmmer;

# parse(\*STDIN) would parse BLAST output
# coming in via stdin.
#
sub parse (*) {
    my $fh = shift;
    my $parsing_header  = 1;
    my $parsing_seqs    = 0;
    my $parsing_domains = 0;
    my $parsing_alis    = 0;
    my $target;
    my $firstchunk;

    # Initialize everything... so we can call the parser
    # repeatedly, one call per hmmsearch output.
    #
    # This section also documents what's available in
    # variables within the package.
    # 
    # Arrays are indexed [0..nhits-1] or [0..nali-1]
    #
    $queryname      = "";	# Name of the query sequence
    $querydesc      = "";	# Description line for the query (or "")
    $querylen       = 0;	# Length of the query in residues
    $db             = "";	# Name of the target database file
    $db_nseq        = 0;	# Number of sequences in the target database
    $db_nletters    = "";	# Number of residues in the target database
                                # (stored as a string so we can have large #'s)

				# The top hit list (still in rank order)
    $nhits          = 0;	# Number of entries in the hit list
    @hit_target     = ();	# Target sequence name (by rank)
    %target_desc    = ();	# Target sequence description (by targ name)
    @hit_bitscore   = ();	# Raw score (by rank)
    @hit_Eval       = ();	# P-val or E-val (by rank)

				# The alignment output (in order it came in)
				# all indexed by # of alignment [0..nali-1]
    $nali           = 0;	# Number of alignments
    @ali_target     = ();	# Target sequence name
    @ali_score      = ();	# Raw score of alignment
    @ali_bitscore   = ();	# bit score
    @ali_evalue     = ();	# E-value
    @ali_pvalue     = ();	# P-value
    @ali_nident     = ();	# Number of identical residues
    @ali_alen       = ();	# Length of alignment
    @ali_identity   = ();	# Percent identity
    @ali_npos       = ();       # Number of positives (similar positions)
    @ali_positive   = ();       # Percent of positions matched or similar
    @ali_qstart     = ();	# Start position on query
    @ali_qend       = ();	# End position on query
    @ali_tstart     = ();	# Start position on target
    @ali_tend       = ();	# End position on target
    @ali_qali       = ();       # Aligned string from query
    @ali_tali       = ();       # Aligned string from target (subject)

    # Now, loop over the lines of our input, and parse 'em.
    #
    while (<$fh>) {
	if ($parsing_header) {
	    if (/^Scores for complete sequences /) { # seq section starts with this line
		$parsing_header  = 0;
		$parsing_seqs    = 1;
		<$fh>;		# This discards the next line:   --- full sequence ----  ---- single best dom ----
		<$fh>;		# and discard another line:      E-value   score   bias    score   bias 
		<$fh>;		# and discard another line:      ------- ------- ------  ------- ------ ----------
		next;
	    } elsif (/Query:\s*(\S+)\s*(.*)$/) {
		$queryname = $1;
		$querydesc = $2;
	    }
	}

	elsif ($parsing_seqs) {
	    if (/^\s*$/) { 	# seq section ends with a blank line
		$parsing_seqs    = 0;
		$parsing_domains = 1;
		next;
	    } elsif (/^\s*(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\d+)\s+(\S+)\s+(.+)$/) {
		#        evalue   score    bias evalue score  bias   exp      N     target  desc
                #           1       2       3                         4       5       6       7
		$hit_target[$nhits]    = $6;
		$target_desc{$1}       = $7;
		$hit_bitscore[$nhits]  = $2;
		$hit_Eval[$nhits]      = $1;
		$nhits++;
	    }
	}
    }
}


sub profmark_out {
    my $ofh = shift;
    my $i;

    for ($i = 0; $i < $nhits; $i++) {
	printf $ofh "%g\t%.1f\t%s\t%s\n", $hit_Eval[$i], $hit_bitscore[$i], $hit_target[$i], $queryname;
    }
}
1;
__END__
