############################################################################
# demotic_fasta package   
#    Parses fasta or ssearch output, stores extracted information in convenient vars.
#    SRE, Wed Jun 25 13:41:41 2003
############################################################################
#  CVS $Id$
############################################################################

package demotic_fasta;

# parse(\*STDIN) would parse ssearch output
# coming in via stdin.
#
sub parse (*) {
    my $fh = shift;
    my $parsing_header  = 1;
    my $parsing_hitlist = 0;
    my $parsing_alilist = 0;
    my $target;
    my $alilinecount = 0;

    # Initialize everything... so we can call the parser
    # repeatedly, one call per ssearch output.
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
    %target_len     = ();	# Length of target sequence
    @hit_score      = ();	# Raw score (by rank)
    @hit_bitscore   = ();	# Bit score (by rank)
    @hit_Eval       = ();	# E-value (by rank)

				# The alignment output (in order it came in)
				# all indexed by # of alignment [0..nali-1]
    $nali           = 0;	# Number of alignments
    @ali_target     = ();	# Target sequence name
    @ali_score      = ();	# Smith/Waterman raw score of alignment
    @ali_bitscore   = ();	# bit score
    @ali_evalue     = ();	# E-value
    @ali_nident     = ();	# Number of identical residues
    @ali_alen       = ();	# Length of alignment (overlap)
    @ali_identity   = ();	# Percent identity
#    @ali_npos       = (); # Number of positives (similar positions)
#    @ali_positive   = (); # Percent of positions matched or similar
    @ali_qstart     = ();	# Start position on query
    @ali_qend       = ();	# End position on query
    @ali_tstart     = ();	# Start position on target
    @ali_tend       = ();	# End position on target
    @ali_qali       = (); # Aligned string from query
    @ali_tali       = (); # Aligned string from target (subject)
    
    # Now, loop over the lines of our input, and parse 'em.
    #
    while (<$fh>) {
	if ($parsing_header) {
	    if (/^The best scores are:/) { # start of hit list
		$parsing_header  = 0;
		$parsing_hitlist = 1;
		next;
	    } elsif (/^\s+\d+>>>\s*(\S*)\s*(.*)\s*-\s*(\d+) nt$/) { # allows blank query
		$queryname = $1;
		$querydesc = $2;
		$querylen  = $3;
		if ($queryname eq "") { 
		    $queryname = "unnamed_query";
		}
	    } elsif (/^\s+(\d+)\s+residues in\s+(\d+)\s+sequences\s*$/) {
		$db_nletters = $1;
		$db_nseq     = $2;
	    }
	} 
	elsif ($parsing_hitlist) {
	    if (/^\s*$/) {	# blank line marks end of hit list, start of alignments
		$parsing_hitlist = 0;
		$parsing_alilist = 1;
		next;
	    } elsif (/^(\S+)\s*(.*\S?)\s*\(\s*(\d+)\)\s+(\d+)\s+(\S+)\s+(\S+)\s*$/) {
		$hit_target[$nhits]    = $1;
		$target_desc{$1}       = $2;
	        $target_len{$1}        = $3;
		$hit_score[$nhits]     = $4;
		$hit_bitscore[$nhits]  = $5;
		$hit_Eval[$nhits]      = $6;
		$nhits++;
	    }
	}
	elsif ($parsing_alilist) {
	    if (/^>>(\S+)\s*(.*)\s+\((\d+) \S\S\)\s*$/) {  # the \S\S is either nt or aa
		$target = $1;
		$target_desc{$target} = $2;
		if ($3 != $target_len{$target}) { die "can't happen.", "1)", $3, "2)", $target_len{$target}; }
	    } 
	    elsif (/^ s-w opt:\s+(\d+)\s+Z-score:\s*(\S+)\s+bits:\s+(\S+)\s+E\(\):\s+(\S+)\s*$/) {  # SSEARCH
		$nali++;
		$ali_target[$nali-1]   = $target;
		$ali_score[$nali-1]    = $1;
		$ali_bitscore[$nali-1] = $3;
		$ali_evalue[$nali-1]   = $4;
	    } 
	    elsif (/^ initn:\s*\d+\s*init1:\s*\d+\s*opt:\s*(\d+)\s*Z-score:\s*(\S+)\s*bits:\s*(\S+)\s*E\(\):\s*(\S+)\s*$/) { # FASTA
		$nali++;
		$ali_target[$nali-1]   = $target;
		$ali_score[$nali-1]    = $1;
		$ali_bitscore[$nali-1] = $3;
		$ali_evalue[$nali-1]   = $4;
	    }		
	    elsif (/^Smith-Waterman score:\s+(\d+);\s+(\S+)% identity \(\S+% ungapped\) in (\d+) nt overlap \((\d+)-(\d+):(\d+)-(\d+)\)\s*/) {
		$ali_identity[$nali-1]   = $2;
		$ali_alen[$nali-1]       = $3;
		$ali_qstart[$nali-1]     = $4;
		$ali_qend[$nali-1]       = $5;
		$ali_tstart[$nali-1]     = $6;
		$ali_tend[$nali-1]       = $7;
#		$ali_nident[$nali-1]     = $1;
#		$ali_npos[$nali-1]       = $4;
#		$ali_positive[$nali-1]   = $5;
		$alilinecount            = 0;
	    } 
	    elsif (/^\S+\s+(\S+)\s*$/) { # only ali lines are right-flushed
		if ($alilinecount % 2 == 0) {
		    $ali_qali[$nali-1]  .= $1; 
		} else {
		    $ali_qali[$nali-1]  .= $1; 
		}
		$alilinecount++;
	    }
	}

    } # this closes the loop over lines in the input stream.
}

sub profmark_out {
    my $ofh = shift;
    my $i;

    for ($i = 0; $i < $nhits; $i++) {
	printf $ofh "%g\t%.1f\t%s\t%s\n", $hit_Eval[$i], $hit_bitscore[$i], $hit_target[$i], $queryname;
    }
}


sub exblxout {
    my $ofh     = shift;
    my $i;
    
    for ($i = 0; $i <= $nali-1; $i++) {
	printf $ofh "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%s\n",
	$ali_evalue[$i],
	$ali_identity[$i],
	$ali_tstart[$i],
	$ali_tend[$i],
	$ali_target[$i],
	$ali_qstart[$i],
	$ali_qend[$i],
	$queryname;
    }
}


sub gffout {
    my $ofh     = shift;
    my $source  = shift;
    my $feature = shift;
    my $i;
    my $strand;
    
    for ($i = 0; $i <= $nali-1; $i++) {
	if ($ali_qstart[$i] > $ali_qend[$i]) { $strand = "-"; }
	else { $strand = "+"; } 

	printf $ofh "%s\t%s\t%s\t%d\t%d\t%.1f\t%s\t.\tgene \"%s\"\n",
	$ali_target[$i],
	$source,
	$feature,
	$ali_tstart[$i],
	$ali_tend[$i],
	$ali_bitscore[$i],
	$strand,
	$queryname;
    }
}

1;
__END__
