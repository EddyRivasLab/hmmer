############################################################################
# demotic_wublast package   
#    Parses blast output, stores extracted information in convenient vars.
#    Despite the name, this works on either NCBI-BLAST or WU-BLAST.
#    (SRE originated, 10/2000)
############################################################################
# CVS information
#  Revision $Revision: 1.8 $
#  Last author: $Author: eddy $
#  Last modification: $Date: 2002/06/21 14:54:57 $
############################################################################

package demotic_wublast;

# parse(\*STDIN) would parse BLAST output
# coming in via stdin.
#
sub parse (*) {
    my $fh = shift;
    my $parsing_header  = 1;
    my $parsing_hitlist = 0;
    my $parsing_alilist = 0;
    my $is_wublast      = 0;
    my $target;
    my $firstchunk;

    # Initialize everything... so we can call the parser
    # repeatedly, one call per BLAST output.
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
    @hit_Eval       = ();	# E-val (by rank)

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
    @ali_npos       = (); # Number of positives (similar positions)
    @ali_positive   = (); # Percent of positions matched or similar
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
	    if (/^Sequences producing /) { # wu and ncbi share this
		$parsing_header  = 0;
		$parsing_hitlist = 1;
		<$fh>;		# This discards the next (blank) line (ncbi, wu)
		next;
	    } elsif (/^Query=\s*(\S*)\s*(.*)\s*$/) { # allows blank query
		$queryname = $1;
		$querydesc = $2; chomp $querydesc;
		if ($queryname eq "") { 
		    $queryname = "unnamed_query";
		}
		while (1) {
		    $_ = <$fh>; # perl quirk: unless part of a while()
		                # loop, line must be explicitly assigned  
                                # to $_ or else it will be lost.
		    if (/^\s+\((\d+) letters/) {
			$querylen  = $1; 
			last;
		    } elsif (/^\s*( .+)\s*$/) { 
			$querydesc .= $1; chomp $querydesc;
		    }
		}
	    } elsif (/^Database:\s*(.+)\s*$/) {
		$db  = $1;
		$_ = <$fh>;
		if (/^\s+(\d+) sequences; (\S+) total letters/) {
		    $db_nseq     = $1;
		    $db_nletters = $2;
		}
	    } elsif (/^Copyright.+Washington University/) {
		$is_wublast = 1;
	    }
	} 
	elsif ($parsing_hitlist) {
	    if (/^\s*$/) { 
		$parsing_hitlist = 0;
		$parsing_alilist = 1;
		next;
	    } elsif (/^(\S+)\s+(.+)\s+(\d+)\s+(\S+)/) {
		$hit_target[$nhits]    = $1;
		$target_desc{$1}             = $2;
		$hit_bitscore[$nhits] = $3;
		if ($is_wublast) { $hit_Eval[$nhits] = -1.0 * log(1.0 - $4); } # conversion from P-value
		else             { $hit_Eval[$nhits] = $4; }

		$nhits++;
	    }
	}
	elsif ($parsing_alilist) {
	    if (/^>(\S+)\s*(.*)$/) {
		$target = $1;
		$target_desc{$target} = $2;

		$_ = <$fh>; 
		if (/^\s+Length = (\S+)/) { 
		    $target_len{$target} = $1;
		} 
	    } 
	    elsif (/^ Score =\s+(\d+) \((\S+) bits\), Expect = (\S+),/) { # WU
		$nali++;
		$ali_target[$nali-1]   = $target;
		$ali_score[$nali-1]    = $1;
		$ali_bitscore[$nali-1] = $2;
		$ali_evalue[$nali-1]   = $3;
	    } 
	    elsif (/^ Score =\s+(\S+) bits \((\S+)\), Expect = (\S+)/) { # NCBI
		$nali++;
		$ali_target[$nali-1]   = $target;
		$ali_bitscore[$nali-1] = $1;
		$ali_score[$nali-1]    = $2;
		$ali_evalue[$nali-1]   = $3;
	    }
	    elsif (/^ Identities = (\d+)\/(\d+) \((\d+)%\).+Positives = (\d+).+\((\d+)%/) { # NCBI or WU
		$ali_nident[$nali-1]     = $1;
		$ali_alen[$nali-1]       = $2;
		$ali_identity[$nali-1]   = $3;
		$ali_npos[$nali-1]       = $4;
		$ali_positive[$nali-1]   = $5;
		$firstchunk = 1;
	    } 
	    elsif (/^Query:\s+(\d+)\s+(\S+)\s+(\d+)\s*$/) {
		if ($firstchunk) { $ali_qstart[$nali-1] = $1; }
		$ali_qali[$nali-1]  .= $2;
		$ali_qend[$nali-1]   = $3;
	    } 
	    elsif (/^Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)\s*$/) {
		if ($firstchunk) { $ali_tstart[$nali-1] = $1; }
		$ali_tali[$nali-1]  .= $2;
		$ali_tend[$nali-1]   = $3;
		$firstchunk = 0;
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
