#! /usr/bin/perl

package h3;

sub ParseDomTbl {
    my ($domtblfile) = @_;
    my (@fields);

    $ndomtbl  = 0;
    @tname    = ();
    @tacc     = ();
    @qname    = ();
    @qlen     = ();
    @seqE     = ();
    @seqsc    = ();
    @seqbias  = ();
    @domidx   = ();
    @ndom     = ();
    @cE       = ();
    @iE       = ();
    @domsc    = ();
    @dombias  = ();
    @hmmi     = ();
    @hmmj     = ();
    @iali     = ();
    @jali     = ();
    @ienv     = ();
    @jenv     = ();
    @accuracy = ();
    @tdesc    = ();    

    if (! open(DOMFILE, $domtblfile)) { print "FAIL: couldn't open first domain table file"; exit 1 ; }
    while (<DOMFILE>)
    {
	if (/^\#/) { next; }
	chop;
	@fields = split(' ', $_, 23);

	$tname[$ndomtbl]   = $fields[0];
	$tacc[$ndomtbl]    = $fields[1];
	$tlen[$ndomtbl]    = $fields[2];
	$qname[$ndomtbl]   = $fields[3];
	$qacc[$ndomtbl]    = $fields[4];
	$qlen[$ndomtbl]    = $fields[5];
	$seqE[$ndomtbl]    = $fields[6];
	$seqsc[$ndomtbl]   = $fields[7];
	$seqbias[$ndomtbl] = $fields[8];
	$domidx[$ndomtbl]  = $fields[9];
	$ndom[$ndomtbl]    = $fields[10];
	$cE[$ndomtbl]      = $fields[11];
	$iE[$ndomtbl]      = $fields[12];
	$domsc[$ndomtbl]   = $fields[13];
	$dombias[$ndomtbl] = $fields[14];
	$hmmi[$ndomtbl]    = $fields[15];
	$hmmj[$ndomtbl]    = $fields[16];
	$iali[$ndomtbl]    = $fields[17];
	$jali[$ndomtbl]    = $fields[18];
	$ienv[$ndomtbl]    = $fields[19];
	$jenv[$ndomtbl]    = $fields[20];
	$accuracy[$ndomtbl]= $fields[21];
	$tdesc[$ndomtbl]   = $fields[22];
	$ndomtbl++;
    }
    close DOMFILE;
    1;
}

1;
