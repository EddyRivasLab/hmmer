#! /usr/bin/perl


# Example of a HMMER benchmark:
#   ./pmark-master.pl h3-results    100 pmark ./pmark-h3
#   ./pmark-master.pl h2-results-ls 100 pmark ./pmark-h2-ls
#   ./pmark-master.pl h2-results-fs 100 pmark ./pmark-h2-fs
#

$resultdir    = shift;
$ncpu         = shift;
$benchmark    = shift;
$pmark        = shift;

$tbl          = "$benchmark.tbl";
$msafile      = "$benchmark.msa";
$fafile       = "$benchmark.fa";

if (-e $resultdir) { die("$resultdir exists");}
system("mkdir $resultdir");

# Suck in the master table
open(BENCHMARK_TBL, $tbl) || die;
$n = 0;
while (<BENCHMARK_TBL>) 
{
    ($msaname[$n], $pid, $L, $nseq) = split;
    $alen{$msaname[$n]} = $L;
    $n++;
}
close BENCHMARK_TBL;

# Sort it by alen - this helps load balance.
sub by_alen { $alen{$b} <=> $alen{$a} }
@sorted_msaname = sort by_alen @msaname;

# Create <ncpu> subtables.
for ($i = 0; $i < $n; $i++)
{
    $subtbl[$i % $ncpu] .= $sorted_msaname[$i];
    $subtbl[$i % $ncpu] .= "\n";
}

# Output the <ncpu> subtables
for ($i = 0; $i < $ncpu; $i++)
{
    open(SUBTBL, ">$resultdir/tbl.$i") || die ("Failed to create $resultdir/tbl.$i");
    print SUBTBL $subtbl[$i];
    close SUBTBL;
}

# Submit all the individual profmark jobs
for ($i = 0; $i < $ncpu; $i++)
{
   system("qsub -V -cwd -b y -N tbl$i -j y -o $resultdir/tbl$i.sge '$pmark $resultdir $resultdir/tbl.$i $msafile $fafile $resultdir/tbl$i.out'");
}

