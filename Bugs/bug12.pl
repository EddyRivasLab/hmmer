#! /usr/local/bin/perl

# Setting a --domE or --domT threshold that removes the second
# domain from:
#    rrm        1/3     151   235 ..     1    72 []    82.9  6.6e-21
#    rrm        2/3     250   322 ..     1    72 []    77.0  3.9e-19
#    rrm        3/3     404   475 ..     1    72 []   106.3    6e-28
# should not remove domain 3 too.
# SRE, Fri Dec 15 10:04:41 2000
#

$hmmpfam = "../binaries/hmmpfam";
print "Bug 12 ...\t";

$output = `$hmmpfam bug12-minipfam bug12-ELAV_DROME`;
if (! &hmmpfam_output_has_domain($output, "rrm", 151, 235)) { print("OOPS. er, something else is wrong (1).\n"); exit 1;}
if (! &hmmpfam_output_has_domain($output, "rrm", 250, 322)) { print("OOPS. er, something else is wrong (2).\n"); exit 1;}
if (! &hmmpfam_output_has_domain($output, "rrm", 404, 475)) { print("OOPS. er, something else is wrong (3).\n"); exit 1;}

$output = `$hmmpfam --domE 9e-24 bug12-minipfam bug12-ELAV_DROME`;
if (! &hmmpfam_output_has_domain($output, "rrm", 151, 235)) { print("OOPS. er, something else is wrong (4).\n"); exit 1;}
if (! &hmmpfam_output_has_domain($output, "rrm", 404, 475)) { print("DETECTED! (1)\n"); exit 1;}

$output = `$hmmpfam --domT 80 bug12-minipfam bug12-ELAV_DROME`;
if (! &hmmpfam_output_has_domain($output, "rrm", 151, 235)) { print("OOPS. er, something else is wrong. (5)\n"); exit 1;}
if (! &hmmpfam_output_has_domain($output, "rrm", 404, 475)) { print("DETECTED! (2)\n"); exit 1;}

print "not detected.\n";
exit 0;

sub 
hmmpfam_output_has_domain {
    my $output       = shift;
    my $domain_name  = shift;
    my $domain_start = shift;
    my $domain_end   = shift;

    foreach $line (split(/^/, $output)) {
	if ($line =~ /^Parsed for domains/) { $in_domains = 1; }
	if ($in_domains &&
	    $line =~ /^${domain_name}\s+\S+\s+${domain_start}\s+${domain_end}\s/)
	{ return 1; }
	if ($line =~ /^Alignments of top-scoring/) { break; }
    }
    return 0;
}
