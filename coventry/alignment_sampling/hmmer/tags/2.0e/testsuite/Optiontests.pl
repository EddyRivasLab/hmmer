#! /usr/local/bin/perl

@tests = (
	  "hmmbuild -F Optiontests.hmm  Optiontests.slx",   # Make a protein HMM
	  "hmmbuild -F Optiontests.nhmm Optiontests.nslx",  # Make a DNA HMM
	  "hmmalign -h", 
	  "hmmalign Optiontests.hmm Optiontests.fa",
	  "hmmalign -m Optiontests.hmm Optiontests.fa",
	  "hmmalign -o tmp Optiontests.hmm Optiontests.fa",
	  "hmmalign -q Optiontests.hmm Optiontests.fa",
	  "hmmbuild -h",
	  "hmmbuild tmp.hmm Optiontests.slx",
	  "hmmbuild -F tmp.hmm Optiontests.slx",            # Need -F to force 
	  "hmmbuild -n foo -F tmp.hmm Optiontests.slx",
	  "hmmbuild -o foo -F tmp.hmm Optiontests.slx",
	  "hmmbuild -A tmp.hmm Optiontests.slx",
	  "hmmbuild -f -F tmp.hmm Optiontests.slx",
	  "hmmbuild -g -F tmp.hmm Optiontests.slx",
	  "hmmbuild -s -F tmp.hmm Optiontests.slx",
	  "hmmbuild --fast -F tmp.hmm Optiontests.slx",
	  "hmmbuild --hand -F tmp.hmm Optiontests.slx",
	  "hmmbuild --null  Optiontests.rnd -F tmp.hmm Optiontests.slx",
	  "hmmbuild --pam   Optiontests.pam -F tmp.hmm Optiontests.slx",
	  "hmmbuild --prior Optiontests.pri -F tmp.hmm Optiontests.slx",
	  "hmmbuild --wblosum -F tmp.hmm Optiontests.slx",
	  "hmmbuild --wgsc -F tmp.hmm Optiontests.slx",
	  "hmmbuild --wme -F tmp.hmm Optiontests.slx",
	  "hmmbuild --wvoronoi -F tmp.hmm Optiontests.slx",
	  "hmmbuild --wnone -F tmp.hmm Optiontests.slx",
	  "hmmbuild --noeff -F tmp.hmm Optiontests.slx",
	  "hmmbuild --amino -F tmp.hmm Optiontests.slx",
	  "hmmbuild --nucleic -F tmp.hmm Optiontests.nslx",
	  "hmmbuild --archpri 0.9 -F tmp.hmm Optiontests.slx",
	  "hmmbuild --binary -F tmp.hmm Optiontests.slx",
	  "hmmbuild --cfile tmp -F tmp.hmm Optiontests.slx",
	  "hmmbuild --gapmax 0.6 --fast -F tmp.hmm Optiontests.slx",
	  "hmmbuild --idlevel 0.5 -F tmp.hmm Optiontests.slx",
	  "hmmbuild --pamwgt 10 --pam Optiontests.pam -F tmp.hmm Optiontests.slx",
	  "hmmbuild --swentry 0.3 -F -s tmp.hmm Optiontests.slx",
	  "hmmbuild --swexit 0.3 -F -s tmp.hmm Optiontests.slx",
	  "hmmbuild --verbose -F tmp.hmm Optiontests.slx",
	  "hmmcalibrate -h",
	  "hmmcalibrate Optiontests.hmm",
	  "hmmcalibrate --fixed 15 Optiontests.hmm",
	  "hmmcalibrate --mean  25 Optiontests.hmm",
	  "hmmcalibrate --histfile tmp --fixed 15 Optiontests.hmm",
	  "hmmcalibrate --num 4500 --fixed 15 Optiontests.hmm",
	  "hmmcalibrate --sd 50 --mean  25   Optiontests.hmm",
	  "hmmcalibrate --seed 666 --fixed 15 Optiontests.hmm",
	  "hmmconvert -h",
	  "hmmconvert Optiontests.hmm tmp2.hmm",
	  "hmmconvert -F Optiontests.hmm tmp2.hmm",
	  "hmmconvert -a -F Optiontests.hmm tmp2.hmm",
	  "hmmconvert -A Optiontests.hmm tmp2.hmm",     # order sensitive. tmp2.hmm must be HMM
	  "hmmconvert -b -F Optiontests.hmm tmp2.hmm",
	  "hmmconvert -p -F Optiontests.hmm tmp2.hmm",
	  "hmmconvert -P -F Optiontests.hmm tmp2.hmm",
	  "hmmemit -h",
	  "hmmemit Optiontests.hmm",
	  "hmmemit -a Optiontests.hmm",
	  "hmmemit -n 6 Optiontests.hmm",
	  "hmmemit -o tmp Optiontests.hmm",
	  "hmmemit -q Optiontests.hmm",
	  "hmmemit --seed 666 Optiontests.hmm",
	  "hmmpfam -h",
	  "hmmpfam -n Optiontests.nhmm Optiontests.nfa",
	  "hmmpfam -A 0 Optiontests.hmm Optiontests.fa",
	  "hmmpfam -E 1 Optiontests.hmm Optiontests.fa",
	  "hmmpfam -T 1 Optiontests.hmm Optiontests.fa",
	  "hmmpfam -Z 10 Optiontests.hmm Optiontests.fa",
	  "hmmpfam --domE 1 Optiontests.hmm Optiontests.fa",
	  "hmmpfam --domT 1 Optiontests.hmm Optiontests.fa",
	  "hmmpfam --forward Optiontests.hmm Optiontests.fa",
	  "hmmpfam --null2 Optiontests.hmm Optiontests.fa",
	  "hmmpfam --xnu Optiontests.hmm Optiontests.fa",
	  "hmmsearch -h",
	  "hmmsearch -A 0 Optiontests.hmm Optiontests.fa",
	  "hmmsearch -E 1 Optiontests.hmm Optiontests.fa",
	  "hmmsearch -T 1 Optiontests.hmm Optiontests.fa",
	  "hmmsearch -Z 10 Optiontests.hmm Optiontests.fa",
	  "hmmsearch --domE 1 Optiontests.hmm Optiontests.fa",
	  "hmmsearch --domT 1 Optiontests.hmm Optiontests.fa",
	  "hmmsearch --forward Optiontests.hmm Optiontests.fa",
	  "hmmsearch --null2 Optiontests.hmm Optiontests.fa",
	  "hmmsearch --xnu Optiontests.hmm Optiontests.fa",
	  );


while ($testline = shift(@tests))
{
    $status = system("../src/$testline 2>&1 > tmp.out");
    if ($status > 0) {
	print "failure: $testline\n";
	$failed++;
    }
    $total++;
}

$passed = $total - $failed;
printf "Option tests: %d. Passed: %d. Failed: %d\n", $total, $passed, $failed;

unlink "tmp";
unlink "tmp.out";
unlink "tmp.hmm";
unlink "tmp2.hmm";
unlink "Optiontests.hmm";
unlink "Optiontests.nhmm";

