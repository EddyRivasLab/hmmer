#! /usr/local/bin/perl

@tests = (
	  "hmmbuild --informat selex -F Optiontests.hmm  Optiontests.slx",   # Make a protein HMM
	  "hmmbuild --informat selex -F Optiontests.nhmm Optiontests.nslx",  # Make a DNA HMM
	  "hmmalign -h", 
	  "hmmalign Optiontests.hmm Optiontests.fa",
	  "hmmalign -m Optiontests.hmm Optiontests.fa",
	  "hmmalign -o tmp Optiontests.hmm Optiontests.fa",
	  "hmmalign -q Optiontests.hmm Optiontests.fa",
	  "hmmalign --withali Optiontests.slx Optiontests.hmm Optiontests.fa",
	  "hmmalign --mapali Optiontests.slx Optiontests.hmm Optiontests.fa",
	  "hmmbuild -h",
	  "hmmbuild --informat selex tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex -F tmp.hmm Optiontests.slx",            # Need -F to force 
	  "hmmbuild --informat selex -n foo -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex -o tmp -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex -A tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex -f -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex -g -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex -s -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --fast -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --hand -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --null ../tutorial/amino.null -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --pam   Optiontests.pam -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --prior ../tutorial/amino.pri -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --wblosum -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --wgsc -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --wme -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --wvoronoi -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --wnone -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --noeff -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --amino -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --nucleic -F tmp.hmm Optiontests.nslx",
	  "hmmbuild --informat selex --archpri 0.9 -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --binary -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --cfile tmp -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --gapmax 0.6 --fast -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --idlevel 0.5 -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --pamwgt 10 --pam Optiontests.pam -F tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --swentry 0.3 -F -s tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --swexit 0.3 -F -s tmp.hmm Optiontests.slx",
	  "hmmbuild --informat selex --verbose -F tmp.hmm Optiontests.slx",
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
	  "hmmindex -h",
	  "hmmindex Optiontests.hmm",
	  "hmmfetch -h",
	  "hmmfetch Optiontests.hmm Optiontests",
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


unlink "tmp.hmm";
while ($testline = shift(@tests))
{
    $status = system("../binaries/$testline 2>&1 > tmp.out");
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
unlink "Optiontests.hmm.ssi";

