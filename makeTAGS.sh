#! /bin/sh

# Default: index both hmmer and easel. To do easel, need `find -L` to
# follow symlinks.
#
# `-x`: only index hmmer. Higher-level ~/src/makeTAGS.sh calls us with
# that option, for indexing entire stack of lab code (easel, hmmer3,
# hmmer4, infernal) without redundancy.
#

while [[ "$#" -gt 0 ]]; do case $1 in
  -x) optx=1;;
  *) echo "Unknown option: $1"; exit 1;;
esac; shift; done

if [ $optx ]; then
    opt=""
    excl=" -path ./easel -prune -or"
    out="-o TAGS.part"
else
    opt="-L"
    excl=""
    out=""
fi    


etags    $out configure.ac
etags -a $out INSTALL
etags -a $out LICENSE

find $opt . $excl -name "*.c"   -print -or -name "*.h"  -print | xargs etags -a $out
find $opt . $excl -name "*.pl"  -print -or -name "*.pm" -print | xargs etags -a $out
find $opt . $excl -name "*.py"  -print                         | xargs etags -a $out
find $opt . $excl -name "*.sh"  -print                         | xargs etags -a $out
find $opt . $excl -name "*.md"  -print                         | xargs etags -a $out
find $opt . $excl -name "*.tex" -print                         | xargs etags -a $out
find $opt . $excl -name "*.man" -print                         | xargs etags -a $out
find $opt . $excl -name "*.in"  -print                         | xargs etags -a $out
find $opt . $excl -name "*.sqc" -print                         | xargs etags -a $out
find $opt . $excl -name "*README" -print                       | xargs etags -a $out


