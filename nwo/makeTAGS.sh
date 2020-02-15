#! /bin/sh

# This is a temporary makeTAGS.sh that only indexes nwo.
# It doesn't index easel, nor ../src.

etags Makefile.in

find . -name "*.c"   -print -or -name "*.h"  -print | xargs etags -a
find . -name "*.pl"  -print -or -name "*.pm" -print | xargs etags -a
find . -name "*.py"  -print                         | xargs etags -a
find . -name "*.sh"  -print                         | xargs etags -a
find . -name "*.md"  -print                         | xargs etags -a
find . -name "*.tex" -print                         | xargs etags -a
find . -name "*.man" -print                         | xargs etags -a
find . -name "*.in"  -print                         | xargs etags -a
find . -name "*.sqc" -print                         | xargs etags -a
find . -name "*README"    -print                    | xargs etags -a


