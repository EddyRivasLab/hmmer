#! /bin/sh

etags    configure.ac
etags -a INSTALL
etags -a LICENSE

# Recursively add all .c, .h, .pl, *.tex, *.man
# `-L`:  follow symlinks, e.g. lib/easel
#
find -L . -name "*.c"   -print -or -name "*.h"  -print | xargs etags -a
find -L . -name "*.pl"  -print -or -name "*.pm" -print | xargs etags -a
find -L . -name "*.py"  -print                         | xargs etags -a
find -L . -name "*.sh"  -print                         | xargs etags -a
find -L . -name "*.md"  -print                         | xargs etags -a
find -L . -name "*.tex" -print                         | xargs etags -a
find -L . -name "*.man" -print                         | xargs etags -a
find -L . -name "*.in"  -print                         | xargs etags -a
find -L . -name "*.sqc" -print                         | xargs etags -a
find -L . -name "*README"    -print                    | xargs etags -a


