## Port of HMMER to Windows using MinGW

The following recipie describes how to build native executables for windows using 
[MSYS2](https://www.msys2.org/) and [MInGW](https://www.mingw-w64.org/).

Prior to building install MSYS2, I used [chocolately](https://chocolatey.org/) 
(`choco.exe install msys2`), but you can also use an
[installer](https://repo.msys2.org/distrib/x86_64/)


In the MSYS2 shell install these packages (picking all groups)

```console
pacman -S base-devel gcc vim cmake 
pacman -S mingw-w64-x86_64-toolchain
pacman -S diffutils
```

There is no mman library for windows, but a wrapper is available for MinGW.  
Download it from [github](https://github.com/alitrack/mman-win32) and build 
in MSYS2 (here installing it to c:\local).

```console
./configure --prefix=c:/local
mingw32-make.exe dev
mingw32-make.exe install
```

Configure HMMR
```console
CPPFLAGS="-Ic:/local/include" LDFLAGS=" -Lc:/local/lib" CFLAGS="-g -O0" ./configure --prefix=c:/local
```

If `configure` is not present it will need to be generated on a LINUX system using `autoconf`.

Check that configure has picked up the mman library.

```console
HMMER configuration:
   compiler:             gcc -g -O0   -pthread
   host:                 x86_64-w64-mingw32
   linker:               -L/c/local/lib
   libraries:            -lmman   -lpthread
   DP implementation:    sse
```

Make, test and install.  

```console
mingw32-make.exe dev           # build everything- there should be no warnings
rm -rf src/impl                
cp -r src/impl_sse src/impl    # No symlinks on windows, so this needs to be copied
mingw32-make.exe check         # Run the tests
mings32-make.exe install       # install if you're happy with test results
```

### Porting notes

Windows (and MinGW) do not support signals and sockets.  Programs that use sockets will not work:

Use of shell commands is not supported.  For example, code like 

```C
popen("cat esltmpfile 2>/dev/null");
system("gzip -c eslfile 2>/dev/null > eslfile.gz");
```

Will compile but not run (at least not in a any portable fashion).
I've disabled unit tests that use this type of code.

The temporary files created by `esl_tmpfile` will persist on the filesystem in TEMP.
In Windows it is not possible to delete a file that is in use.

### Failed tests

All C unit tests should pass.

Many of the perl test script use shell

i1-degen-residues.pl fails with `FAIL: reformat changed .dna test` L46 failed system command.

esl-afetch.itest.pl fails with `FAIL: esl-afetch fetched incorrectly at ./esl-afetch.itest.pl line 61.` L61 pattern match fails because of windows CRLF




### Failed tests

## HMMER - biological sequence analysis using profile HMMs

[![](https://travis-ci.org/EddyRivasLab/hmmer.svg?branch=develop)](https://travis-ci.org/EddyRivasLab/hmmer)
![](http://img.shields.io/badge/license-BSD-brightgreen.svg)

[HMMER](http://hmmer.org) searches biological sequence databases for
homologous sequences, using either single sequences or multiple
sequence alignments as queries. HMMER implements a technology called
"profile hidden Markov models" (profile HMMs). HMMER is used by many
protein family domain databases and large-scale annotation pipelines,
including [Pfam](http://pfam.xfam.org) and other members of the
[InterPro Consortium](http://www.ebi.ac.uk/interpro/).

To obtain HMMER releases, please visit [hmmer.org](http://hmmer.org).

To participate in HMMER development, visit us at
[github](https://github.com/EddyRivasLab/hmmer).  HMMER development
depends on the Easel library, also at
[github](https://github.com/EddyRivasLab/easel).


### to download and build the current source code release:

```
   % wget http://eddylab.org/software/hmmer/hmmer.tar.gz
   % tar zxf hmmer.tar.gz
   % cd hmmer-3.3.2
   % ./configure --prefix /your/install/path
   % make
   % make check                 # optional: run automated tests
   % make install               # optional: install HMMER programs, man pages
   % (cd easel; make install)   # optional: install Easel tools
``` 

Executable programs will be installed in `/your/install/path/bin`. If
you leave this optional `./configure` argument off, the default prefix
is `/usr/local`.

Files to read in the source directory:

   * INSTALL - brief installation instructions.
   * Userguide.pdf - the HMMER User's Guide.
 
To get started after installation, see the Tutorial section in the
HMMER User's Guide (Userguide.pdf).



### to clone a copy of HMMER3 source from github:

The tarball way, above, is a better way to install HMMER (it includes
a precompiled Userguide.pdf, for example), but you can also clone our
github repo. You need to clone both the HMMER and Easel repositories,
as follows:

```
   % git clone https://github.com/EddyRivasLab/hmmer
   % cd hmmer
   % git clone https://github.com/EddyRivasLab/easel
   % autoconf
```

and to build:

```bash
   % ./configure
   % make
```

Our [git workflow](https://github.com/EddyRivasLab/hmmer/wiki/Git-workflow)
includes three main branches:

 * **master** is the stable branch for HMMER3 releases (including when
   H3 is released as a library inside Infernal)
 * **develop** is the HMMER3 development branch
 * **h4-develop** is the HMMER4 development branch.

To build the most recent official release, leave both HMMER and Easel
on their default **master** branch.  To contribute to HMMER3
development, you want to be on the **develop** branches. If you want
to send us a pull request on GitHub, please base your changes on our
**develop** branches.


### to report a problem:

Visit our
[issues tracking page at github](https://github.com/EddyRivasLab/hmmer/issues).

