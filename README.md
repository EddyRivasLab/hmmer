# Port of HMMER to Windows using MinGW

This project is a port of HMMER source code to the [MInGW](https://www.mingw-w64.org) compiler.  

With the exception of `hmmpgmd` and `hmmpgmd_shard` the Windows native executables pass all tests.

If prefer not to do your own build the resulting executables are in the [binaries](/binaries) directory.

The source code in this project was created from [hmmer](https://github.com/jones-gareth/hmmer/tree/mingw-build)
and [easel](https://github.com/jones-gareth/easel/tree/mingw-build) forks of the develop branches 
of [HMMER](https://github.com/EddyRivasLab/hmmer/tree/develop).

## Build Recipe

The following recipe describes how to build native executables for windows using 
[MSYS2](https://www.msys2.org/) and [MInGW](https://www.mingw-w64.org/). At the time
 of writing the compiler was `gcc version 11.2.0 (Rev1, Built by MSYS2 project)`

**Prerequisites**

Prior to building install MSYS2, For example use [chocolately](https://chocolatey.org/) 
(`choco.exe install msys2`), or an
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
mingw32-make.exe 
mingw32-make.exe install
```

**Download code**

```console
git clone https://github.com/jones-gareth/hmmer
cd hmmer
git clone https://github.com/jones-gareth/easel
```

**Configure HMMER**

If `configure` is not present it will need to be generated on a LINUX system using `autoconf`.

Run configure (here installing to c:\local and looking for mman header files and library in c:\local):

```console
CPPFLAGS="-Ic:/local/include" LDFLAGS=" -static -Lc:/local/lib" CFLAGS="-g -O0" ./configure --prefix=c:/local
```
Check that configure has picked up the mman library.

```console
HMMER configuration:
   compiler:             gcc -g -O0   -pthread
   host:                 x86_64-w64-mingw32
   linker:                -static -Lc:/local/lib
   libraries:            -lmman   -lpthread
   DP implementation:    sse
```

**Make, test and install**  

```console
mingw32-make.exe clean
mingw32-make.exe dev           # build everything- there should be no warnings
rm -rf src/impl                
cp -r src/impl_sse src/impl    # No symlinks on windows, so this needs to be copied
mingw32-make.exe check         # Run the tests. 2 will fail
mingw32-make.exe install       # install if you're happy with test results
```

### Porting notes

The presence of carriage returns (`\r`) in `config.sub` or `config.guess` will cause `autoconf` to fail.
Git operations may result in the insertion of `\r` depending on the `core.autocrlf` setting.

Windows (and MinGW) do not support POSIX signals or sockets.  Lack of socket support
mean that the  client server daemon 
commands `hmmpgmd` and `hmmpgmd_shard` do not function (the executables are compiled, but they
won't work.)

Use of shell commands in C system calls is not properly supported.  For example, code like 

```C
popen("cat esltmpfile 2>/dev/null");
system("gzip -c eslfile 2>/dev/null > eslfile.gz");
```

Will compile but not generally run (at least not in any portable fashion).
I  disabled unit tests that use this type of code.

1. Case 5 in esl_buffer_utest::main
2. esl_buffer_utest::utestOpenPipe
3. The gzip output case in esl_buffer_utest::utestSetOffset
4. esl_buffer_utest::utest_halfnewline (requires signals)

The temporary files created by `esl_tmpfile` will persist on the filesystem in TEMP.
In Windows it is not possible to delete a file that is in use.

### Failed tests

All tests except `hmmpgmd_ga` and `hmmpgmd_shard_ga` should pass.
Since the `hmmpgmd` and `hmmpgmd_shard` commands are broken that is expected.

```console
    exercise  297 [           hmmpgmd_ga] ...     FAILED [command failed]
    exercise  299 [     hmmpgmd_shard_ga] ...     FAILED [command failed]

2 of 307 exercises at level <= 2 FAILED.
```

The tests use the perl and python executables installed in the MSYS2 shell.

Test code at a level > 2 has not been ported and may not run correctly (because of windows text mode,
the test scripts frequently need changes even when the commands function correctly).

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

