## Port of HMMER to Windows using MinGW

The following recipie describes how to build native executables for windows using 
[MSYS2](https://www.msys2.org/) and [MInGW](https://www.mingw-w64.org/). At the time
 of writing the compiler was `gcc version 11.2.0 (Rev1, Built by MSYS2 project)`

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
mingw32-make.exe dev
mingw32-make.exe install
```

Configure HMMR

If `configure` is not present it will need to be generated on a LINUX system using `autoconf`.

Run configure (here installing to c:\local and looking for mman header files and library in c:\local):

```console
CPPFLAGS="-Ic:/local/include" LDFLAGS=" -Lc:/local/lib" CFLAGS="-g -O0" ./configure --prefix=c:/local
```
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

Will compile but not generally run (at least not in a any portable fashion).
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
the test scripts frequently need changes even when the commands function correcty).


