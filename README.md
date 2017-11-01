## HMMER - biological sequence analysis using profile HMMs

[HMMER](http://hmmer.org) searches biological sequence databases for homologous sequences,
using either single sequences or multiple sequence alignments as
queries. HMMER implements a technology called "profile hidden Markov
models" (profile HMMs). HMMER is used by many protein family domain
databases and large-scale annotation pipelines, including
[Pfam](http://pfam.xfam.org) and other members of the
[InterPro Consortium](http://www.ebi.ac.uk/interpro/).

To obtain HMMER releases, please visit [hmmer.org](http://hmmer.org).

To participate in HMMER development, visit us at
[github](https://github.com/EddyRivasLab/hmmer).  HMMER development
depends on the Easel library, also at
[github](https://github.com/EddyRivasLab/easel).

We don't intend github to be a source of "official" HMMER release
code.  Our release tarballs bundle up some convenient stuff that you
have to be able to create for yourself if you're trying to do it from
our github repository.


### to download and build the current source code release:

Various tarballs of pre-compiled executables and source code for a
variety of systems are available at [hmmer.org](http://hmmer.org), but
it is also straightforward to download and build from source on most
systems:

```bash
   % wget http://eddylab.org/software/hmmer3/3.2/hmmer-3.2.tar.gz
   % tar zxf hmmer-3.2.tar.gz
   % cd hmmer-3.2
   % ./configure --prefix /your/install/path
   % make
   % make check
   % make install
``` 

Executable programs will be installed in `/your/install/path/bin`.

Files to read in the top-level source directory:

   * INSTALL - brief installation instructions.
   * Userguide.pdf - the HMMER User's Guide.
 
To get started after installation, see the Tutorial section in the
HMMER User's Guide (Userguide.pdf).


### to clone a copy of HMMER3 source from github:

You need to clone both the HMMER and Easel repositories, as follows:

```bash
   % git clone https://github.com/EddyRivasLab/hmmer
   % cd hmmer
   % git clone https://github.com/EddyRivasLab/easel
   % git checkout develop
   % (cd easel; git checkout develop)
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

To contribute to HMMER3 development, you want to be on the **develop**
branch.


### to report a problem:

Visit our
[issues tracking page at github](https://github.com/EddyRivasLab/hmmer/issues).

