### **HMMER: profile HMMs for biological sequence analysis**

HMMER searches biological sequence databases for homologous sequences,
using either single sequences or multiple sequence alignments as
queries. HMMER implements a technology called "profile hidden Markov
models" (profile HMMs). HMMER is the software engine underlying many
protein family domain databases and large-scale annotation pipelines,
including the Pfam and SMART databases.

Other files to read in the top-level source directory:

  * INSTALL -  brief installation instructions.
  * Userguide.pdf - the HMMER User's Guide.
  * LICENSE - copyright and license information.

To get started after installation, see the Tutorial section in the
HMMER User's Guide.

________________________________________________________________

To obtain HMMER releases, please visit:  http://hmmer.org

To participate in HMMER development, visit us at github: https://github.com/EddyRivasLab/hmmer

HMMER depends on the Easel library, also hosted at github:  https://github.com/EddyRivasLab/easel

To clone your own copy of HMMER source code for the first time:

```bash
   $ git clone https://github.com/EddyRivasLab/hmmer
   $ cd hmmer
   $ git clone https://github.com/EddyRivasLab/easel
   $ ln -s easel/aclocal.m4 aclocal.m4
   $ autoconf
   $ (cd easel; autoconf)
```

and to build:

```bash
   $ ./configure
   $ make
```

We follow a "git flow" workflow. The "master" branch is a stable
release branch. The "develop" branch is for integrating feature
branches. Most new features are developed on specific feature
branches.



