#! /usr/bin/perl

# Nightly builds for HMMER3
# 
# Usage:     hmmer_autobuilds.pl  <srcdir>
# Example:   hmmer_autobuilds.pl  ~/nightlies/hmmer/trunk > /tmp/hmmer_autobuilds.log

@buildconfigs = (
    { name => "intel-macosx-gcc",            host => ".",            use_qsub => 0  },
    { name => "intel-macosx-gcc-debug",      host => ".",            use_qsub => 0  },
    { name => "intel-macosx-gcc-mpi",        host => ".",            use_qsub => 0  },
    { name => "intel-linux-gcc",             host => "login-eddy",   use_qsub => 0  },
    { name => "intel-linux-icc-intel64-mpi", host => "login-eddy",   use_qsub => 1  }, # Actually compiles/tests on a cluster node.
    { name => "intel-linux-icc-ia32",        host => "login-eddy",   use_qsub => 0  },
    { name => "intel-linux-gcc-ubuntu32",    host => "cf-ubuntu32",  use_qsub => 0  },
    { name => "intel-freebsd-gcc",           host => "cf-freebsd",   use_qsub => 0  },
    );

$autoconf = "/opt/local/bin/autoconf";                          # on the build master: wol
$qsub     = ". /sge/8.0.1p4/default/common/settings.sh; qsub";  # on the SGE master:   login-eddy

if ($#ARGV+1 != 1) { die "FAIL: incorrect number of command line arguments"; }
$srcdir = shift;
if (! -d $srcdir) { die "FAIL: source working directory $srcdir not found"; }

# First we update in the source working directory.
#
chdir $srcdir || die "FAIL: couldn't cd to $srcdir"; 
system("svn update                > autobuilds.log 2>&1");   if ($?) { die "FAIL: svn update"; }
system("$autoconf                 >> autobuilds.log 2>&1 "); if ($?) { die "FAIL: H3 $autoconf"; }
system("(cd lib/easel; $autoconf) >> autobuilds.log 2>&1");  if ($?) { die "FAIL: esl $autoconf"; }

# Then we try to build on everything
#
foreach $build (@buildconfigs)
{
    $script = $build->{name};
    $cmd = "$srcdir/autobuild/autobuild.pl $srcdir $script";
    if ($build->{use_qsub})    { $cmd = "($qsub -j y -o $srcdir/autobuilds.qsub -N autobuild -sync y -V -cwd -b y '$cmd' > /dev/null 2>&1; cat $srcdir/autobuilds.qsub)"; } # See [qsub].
    if ($build->{host} ne ".") { $cmd = "ssh $build->{host} \"$cmd\""; }

    printf ("%-30s %-20s ", $build->{name}, $build->{host});
    $output = `$cmd`;
    print $output;

    # cleanup
    if ($build->{use_qsub}) { unlink "$srcdir/autobuilds.qsub"; }
}


# [qsub] On using qsub. 
#     Rube Goldberg-esque, sorry.
#     We want to run tests under Intel MPI.
#     The Intel MPI environment runs on our cluster nodes, but not on login-eddy.
#     We don't want to directly ssh into a cluster node; we want to go through qsub.
#     qsub -sync y flag makes qsub wait on the command it runs, and return the exit status of that command.
#     The autobuild.pl script writes its "ok" or "FAIL" message to stdout, and hmmer_autobuilds.pl expects to see it.
#     qsub -o autobuilds.qsub captures the autobuild.pl output and writes it to the file $srcdir/autobuilds.qsub
#     qsub's own output (stdout and stderr) we don't want, so send to /dev/null
#     When the qsub is done, we cat the autobuilds.qsub file, to get it back on stdout where hmmer_autobuilds expects it.
#     
#     ssh is a non-interactive shell. It does not run /etc/profile, so it doesn't get Goran's SGE config set up.
#     So we make a more complicated $qsub definition, which first sources Goran's config file.
#
#     Ta-da, that explains the magic incantation: "($qsub ... `$cmd`; cat autobuilds.qsub)"
