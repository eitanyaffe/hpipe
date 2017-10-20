#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;
use Getopt::Long;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <command> [options]\n";
	print STDERR "command options: \n";
	print STDERR "  start: start hpipe container\n";
	print STDERR "  run: execute step on container\n";
	print STDERR "  stop: stop hpipe container\n";
	print STDERR "options\n";
	print STDERR " -pdir: pipeline dir\n";
	print STDERR " -cdir: config dir\n";
	print STDERR " -cfg: config file\n";
	print STDERR " -step: step\n";
	exit 1;
}

my $command = $ARGV[0];
my $pdir="pipeline";
my $cdir="config/ref";
my $cfg="config/ref/n5.cfg";
my $step="pp_basic";

GetOptions (
    "pdir=s"   => \$pdir,
    "cdir=s"   => \$cdir,
    "cfg=s"   => \$cfg,
    "step=s"   => \$step
    );

if ($command eq "start") {
    system("bash ./docker/drun.sh $pdir $cdir") == 0 or die;
    system("bash ./docker/dexec.sh $cdir make c=$cfg init") == 0 or die;
}

if ($command eq "stop") {
    system("bash ./docker/dremove.sh $cdir") == 0 or die;
}

if ($command eq "run") {
    system("bash ./docker/dexec.sh $cdir make c=$cfg $step") == 0 or die;
}
