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
	print STDERR " -c: config file\n";
	print STDERR " -s: step\n";
	exit 1;
}

my $command = $ARGV[0];
my $cfg="config/ref/n5.cfg";
my $step="pp_basic";

GetOptions (
    "c=s" => \$cfg,
    "s=s" => \$step);

my $dir = dirname($cfg);
my $fn = basename($cfg);
print "config dir: $dir\n";
print "config file: $fn\n";

if ($command eq "start") {
    system("bash ./docker/drun.sh $dir") == 0 or die;
    system("bash ./docker/dexec.sh $dir make c=config/$fn init") == 0 or die;
}

if ($command eq "stop") {
    system("bash ./docker/dremove.sh $dir") == 0 or die;
}

if ($command eq "run") {
    system("bash ./docker/dexec.sh $dir make c=config/$fn $step") == 0 or die;
}

if ($command eq "dryrun") {
    system("bash ./docker/dexec.sh $dir make c=config/$fn $step -n | grep START") == 0 or die;
}

if ($command eq "debug") {
    system("bash ./docker/dexec.sh $dir bash") == 0 or die;
}
