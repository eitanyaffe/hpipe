#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;
use Getopt::Long;
use Env;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <command> [options]\n";
	print STDERR "command options: \n";
	print STDERR "  start: start hpipe container\n";
	print STDERR "  stop: stop hpipe container\n";
	print STDERR "  restart: stop and start hpipe container\n";
	print STDERR "  status: check status of docker container\n";
	print STDERR "  step: execute step\n";
	print STDERR "  drystep: show which targets will be generated if step is executed\n";
	print STDERR "options\n";
	print STDERR " -c: config file\n";
	print STDERR " -s: step\n";
	exit 1;
}

my $command = $ARGV[0];
my $cfg = defined($ENV{HPIPE_CONFIG}) ? $ENV{HPIPE_CONFIG} : "config/example/example.cfg";
my $step = defined($ENV{HPIPE_STEP}) ? $ENV{HPIPE_STEP} : "pp_basic";

GetOptions (
    "c=s" => \$cfg,
    "s=s" => \$step);

my $dir = dirname($cfg);
my $fn = basename($cfg);
print "config file: $cfg\n";

if ($command eq "start") {
    msystem("bash ./docker/drun.sh $dir");
    msystem("bash ./docker/dexec.sh $dir make c=config/$fn init");
}

if ($command eq "stop") {
    msystem("bash ./docker/dremove.sh $dir");
}

if ($command eq "restart") {
    msystem("bash ./docker/dremove.sh $dir");
    msystem("bash ./docker/drun.sh $dir");
    msystem("bash ./docker/dexec.sh $dir make c=config/$fn init");
}

if ($command eq "step") {
    print "step: $step\n";
    msystem("bash ./docker/dexec.sh $dir make c=config/$fn $step");
}

if ($command eq "drystep") {
    print "step: $step\n";
    msystem("bash ./docker/dexec.sh $dir make c=config/$fn $step -n | grep START");
}

if ($command eq "debug") {
    msystem("bash ./docker/dexec.sh $dir bash");
}

if ($command eq "status") {
    msystem("bash ./docker/dstatus.sh $dir");
}

sub msystem {
    my $cmd = shift;
    print "# $cmd\n";
    system($cmd) == 0 or die;
}
