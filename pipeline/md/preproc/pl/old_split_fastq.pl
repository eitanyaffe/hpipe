#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print "usage: $0 <output dir> <max number of reads> <should trim> <trim read offset> <trim read length> <input fastq files>\n";
	exit 1;
}

my $odir = $ARGV[0];
my $max_reads = $ARGV[1];
my $trim = $ARGV[2] eq "T";
my $offset = $ARGV[3];
my $rlen = $ARGV[4];
shift; shift; shift; shift; shift;
my @files = @ARGV;

($#files > 0) or die "no fastq files specified";

my $lines_per_split = $max_reads * 4;

print STDERR "input files:\n", join("\n", @files), "\n";

print sprintf("Splitting with <%d sequences per files, trimming with offset=%d, length=%d\n", $max_reads, $offset, $rlen) if ($trim);
print sprintf("Splitting with <%d sequences per files, no trimming\n", $max_reads) if (!$trim);

my $line;
foreach my $file (@files)
{
	# get file basename
	my $base_name = fileparse($file);

	my $prefix = $odir."/".$base_name.".";

	my $command;
	if ($trim)
	{
	    $command = "cat $file | perl -ane 'BEGIN{\$c=1}if(\$c % 4 == 2 or \$c % 4 == 0){print substr(\$F[0], $offset, $rlen).\"\\n\";}else{print \"\$F[0]\\n\";}\$c++;' | split -a 4 -d -l $lines_per_split - $prefix";
	}
	else
	{
	    $command = "split -a 4 -d -l $lines_per_split $file $prefix";
	}

	print STDERR "splitting input file:\n$command\n";
	system($command) == 0 or die;
}
