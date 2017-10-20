#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig table> <input dir> <close threshold> <far threshold> <binsize> <max reads> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $idir = $ARGV[1];
my $close_t = $ARGV[2];
my $far_t = $ARGV[3];
my $binsize = $ARGV[4];
my $mreads = $ARGV[5];
my $ofn = $ARGV[6];

my %contigs;

print "close threshold: $close_t\n";
print "far threshold: $close_t\n";
print "binsize: $binsize\n";
print "max reads: $mreads\n" if ($mreads != -1);

###############################################################################################
# setup contig bins
###############################################################################################

print STDERR "reading contig table: $ifn\n";
open(IN, $ifn) || die $ifn;
my $header = <IN>;
my %h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $length = $f[$h{length}];

    $contigs{$contig} = {};
    $contigs{$contig}->{length} = $length;

    # intra
    $contigs{$contig}->{intra_close_count} = 0;
    $contigs{$contig}->{intra_far_count} = 0;

    $contigs{$contig}->{intra_close_area} = $length * $close_t - 0.5 * $close_t**2;
    $contigs{$contig}->{intra_far_area} = 0.5 * ($length - $far_t)**2;

    # inter
    $contigs{$contig}->{contig_bins} = {};
}

###############################################################################################
# traverse all reads
###############################################################################################

print STDERR "Input dir: $idir\n";
my @ifns = <$idir/P*>;

print STDERR "number of input files: ", scalar(@ifns), "\n";
scalar(@ifns) > 0 or die;

my $count = 0;
foreach my $ifn (@ifns) {
#    print STDERR "traversing file: $ifn\n";
    print STDERR ".";
    open(IN, $ifn) || die $ifn;
    my $header = <IN>;
    my %h = parse_header($header);
    while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $contig1 = $f[$h{contig1}];
	my $contig2 = $f[$h{contig2}];
	my $coord1 = $f[$h{back_coord1}];
	my $coord2 = $f[$h{back_coord2}];

	$count++;
	last if ($mreads != -1 && $count > $mreads);

	next if (!defined($contigs{$contig1}) or !defined($contigs{$contig2}));

	if ($contig1 eq $contig2) {
	    my $contig = $contig1;
	    my $dist = abs($coord1 - $coord2);
	    $contigs{$contig}->{intra_close_count}++ if ($dist < $close_t);
	    $contigs{$contig}->{intra_far_count}++ if ($dist > $far_t);
	} else {
	    my $bin1 = $contig1."_".int($coord1 / $binsize);
	    my $bin2 = $contig2."_".int($coord2 / $binsize);
	    $contigs{$contig1}->{contig_bins}->{$bin2} = 0;
	    $contigs{$contig2}->{contig_bins}->{$bin1} = 0;
	}
    }
    close(IN);
    last if ($mreads != -1 && $count > $mreads);
}
print STDERR "\n";

print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "contig\tlength\tintra_close_count\tintra_close_area\tintra_far_count\tintra_far_area\tinter_count\n";
foreach my $contig (keys %contigs) {
    defined($contigs{$contig}) or die;
    print OUT
	$contig, "\t",
	$contigs{$contig}->{length}, "\t",
	$contigs{$contig}->{intra_close_count}, "\t",
	$contigs{$contig}->{intra_close_area}, "\t",
	$contigs{$contig}->{intra_far_count}, "\t",
	$contigs{$contig}->{intra_far_area}, "\t",
	scalar(keys %{$contigs{$contig}->{contig_bins}}), "\n";
}
close(OUT);

######################################################################################################
# Subroutines
######################################################################################################


sub parse_header
{
	my ($header) = @_;
	chomp($header);
	my @f = split("\t", $header);
	my %result;
	for (my $i = 0; $i <= $#f; $i++) {
		$result{$f[$i]} = $i;
	}
	return %result;
}
