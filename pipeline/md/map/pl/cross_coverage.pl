#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig table> <input dir> <binsize> <min distance> <max distance> <output dir>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $idir = $ARGV[1];
my $binsize = $ARGV[2];
my $min_d = $ARGV[3];
my $max_d = $ARGV[4];
my $odir = $ARGV[5];

my %contigs;

print STDERR "binsize=$binsize\n";
print STDERR "min distance=$min_d\n";
print STDERR "max distance=$max_d\n";

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
    my $contig_length = $f[$h{length}];
    my $nbins = ceil($contig_length / $binsize) + 1;
    $contigs{$contig} = {};
    $contigs{$contig}->{length} = $contig_length;
    $contigs{$contig}->{bins} = [(0) x $nbins];
}

###############################################################################################
# traverse all reads
###############################################################################################

print STDERR "Input dir: $idir\n";
my @ifns = <$idir/*>;

print STDERR "number of input files: ", scalar(@ifns), "\n";

foreach my $ifn (@ifns) {
    next if (basename($ifn) eq "files");
#    print STDERR "traversing file: $ifn\n";
    print STDERR ".";
    open(IN, $ifn) || die $ifn;
    my $header = <IN>;
    my %h = parse_header($header);
    while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $contig1 = $f[$h{contig1}];
	my $coord1 = $f[$h{coord1}];
	
	my $contig2 = $f[$h{contig2}];
	my $coord2 = $f[$h{coord2}];

	next if ($contig1 ne $contig2);
	my $contig = $contig1;
	next if (!defined($contigs{$contig}));

	my ($coord_l, $coord_r) = $coord1 < $coord2 ? ($coord1,$coord2) : ($coord2,$coord1);
	
	my $dist = $coord_r - $coord_l;
	next if ($dist > $max_d || $dist < $min_d);
	
	my $bin_l = floor(($coord_l-1+$min_d/2) / $binsize);
	my $bin_r = floor(($coord_r-1-$min_d/2) / $binsize);
	next if ($bin_l > $bin_r);
	
	for (my $bin=$bin_l; $bin <= $bin_r; $bin++) {
	    $contigs{$contig}->{bins}->[$bin] += 1;
	}
    }
    close(IN);
}
print STDERR "\n";

print STDERR "generating files in: $odir\n";
foreach my $contig (keys %contigs) {
    my $ofn = $odir."/".$contig;
    open(OUT, ">", $ofn) || die $ofn;
    print OUT "coord\tcount\n";
    my $n = scalar(@{$contigs{$contig}->{bins}});
    for (my $i = 0; $i < $n; $i++) {
	my $coord = $i*$binsize;
	my $count = $contigs{$contig}->{bins}->[$i];
	my $length = $contigs{$contig}->{length};
	next if ($coord <= $max_d || $coord >= ($length-$max_d));
	print OUT "$coord\t$count\n";
    }
    close(OUT);
}

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
