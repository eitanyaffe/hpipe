#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig table> <anchor table> <input dir> <binsize> <output dir>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ifn_anchors = $ARGV[1];
my $idir = $ARGV[2];
my $binsize = $ARGV[3];
my $odir = $ARGV[4];

my %contigs;

print STDERR "binsize=$binsize\n";

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
    $contigs{$contig} = [(0) x $nbins];
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
	my $contig = $f[$h{contig}];
	my $coord = $f[$h{coord}];
	my $bin = floor(($coord-1) / $binsize);

	# sanity checks
	next if (!defined($contigs{$contig}));
	$bin >= 0 && $bin < scalar(@{$contigs{$contig}}) or die $line, "\n", $contig, ":", $coord;

	$contigs{$contig}->[$bin] += 1;
    }
    close(IN);
}
print STDERR "\n";

print STDERR "generating files in: $odir\n";
foreach my $contig (keys %contigs) {
    my $ofn = $odir."/".$contig;
    open(OUT, ">", $ofn) || die $ofn;
    print OUT "coord\tcount\n";
    my $n = scalar(@{$contigs{$contig}});
    for (my $i = 0; $i < $n; $i++) {
	my $coord = $i*$binsize;
	my $count = $contigs{$contig}->[$i];
	print OUT "$coord\t$count\n" if ($count != 0);
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
