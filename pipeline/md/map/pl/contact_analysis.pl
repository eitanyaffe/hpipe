#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig table> <input dir> <max reads> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $idir = $ARGV[1];
my $mreads = $ARGV[2];
my $ofn = $ARGV[3];

my %contigs;


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
    $contigs{$contig} = {};
    $contigs{$contig}->{length} = $f[$h{length}];
    $contigs{$contig}->{contigs} = {};
    $contigs{$contig}->{intra_count} = 0;
    $contigs{$contig}->{inter_count} = 0;
}

###############################################################################################
# traverse all reads
###############################################################################################

print STDERR "Input dir: $idir\n";
my @ifns = <$idir/P*>;

print STDERR "number of input files: ", scalar(@ifns), "\n";

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
	$count++;

	last if ($mreads != -1 && $count > $mreads);
	next if (!defined($contigs{$contig1}) or !defined($contigs{$contig2}));

	if ($contig1 eq $contig2) {
	    $contigs{$contig1}->{intra_count}++;
	} else {
	    $contigs{$contig1}->{inter_count}++;
	    $contigs{$contig2}->{inter_count}++;
	    $contigs{$contig1}->{contigs}->{$contig2} = 0;
	    $contigs{$contig2}->{contigs}->{$contig1} = 0;
	}
    }
    close(IN);
    last if ($mreads != -1 && $count > $mreads);
}
print STDERR "\n";

print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "contig\tlength\tintra_count\tinter_count\tneighbour_count\n";
foreach my $contig (keys %contigs) {
    defined($contigs{$contig}) or die;
    print OUT
	$contig, "\t",
	$contigs{$contig}->{length}, "\t",
	$contigs{$contig}->{intra_count}, "\t",
	$contigs{$contig}->{inter_count}, "\t",
	scalar(keys %{$contigs{$contig}->{contigs}}), "\n";
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
