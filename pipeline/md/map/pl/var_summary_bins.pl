#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig table> <nt_table> <link_table> <variant> <percentage cutoff> <count cutoff> <ofn>\n";
	exit 1;
}

my $ifn_contig = $ARGV[0];
my $ifn_nt = $ARGV[1];
my $ifn_link = $ARGV[2];
my $min_perc = $ARGV[3];
my $min_count = $ARGV[4];
my $ofn = $ARGV[5];

###############################################################################################
# setup contig bins
###############################################################################################

my %contigs;
print STDERR "reading contig table: $ifn_contig\n";
open(IN, $ifn_contig) || die $ifn_contig;
my $header = <IN>;
my %h = parse_header($header);

my ($min_bin,$max_bin) = (1,100);
my @bins = $min_bin..$max_bin;

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $length = $f[$h{length}];
    $contigs{$contig} = {};
    $contigs{$contig}->{length} = $length;
    $contigs{$contig}->{counts} = {};
    for my $bin (@bins) {
	$contigs{$contig}->{counts}->{$bin} = 0;
    }
}

###############################################################################################
# traverse all reads
###############################################################################################

foreach my $ifn (($ifn_nt, $ifn_link)) {
    print STDERR "reading table: $ifn\n";
    open(IN, $ifn) || die $ifn;
    my $header = <IN>;
    my %h = parse_header($header);
    my $lcount = 0;
    while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $contig = $f[$h{contig}];
	my $coord = $f[$h{coord}];
	my $type = $f[$h{type}];
	my $count = $f[$h{count}];
	my $perc = $f[$h{percent}];
	next if (!defined($contigs{$contig}));
	next if ($count < $min_count || $perc < $min_perc);
	next if ($type eq "REF");
	my $bin = int($perc);
	next if ($bin < $min_bin) || ($bin > $max_bin);
	$contigs{$contig}->{counts}->{$bin}++;

	$lcount++;
    }
    close(IN);
}

open(OUT, ">", $ofn) || die $ofn;
print STDERR "generating file: $ofn\n";
print OUT "contig\tlength";
for my $bin (@bins) {
    print OUT "\t", $bin;
}
print OUT "\n";
foreach my $contig (sort keys %contigs) {
    my $length = $contigs{$contig}->{length};
    print OUT $contig, "\t", $length;
    for my $bin (@bins) {
	my $count = $contigs{$contig}->{counts}->{$bin};
	print OUT "\t", $count;
    }
    print OUT "\n";
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
