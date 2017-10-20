#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <observed table> <gene exp> <contig exp> <gene bins> <contig bins> <anchor bins> <ofn>\n";
	exit 1;
}

my $ifn_o = $ARGV[0];
my $ifn_contigs = $ARGV[1];
my $ifn_contig_bins = $ARGV[2];
my $ifn_anchor_bins = $ARGV[3];
my $ofn = $ARGV[4];

#######################################################################################
# read gene bins
#######################################################################################

my %contigs;
print STDERR "reading contig bin table: $ifn_contig_bins\n";
open(IN, $ifn_contig_bins) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    $contigs{$f[$h{bin}]} = $f[$h{contig}];
}
close(IN);

my %anchors;
print STDERR "reading anchor bin table: $ifn_anchor_bins\n";
open(IN, $ifn_anchor_bins) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    $anchors{$f[$h{bin}]} = $f[$h{anchor}];
}
close(IN);

#######################################################################################
# go over contig exp table
#######################################################################################

my %matrix;

print STDERR "reading contig expected: $ifn_contigs\n";
open(IN, $ifn_contigs) || die;
$header = <IN>;
%h = parse_header($header);

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig_bin = $f[$h{contig_bin1}];
    my $anchor_bin = $f[$h{anchor_bin2}];

    defined($contigs{$contig_bin}) or die;
    defined($anchors{$anchor_bin}) or die;

    my $contig = $contigs{$contig_bin};
    my $anchor = $anchors{$anchor_bin};

    my $exp = $f[$h{value}];

    $matrix{$anchor} = {} if (!defined($matrix{$anchor}));
    if (!defined($matrix{$anchor}->{$contig})) {
	$matrix{$anchor}->{$contig} = {};
	$matrix{$anchor}->{$contig}->{contig_expected} = $exp;
	$matrix{$anchor}->{$contig}->{any_observed} = "F";
	$matrix{$anchor}->{$contig}->{contig_anchor} = 0;
	$matrix{$anchor}->{$contig}->{contig_total_count} = 0;
	$matrix{$anchor}->{$contig}->{contig_coverage} = 0;
	$matrix{$anchor}->{$contig}->{anchor_contig_count} = 0;
	$matrix{$anchor}->{$contig}->{fend_count} = 0;
    }
}
close(IN);


#######################################################################################
# append observed info
#######################################################################################

print STDERR "appending observed table: $ifn_o\n";
open(IN, $ifn_o) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $anchor = $f[$h{anchor}];
    defined($matrix{$anchor}) && defined($matrix{$anchor}->{$contig}) or die;
    $matrix{$anchor}->{$contig}->{contig_anchor} = $f[$h{contig_anchor}];
    $matrix{$anchor}->{$contig}->{contig_total_count} = $f[$h{contig_total_count}];
    $matrix{$anchor}->{$contig}->{contig_coverage} = $f[$h{contig_coverage}];
    $matrix{$anchor}->{$contig}->{anchor_contig_count} = $f[$h{anchor_contig_count}];
    $matrix{$anchor}->{$contig}->{fend_count} = $f[$h{fend_count}];
    $matrix{$anchor}->{$contig}->{any_observed} = $f[$h{contig_total_count}] > 0 ? "T" : "F";
}
close(IN);

#######################################################################################
# output complete table
#######################################################################################

print STDERR "creating ofn: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT "contig\tcontig_anchor\tanchor\tcontig_total_count\tcontig_coverage\tanchor_contig_count\tfend_count\tcontig_expected\tany_observed\n";

foreach my $anchor (keys %matrix) {
foreach my $contig (keys %{$matrix{$anchor}}) {
    print OUT $contig, "\t", $matrix{$anchor}->{$contig}->{contig_anchor}, "\t", $anchor, "\t";
    print OUT $matrix{$anchor}->{$contig}->{contig_total_count}, "\t", $matrix{$anchor}->{$contig}->{contig_coverage}, "\t";
    print OUT $matrix{$anchor}->{$contig}->{anchor_contig_count}, "\t", $matrix{$anchor}->{$contig}->{fend_count}, "\t";
    print OUT $matrix{$anchor}->{$contig}->{contig_expected}, "\t", $matrix{$anchor}->{$contig}->{any_observed}, "\n";
} }
close(OUT);

#######################################################################################
# utils
#######################################################################################

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
