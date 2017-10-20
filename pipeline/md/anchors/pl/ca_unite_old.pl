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
# save to memory all needed keys
#######################################################################################

my %matrix;

print STDERR "reading obs into memory: $ifn_o\n";
open(IN, $ifn_o) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $anchor = $f[$h{anchor}];
    $matrix{$anchor} = {} if (!defined($matrix{$anchor}));
    $matrix{$anchor}->{$contig} = -1;
}
close(IN);

#######################################################################################
# go over contig exp table
#######################################################################################

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

    next if !defined($matrix{$anchor});
    if (defined($matrix{$anchor}->{$contig})) {
	$matrix{$anchor}->{$contig} == -1 or die;
	$matrix{$anchor}->{$contig} = $exp;
    }
}
close(IN);

#######################################################################################
# go again over obs and fill in exp
#######################################################################################

print STDERR "going over: $ifn_o\n";
open(IN, $ifn_o) || die;
$header = <IN>;
%h = parse_header($header);

print STDERR "creating ofn: $ofn\n";
open(OUT, ">", $ofn) || die;
chomp($header);
print OUT $header, "\tcontig_expected\n";

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $anchor = $f[$h{anchor}];
    defined($matrix{$anchor}) or die;
    defined($matrix{$anchor}->{$contig}) or die;
    $matrix{$anchor}->{$contig} != -1 or die "anchor: $anchor, contig=$contig";

    my $contig_exp = $matrix{$anchor}->{$contig};

    print OUT $line, "\t", $contig_exp, "\n";
}
close(IN);
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
