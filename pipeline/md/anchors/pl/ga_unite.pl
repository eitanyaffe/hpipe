#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <observed table> <gene exp> <contig exp> <gene bins> <contig bins> <anchor bins> <ofn>\n";
	exit 1;
}

my $ifn_o = $ARGV[0];
my $ifn_genes = $ARGV[1];
my $ifn_contigs = $ARGV[2];
my $ifn_gene_bins = $ARGV[3];
my $ifn_contig_bins = $ARGV[4];
my $ifn_anchor_bins = $ARGV[5];
my $ofn = $ARGV[6];

#######################################################################################
# read gene bins
#######################################################################################

my %genes;
print STDERR "reading contig bin table: $ifn_gene_bins\n";
open(IN, $ifn_gene_bins) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    $genes{$f[$h{bin}]} = $f[$h{gene}];
}
close(IN);

my %contigs;
print STDERR "reading contig bin table: $ifn_contig_bins\n";
open(IN, $ifn_contig_bins) || die;
$header = <IN>;
%h = parse_header($header);
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

    my $gene = $f[$h{gene}];
    my $contig = $f[$h{gene_contig}];
    my $anchor = $f[$h{anchor}];

    if (!defined($matrix{$anchor})) {
	$matrix{$anchor} = {};
	$matrix{$anchor}->{genes} = {};
	$matrix{$anchor}->{contigs} = {};
    }
    $matrix{$anchor}->{genes}->{$gene} = -1;
    $matrix{$anchor}->{contigs}->{$contig} = -1;
}
close(IN);

#######################################################################################
# go over gene exp table
#######################################################################################

print STDERR "reading gene expected: $ifn_genes\n";
open(IN, $ifn_genes) || die;
$header = <IN>;
%h = parse_header($header);

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $gene_bin = $f[$h{gene_bin1}];
    my $anchor_bin = $f[$h{anchor_bin2}];

    defined($genes{$gene_bin}) or die;
    defined($anchors{$anchor_bin}) or die;

    my $gene = $genes{$gene_bin};
    my $anchor = $anchors{$anchor_bin};

    my $exp = $f[$h{value}];

    next if !defined($matrix{$anchor});
    if (defined($matrix{$anchor}->{genes}->{$gene})) {
	$matrix{$anchor}->{genes}->{$gene} == -1 or die;
	$matrix{$anchor}->{genes}->{$gene} = $exp;
    }
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
    if (defined($matrix{$anchor}->{contigs}->{$contig})) {
	$matrix{$anchor}->{contigs}->{$contig} == -1 or die;
	$matrix{$anchor}->{contigs}->{$contig} = $exp;
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
print OUT $header, "\tgene_expected\tcontig_expected\n";

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $gene = $f[$h{gene}];
    my $contig = $f[$h{gene_contig}];
    my $anchor = $f[$h{anchor}];
    defined($matrix{$anchor}) or die;
    (defined($matrix{$anchor}->{contigs}->{$contig}) && $matrix{$anchor}->{contigs}->{$contig} != -1) or die;
    (defined($matrix{$anchor}->{genes}->{$gene}) && $matrix{$anchor}->{genes}->{$gene} != -1) or die;

    my $gene_exp = $matrix{$anchor}->{genes}->{$gene};
    my $contig_exp = $matrix{$anchor}->{contigs}->{$contig};

    print OUT $line, "\t", $gene_exp, "\t", $contig_exp, "\n";
}
close(IN);


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
