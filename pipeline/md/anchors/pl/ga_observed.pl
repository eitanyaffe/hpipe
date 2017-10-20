#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <input fends> <input mat> <bins per contig> <output file>\n";
	exit 1;
}

my $fends_ifn = $ARGV[0];
my $mat_ifn = $ARGV[1];
my $bins_per_contig = $ARGV[2];
my $ofn = $ARGV[3];

#############################################################################################
# fend file
#############################################################################################

my %fends;
my %genes;
my %contigs;
my %anchors;

open(IN, $fends_ifn) || die;
print STDERR "Reading file $fends_ifn into hash...\n";
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $fend = $f[$h{fend}];
    my $anchor = $f[$h{anchor}];
    my $contig = $f[$h{contig}];
    my $gene = $f[$h{gene}];
    my $coord = $f[$h{coord}];

    !defined($fends{$fend}) or die "non-unique fend";
    $fends{$fend} = {};
    $fends{$fend}->{anchor} = $anchor;
    $fends{$fend}->{contig} = $contig;
    $fends{$fend}->{gene} = $gene;

    # save all anchors
    $anchors{$anchor} = 0;

    if (substr($gene,0,1) ne '_') {
	if (!defined($genes{$gene})) {
	    $genes{$gene} = {};
	    $genes{$gene}->{contig} = $contig;
	    $genes{$gene}->{anchor} = $anchor;
	} else {
	    $genes{$gene}->{contig} eq $contig or die;
	    $genes{$gene}->{anchor} == $anchor or die;
	}
    }

    if (!defined($contigs{$contig})) {
	$contigs{$contig} = {};
	$contigs{$contig}->{fend_count} = 0;
    }

    $fends{$fend}->{contig_index} = $contigs{$contig}->{fend_count};
    $contigs{$contig}->{fend_count}++;
}
close(IN);

#############################################################################################
# mat file
#############################################################################################

my $appr_lines = apprx_lines($mat_ifn);
print STDERR "traversing file $mat_ifn, with about ".int($appr_lines/1000000)."M lines\n";

our %gene_table;
our %contig_table;

my $lcount = 0;
open(IN, $mat_ifn) || die;
$header = <IN>;
%h = parse_header($header);
my $ucount = 0;
while (my $line = <IN>)
{
    $lcount++;
    print STDERR "line: $lcount\n" if ($lcount % 1000000 == 0);

    chomp $line;
    my @f = split("\t", $line);
    my $fend1 = $f[$h{fend1}];
    my $fend2 = $f[$h{fend2}];
    my $count = $f[$h{count}];

    next if (!defined ($fends{$fend1}) or !defined ($fends{$fend2}));

    my $contig1 = $fends{$fend1}->{contig};
    my $contig2 = $fends{$fend2}->{contig};
    next if ($contig1 eq $contig2);

    my $anchor1 = $fends{$fend1}->{anchor};
    my $anchor2 = $fends{$fend2}->{anchor};

    my $gene1 = $fends{$fend1}->{gene};
    my $gene2 = $fends{$fend2}->{gene};

    defined($contigs{$contig1}) or die;
    defined($contigs{$contig2}) or die;

    $contigs{$contig1}->{nbins} = min($bins_per_contig, $contigs{$contig1}->{fend_count}) if (!defined($contigs{$contig1}->{nbins}));
    $contigs{$contig2}->{nbins} = min($bins_per_contig, $contigs{$contig2}->{fend_count}) if (!defined($contigs{$contig2}->{nbins}));

    if ($anchor1 != 0) {
	my $nbins2 = $contigs{$contig2}->{nbins};
	my $contig_bin2 = int(($fends{$fend2}->{contig_index} / $contigs{$contig2}->{fend_count}) * $nbins2);
	add_gene($gene2, $anchor1) if ((substr($gene2,0,1) ne '_'));
	add_contig($contig2, $contig_bin2, $nbins2, $anchor1, $contig1);
    }
    if ($anchor2 != 0) {
	my $nbins1 = $contigs{$contig1}->{nbins};
	my $contig_bin1 = int(($fends{$fend1}->{contig_index} / $contigs{$contig1}->{fend_count}) * $nbins1);
	add_gene($gene1, $anchor2) if ((substr($gene1,0,1) ne '_'));
	add_contig($contig1, $contig_bin1, $nbins1, $anchor2, $contig2);
    }
    $ucount++;
}

close(IN);
print STDERR "Number of inter contig unique contacts: $ucount\n";

print STDERR "Writing output file: $ofn\n";

open(OUT, ">", $ofn) || die;
print OUT "gene\tgene_contig\tgene_anchor\tanchor\tgene_count\tcontig_total_count\tcontig_coverage\tanchor_contig_count\n";
foreach my $gene (keys %genes) {
foreach my $anchor (sort {$a <=> $b} keys %anchors) {
    next if ($anchor == 0);
    my $contig = $genes{$gene}->{contig};
    my $gene_anchor = $genes{$gene}->{anchor};
    defined($contigs{$contig}) or die;

    # skip if there are no contacts between anchor and all of contig
    next if (!defined($contig_table{$contig}->{$anchor}));

    my $nbins = $contigs{$contig}->{nbins};
    my $contig_total_count = $contig_table{$contig}->{$anchor}->{total};
    my $anchor_contig_count = scalar(keys %{$contig_table{$contig}->{$anchor}->{anchor_contigs}});
    my $contig_coverage = scalar(keys %{$contig_table{$contig}->{$anchor}->{bins}}) / $contig_table{$contig}->{$anchor}->{nbins};

    my $gene_count = (defined($gene_table{$gene}) && defined($gene_table{$gene}->{$anchor})) ? $gene_table{$gene}->{$anchor} : 0;
    print OUT "$gene\t$contig\t$gene_anchor\t$anchor\t$gene_count\t$contig_total_count\t$contig_coverage\t$anchor_contig_count\n";
} }
close(OUT);

######################################################################################################
# Subroutines
######################################################################################################

sub median {
  (sort { $a <=> $b } @_ )[ int( $#_/2 ) ];
}

sub min
{
    my ($a, $b) = @_;
    return $a < $b ? $a : $b;
}

sub apprx_lines
{
	my ($fn) = @_;
	my $tmp = "/tmp/".$$."_apprx_lines.tmp";
	system("head -n 100000 $fn > $tmp");
	my $size_head = -s $tmp;
	my $size_all = -s $fn;
	return (int($size_all/$size_head*100000));
}

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

sub add_gene
{
    my ($gene, $anchor) = @_;
    $gene_table{$gene} = {} if (!defined($gene_table{$gene}));

    $gene_table{$gene}->{$anchor} = 0 if (!defined($gene_table{$gene}->{$anchor}));
    $gene_table{$gene}->{$anchor}++;
}

sub add_contig
{
    my ($contig, $bin, $nbins, $anchor, $anchor_contig) = @_;
    $contig_table{$contig} = {} if (!defined($contig_table{$contig}));
    ($bin >= 0 and $bin < $nbins) or die;

    if (!defined($contig_table{$contig}->{$anchor})) {
	$contig_table{$contig}->{$anchor} = {};

	# all count
	$contig_table{$contig}->{$anchor}->{total} = 0;

	# for the contig side
	$contig_table{$contig}->{$anchor}->{bins} = {};
	$contig_table{$contig}->{$anchor}->{nbins} = $nbins;

	# for the anchor side
	$contig_table{$contig}->{$anchor}->{anchor_contigs} = {};
    }
    $contig_table{$contig}->{$anchor}->{total}++;
    $contig_table{$contig}->{$anchor}->{anchor_contigs}->{$anchor_contig} = 0;
    $contig_table{$contig}->{$anchor}->{bins}->{$bin} = 0
}
