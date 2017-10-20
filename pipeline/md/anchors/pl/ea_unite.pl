#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <observed table> <gene exp> <contig exp> <gene bins> <contig bins> <anchor bins> <cluster table> <ofn>\n";
	exit 1;
}

my $ifn_o = $ARGV[0];
my $ifn_contigs = $ARGV[1];
my $ifn_contig_bins = $ARGV[2];
my $ifn_anchor_bins = $ARGV[3];
my $ifn_cluster = $ARGV[4];
my $ofn = $ARGV[5];

#######################################################################################
# read tables
#######################################################################################

my %cluster_size;
my %contig2cluster;
print STDERR "reading contig cluster table: $ifn_cluster\n";
open(IN, $ifn_cluster) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $cluster = $f[$h{cluster}];
    $contig2cluster{$contig} = $cluster;
    $cluster_size{$cluster} = 0 if (!defined($cluster_size{$cluster}));
    $cluster_size{$cluster}++;
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
# go over observed
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
    my $obs = $f[$h{contig_total_count}];
    next if (!defined($contig2cluster{$contig}) || $contig2cluster{$contig} == -1);
    my $cluster = $contig2cluster{$contig};

    $matrix{$anchor} = {} if (!defined($matrix{$anchor}));
    if (!defined($matrix{$anchor}->{$cluster})) {
	$matrix{$anchor}->{$cluster} = {};
	$matrix{$anchor}->{$cluster}->{obs} = 0;
	$matrix{$anchor}->{$cluster}->{exp} = 0;
    }
    $matrix{$anchor}->{$cluster}->{obs} += $obs;
}
close(IN);

#######################################################################################
# go over exp
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
    my $exp = $f[$h{value}];

    defined($contigs{$contig_bin}) or die;
    defined($anchors{$anchor_bin}) or die;

    my $contig = $contigs{$contig_bin};
    my $anchor = $anchors{$anchor_bin};
    next if (!defined($contig2cluster{$contig}) || $contig2cluster{$contig} == -1);
    my $cluster = $contig2cluster{$contig};

    next if (!defined($matrix{$anchor}) || !defined($matrix{$anchor}->{$cluster}));
    $matrix{$anchor}->{$cluster}->{exp} += $exp;
}
close(IN);

#######################################################################################
# save matrix
#######################################################################################

print STDERR "writing: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT "cluster\tcluster_size\tanchor\tobs\texp\tscore\n";
foreach my $anchor (sort {$a <=> $b} keys %matrix) {
foreach my $cluster (sort {$a <=> $b} keys %{$matrix{$anchor}}) {
    my $obs = $matrix{$anchor}->{$cluster}->{obs};
    my $exp = $matrix{$anchor}->{$cluster}->{exp};
    my $score = log10($obs / $exp);
    print OUT $cluster, "\t", $cluster_size{$cluster}, "\t", $anchor, "\t", $obs , "\t", $exp, "\t", $score, "\n";
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

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}
