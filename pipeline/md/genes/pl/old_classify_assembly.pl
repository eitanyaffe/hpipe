#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig_ifn> <ref_gene_ifn> <pre_gene_ifn> <map_ifn> <id_threshold> <ofn genes> <ofn contigs>\n";
	exit 1;
}

my ($contig_ifn, $ref_gene_ifn, $pre_gene_ifn, $map_ifn, $id_threshold, $ofn_gene, $ofn_contig) = @ARGV;

###################################################################################################################
# ref genes
###################################################################################################################

my %ref_genes;

print STDERR "load ref genes: $ref_gene_ifn\n";
open(IN, $ref_gene_ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $contig = $f[$h{contig}];

    $ref_genes{$gene} = $contig;
}

###################################################################################################################
# contigs
###################################################################################################################

my %contigs;

print STDERR "load contigs: $contig_ifn\n";
open(IN, $contig_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    $contigs{$contig} = {};
}

###################################################################################################################
# predicted genes
###################################################################################################################

my %pre_genes;

print STDERR "load pre genes: $pre_gene_ifn\n";
open(IN, $pre_gene_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $coord = $f[$h{start}];
    my $contig = $f[$h{contig}];

    $pre_genes{$gene} = {};
    $pre_genes{$gene}->{identity} = 0;
    $pre_genes{$gene}->{coord} = $coord;

    next if (!defined($contigs{$contig}));
    $contigs{$contig}->{$coord} = $gene;

}

###################################################################################################################
# ref to pre map
###################################################################################################################

my %back_map;

print STDERR "load ref2pre map: $map_ifn\n";
open(IN, $map_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $identity = $f[$h{identity}];
    my $rgene = $f[$h{gene1}];
    my $pgene = $f[$h{gene2}];

    defined($ref_genes{$rgene}) or die;
    my $rcontig = $ref_genes{$rgene};

    if (!defined($back_map{$pgene})) {
	$back_map{$pgene} = {};
	$back_map{$pgene}->{max_identity} = 0;
	$back_map{$pgene}->{max_rcontig} = {};
	$back_map{$pgene}->{rcontigs} = {};
    }

    # keep all
    if (!defined($back_map{$pgene}->{rcontigs}->{$rcontig}) || $back_map{$pgene}->{rcontigs}->{$rcontig} < $identity) {
	$back_map{$pgene}->{rcontigs}->{$rcontig} = $identity;
    }

    # keep max
    if ($back_map{$pgene}->{max_identity} <= $identity) {
	$back_map{$pgene}->{max_identity} = $identity;
	$back_map{$pgene}->{max_rcontig}->{$rcontig} = 1;
    }

    # keep pre gene identity
    defined($pre_genes{$pgene}) or die;
    $pre_genes{$pgene}->{identity} = $identity if ($pre_genes{$pgene}->{identity} < $identity);
}

###################################################################################################################
# classify all contigs
###################################################################################################################

print STDERR "generating contig output table: $ofn_contig\n";
open(OUT_CONTIG, ">", $ofn_contig) || die;
print OUT_CONTIG "contig\tcount\tref_contig\tsupport\tstrong\tweak\tmiss\tspurious\tmin_identity\n";

print STDERR "generating gene output table: $ofn_gene\n";
open(OUT_GENE, ">", $ofn_gene) || die;
print OUT_GENE "contig\tgene\tcoord\tref\ttype\tmax_identity\tref_identity\n";

foreach my $contig (keys %contigs) {

    # first select popular ref contig
    my %rcontigs;
    my $spurious_count = 0;
    foreach my $coord (sort {$a <=> $b} keys %{$contigs{$contig}}) {
	my $pgene = $contigs{$contig}->{$coord};
	defined($pre_genes{$pgene}) or die;
	if ($pre_genes{$pgene}->{identity} < $id_threshold) {
	    $spurious_count++;
	    next;
	}
	next if (!defined($back_map{$pgene}));
	foreach my $rcontig (keys %{$back_map{$pgene}->{max_rcontig}}) {
	    my $ref_identity = $back_map{$pgene}->{rcontigs}->{$rcontig};
	    next if ($ref_identity < $id_threshold);
	    $rcontigs{$rcontig} = 0 if (!defined($rcontigs{$rcontig}));
	    $rcontigs{$rcontig}++;
	}
    }
    my $max_rcontig = "none";
    my $max_count = 0;
    foreach my $rcontig (keys %rcontigs) {
	my $count = $rcontigs{$rcontig};
	if ($max_count < $count) {
	    $max_count = $count;
	    $max_rcontig = $rcontig;
	}
    }

    my $hit_count = 0;
    my $miss_count = 0;
    my $weak_count = 0;
    my $min_ref_identity = -1;

    my $total = scalar(keys %{$contigs{$contig}});

    foreach my $coord (sort {$a <=> $b} keys %{$contigs{$contig}}) {
	my $pgene = $contigs{$contig}->{$coord};
	my $max_identity = $pre_genes{$pgene}->{identity};
	my $ref_identity = 0;
	my $type = "";
	if ($max_rcontig eq "none") {
	    $type = "no_ref";
	} elsif ($max_identity < $id_threshold) {
	    $type = "spurious";
	} elsif (defined($back_map{$pgene}) && defined($back_map{$pgene}->{rcontigs}->{$max_rcontig})) {

	    $ref_identity = $back_map{$pgene}->{rcontigs}->{$max_rcontig};
	    if ($ref_identity > $id_threshold) {
		$hit_count++;
		$type = "OK";
	    } else {
		$weak_count++;
		$type = "weak";
	    }
	    $min_ref_identity = $ref_identity if ($min_ref_identity > $ref_identity || $min_ref_identity == -1);
	} else {
	    $miss_count++;
	    $type = "missclass";
	}

	print OUT_GENE
	    $contig, "\t", $pgene, "\t", $coord, "\t", $max_rcontig, "\t", $type, "\t",
	    $max_identity, "\t", $ref_identity, "\n";
    }

    print OUT_CONTIG
	$contig, "\t", $total, "\t", $max_rcontig, "\t", $max_count, "\t",
	$hit_count, "\t", $weak_count, "\t", $miss_count, "\t", $spurious_count, "\t",
	$min_ref_identity, "\n";
}

close(OUT_CONTIG);
close(OUT_GENE);


#######################################################################################
# utils
#######################################################################################

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

