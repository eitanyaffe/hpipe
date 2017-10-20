#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util;
use Algorithm::Combinatorics qw(combinations);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <gene_ifn> <gene_cluster_ifn> <gene field> <genome field> <genome_ifn> <aai_ifn> <min tuple size> <max tuple size> <ofn_prefix>\n";
	exit 1;
}

my ($gene_ifn, $gene_cluster_ifn, $gene_field, $genome_field, $genome_ifn, $aai_ifn, $min_tsize, $max_tsize, $ofn_prefix) = @ARGV;

###################################################################################################################
# genomes table
###################################################################################################################

my %genome2index;
my %index2genome;
my $genome_counter = 1;
print STDERR "reading genome table: $genome_ifn\n";
open(IN, $genome_ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $genome = $f[$h{$genome_field}];
    $genome2index{$genome} = $genome_counter;
    $index2genome{$genome_counter} = $genome;
    # print "genome=$genome, ", "index=", $genome_counter, "\n";
    $genome_counter++;
}
close(IN);

###################################################################################################################
# gene table
###################################################################################################################

my %gene_lengths;
print STDERR "reading gene table: $gene_ifn\n";
open(IN, $gene_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $length = $f[$h{length}];
    $gene_lengths{$gene} = $length;
}
close(IN);

###################################################################################################################
# gene cluster table
###################################################################################################################

my %genes;
print STDERR "reading gene cluster table: $gene_cluster_ifn\n";
open(IN, $gene_cluster_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $ogene = $f[$h{gene}];
    my $gene = $f[$h{$gene_field}];
    my $genome = $f[$h{$genome_field}];

    defined($gene_lengths{$ogene}) or die;
    my $length = $gene_lengths{$ogene};

    defined($genome2index{$genome}) or die;
    my $genome_index = $genome2index{$genome};

    if (!defined($genes{$gene})) {
	$genes{$gene} = {};
	$genes{$gene}->{genomes} = {};
	$genes{$gene}->{length} = 0;
    }
    $genes{$gene}->{genomes}->{$genome_index} = 1;
    $genes{$gene}->{length} = $length if ($genes{$gene}->{length} < $length);
}
close(IN);

###################################################################################################################
# read aai table
###################################################################################################################

our %aai_table;
print STDERR "reading aai table: $aai_ifn\n";
open(IN, $aai_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $aai = $f[$h{identity}];
    my $genome1 = $f[$h{set1}];
    my $genome2 = $f[$h{set2}];
    next if (($genome1 eq $genome2) || ($genome1 eq "NONE") || ($genome2 eq "NONE"));
    defined($genome2index{$genome1}) and defined($genome2index{$genome2}) or die;

    my $genome_index1 = $genome2index{$genome1};
    my $genome_index2 = $genome2index{$genome2};

    my $key = genome_pair_key($genome_index1, $genome_index2);
    $aai_table{$key} = $aai;
}
close(IN);

###################################################################################################################
# go over shared genes
###################################################################################################################

# count genes per combo of genomes
our %share_map;
for (my $size = $min_tsize; $size<=$max_tsize; $size++) {
    $share_map{$size} = {};
}

print "number of genomes: ", scalar keys %genome2index, "\n";
print "number of genes: ", scalar keys %genes, "\n";
for (my $size = $min_tsize; $size<=$max_tsize; $size++) {
    print "processing tuples, k=", $size, "\n";
    foreach my $gene (keys %genes) {
	my @genomes = sort {$a <=> $b} keys %{$genes{$gene}->{genomes}};
	my $length = $genes{$gene}->{length};
	next if (scalar(@genomes) < $size);

	# print "cgene: ", $gene, ", genomes=", join(",", @genomes), "\n";
	my $tuple_iter = combinations(\@genomes, $size);

	# process all tuples
	while (my $tuple_ref = $tuple_iter->next) {
	    my $key = genome_tuple_key($tuple_ref);
	    if (!defined($share_map{$size}->{$key})) {
		$share_map{$size}->{$key} = {};
		$share_map{$size}->{$key}->{gene_count} = 0;
		$share_map{$size}->{$key}->{gene_length} = 0;
	    }
	    $share_map{$size}->{$key}->{gene_count}++;
	    $share_map{$size}->{$key}->{gene_length} += $length;
	}
    }
}

###################################################################################################################
# output files
###################################################################################################################

my @all_genomes = sort {$a <=> $b} keys %index2genome;
for (my $size = $min_tsize; $size<=$max_tsize; $size++) {
    my $ofn = $ofn_prefix.".".$size;
    print STDERR "generating output: $ofn\n";
    open(OUT, ">", $ofn) || die;
    for (my $i = 0; $i<$size; $i++) {
	print OUT "genome".($i+1), "\t";
    }
    print OUT "min_aai\tmax_aai\tgene_count\tgene_length\n";

    my $tuple_iter = combinations(\@all_genomes, $size);
    while (my $tuple_ref = $tuple_iter->next) {
	my $key = genome_tuple_key($tuple_ref);
	foreach my $genome_index (@$tuple_ref) {
	    print OUT $index2genome{$genome_index}, "\t";
	}
	my ($min_aai, $max_aai) = get_aai_values($tuple_ref);
	my $gene_count = defined($share_map{$size}->{$key}) ? $share_map{$size}->{$key}->{gene_count} : 0;
	my $gene_length = defined($share_map{$size}->{$key}) ? $share_map{$size}->{$key}->{gene_length} : 0;
	print OUT $min_aai, "\t", $max_aai, "\t", $gene_count, "\t", $gene_length, "\n";
    }
}

#######################################################################################
# utils
#######################################################################################

sub get_aai_values
{
    my ($tuple_ref) = @_;
    my $pair_iter = combinations($tuple_ref, 2);
    my $min_aai = 100;
    my $max_aai = 0;
    while (my $pair_ref = $pair_iter->next) {
	defined($aai_table{genome_tuple_key($pair_ref)}) or die genome_tuple_key($pair_ref);
	my $aai = $aai_table{genome_tuple_key($pair_ref)};
	$min_aai = $aai if ($aai < $min_aai);
	$max_aai = $aai if ($aai > $max_aai);
    }
    return (($min_aai, $max_aai));
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

sub genome_pair_key
{
	my ($genome_index1, $genome_index2) = @_;
	return($genome_index1 < $genome_index2 ? $genome_index1."_".$genome_index2 : $genome_index2."_".$genome_index1);
}

sub genome_tuple_key
{
	my ($genome_ref) = @_;
	my @sorted = sort { $a <=> $b} @$genome_ref;
	return(join("_", @sorted));
}

sub explode_genome_tuple_key
{
    my ($key) = @_;
    return (split("_", $key));
}

