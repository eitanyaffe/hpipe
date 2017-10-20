#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Graph;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <gene_ifn> <gene_map_ifn> <identity_threshold> <coverage_threshold> <name prefix> <ofn>\n";
	exit 1;
}

my ($gene_ifn, $gene_map_ifn, $identity_threshold, $coverage_threshold, $prefix, $ofn) = @ARGV;

print STDERR "identity threshold: $identity_threshold\n";
print STDERR "coverage threshold: $coverage_threshold\n";

###################################################################################################################
# generate in memory the blast graph
###################################################################################################################

my $g = Graph->new( undirected => 1 );

print STDERR "constructing gene graph from map: $gene_map_ifn\n";
open(IN, $gene_map_ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    next if ($f[$h{identity}] < $identity_threshold || 100*$f[$h{coverage}] < $coverage_threshold);
    my $source = $f[$h{source}];
    my $target = $f[$h{target}];
    next if ($source eq $target);
    $g->add_edge($source, $target);
}
close(IN);

###################################################################################################################
# go over gene table
###################################################################################################################

print STDERR "reading gene table: $gene_ifn\n";
open(IN, $gene_ifn) || die;
$header = <IN>;
%h = parse_header($header);

print STDERR "generating output gene cluster matrix: $ofn\n";
open(OUT, ">", $ofn) || die;
chomp($header);
print OUT $header,"\tcgene\n";

my %cgenes;

my $count = 0;
my $cindex = scalar($g->connected_components());
print "number of connected components: $cindex\n";
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $cgene = $g->connected_component_by_vertex($gene);
    $cgene = defined($cgene) ? $prefix."_$cgene" : $prefix."_".$cindex++;
    $cgenes{$cgene} = {} if (!defined($cgenes{$cgene}));
    $cgenes{$cgene}->{$gene} = 0;

    print OUT $line, "\t", $cgene, "\n";
    $count++;
}
close(IN);
close(OUT);

print "number of genes: $count\n";
print "number of gene clusters: ".$cindex, "\n";

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
