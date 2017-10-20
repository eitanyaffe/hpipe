#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Graph;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <cgene_table1_ifn> <cgene_table2_ifn> <map_1_to_2_ifn> <map_2_to_1_ifn> <identity_threshold> <coverage_threshold> <name 1> <name 2> <ofn>\n";
	exit 1;
}

my (
    $cgene_table1_ifn, $cgene_table2_ifn,
    $map_1_to_2_ifn, $map_2_to_1_ifn,
    $identity_threshold, $coverage_threshold,
    $name1, $name2, $ofn) = @ARGV;

print STDERR "identity threshold: $identity_threshold\n";
print STDERR "coverage threshold: $coverage_threshold\n";

###################################################################################################################
# read gene tables into memory
###################################################################################################################

my %genes1;
print STDERR "reading gene table1: $cgene_table1_ifn\n";
open(IN, $cgene_table1_ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $cgene = $f[$h{cgene}];
    $genes1{$gene} = $cgene;
}
close(IN);

my %genes2;
print STDERR "reading gene table2: $cgene_table2_ifn\n";
open(IN, $cgene_table2_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $cgene = $f[$h{cgene}];
    $genes2{$gene} = $cgene;
}
close(IN);

###################################################################################################################
# go over maps
###################################################################################################################

my %map; # from 1 to 2

print STDERR "reading map 1to2: $map_1_to_2_ifn\n";
open(IN, $map_1_to_2_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    next if ($f[$h{identity}] < $identity_threshold || 100*$f[$h{coverage}] < $coverage_threshold);
    my $gene1 = $f[$h{source}];
    my $gene2 = $f[$h{target}];
    defined($genes1{$gene1}) or die;
    defined($genes2{$gene2}) or die;
    my $cgene1 = $genes1{$gene1};
    my $cgene2 = $genes1{$gene2};
    $map{$cgene1} = {} if (!defined($map{$cgene1}));
    $map{$cgene1}->{$cgene2} = 0;
}
close(IN);

print STDERR "reading map 2to1: $map_2_to_1_ifn\n";
open(IN, $map_2_to_1_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    next if ($f[$h{identity}] < $identity_threshold || 100*$f[$h{coverage}] < $coverage_threshold);
    my $gene2 = $f[$h{source}];
    my $gene1 = $f[$h{target}];
    defined($genes2{$gene2}) or die;
    defined($genes1{$gene1}) or die;
    my $cgene1 = $genes1{$gene1};
    my $cgene2 = $genes1{$gene2};
    $map{$cgene1} = {} if (!defined($map{$cgene1}));
    $map{$cgene1}->{$cgene2} = 0;
}
close(IN);

###################################################################################################################
# output graph
###################################################################################################################

print STDERR "generating output graph: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT "cgene_$name1\tcgene_$name2\n";

foreach my $cgene1 (keys %map) {
foreach my $cgene2 (keys %{$map{$cgene1}}) {
    print OUT "$cgene1\t$cgene2\n";
} }
close(OUT);

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
