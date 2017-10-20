#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util;

if ($#ARGV == -1) {
    print STDERR "usage: $0 <gene_ifn> <set_field> <map_ifn> <source_field> <target_field> <ofn>\n";
	exit 1;
}

my ($gene_ifn, $set_field,
    $map_ifn, $source_field, $target_field,
    $ofn) = @ARGV;

my %genes;
my %sets;

print STDERR "reading gene table: $gene_ifn\n";
open(IN, $gene_ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $set = $f[$h{$set_field}];

    # gene
    if (!defined($genes{$gene})) {
	$genes{$gene} = {};
	$genes{$gene}->{best_hit} = "none";
	$genes{$gene}->{best_identity} = 0;
    }

    # set
    $sets{$set} = {} if (!defined($sets{$set}));
    $sets{$set}->{$gene} = 1;
}
close(IN);

print STDERR "reading map: $map_ifn\n";
open(IN, $map_ifn) || die;
$header = <IN>;
%h = parse_header($header);

my %out;

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $source = $f[$h{$source_field}];
    my $target = $f[$h{$target_field}];
    my $identity = $f[$h{identity}];
    my $coverage = $f[$h{coverage}];

    next if (!defined($genes{$source}));

    if ($genes{$source}->{best_identity} < $identity) {
	$genes{$source}->{best_identity} = $identity;
	$genes{$source}->{best_hit} = $target;
    }
}
close(IN);

# write source table
print STDERR "writing table: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT "set\tgene\tbest_hit\tbest_identity\n";
for my $set (keys %sets) {
for my $gene (keys %{$sets{$set}}) {
    my $hit = $genes{$gene}->{best_hit};
    my $identity = $genes{$gene}->{best_identity};
    print OUT "$set\t$gene\t$hit\t$identity\n";
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

sub read_map
{
    my ($ifn) = @_;
}

