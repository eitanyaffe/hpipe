#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <gene/element table> <min identity> <go tree> <ofn.prefix>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $min_identity = $ARGV[1];
my $ifn_go_tree = $ARGV[2];
my $ofn_prefix = $ARGV[3];

###############################################################################################
# go tree
###############################################################################################

my %tree;

print STDERR "reading table: $ifn_go_tree\n";
open(IN, $ifn_go_tree) || die $ifn_go_tree;
my $header = <IN>;
my %h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $id = $f[$h{id}];
    my @parent_ids = split(";", $f[$h{parent_ids}]);

    $tree{$id} = {};
    $tree{$id}->{root} = $f[$h{root}] eq "T";
    $tree{$id}->{desc} = $f[$h{desc}];
    $tree{$id}->{type} = $f[$h{type}];
    $tree{$id}->{genes} = {};
    $tree{$id}->{parent_ids} = {};
    for my $parent_id (@parent_ids) {
	$tree{$id}->{parent_ids}->{$parent_id} = 1;
    }
}
close(IN);
print "GO tree size: ", scalar(keys %tree), "\n";

###############################################################################################
# Gene2GO table
###############################################################################################

print STDERR "reading table: $ifn\n";
open(IN, $ifn) || die $ifn;
$header = <IN>;
%h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $identity = $f[$h{identity}];
    my $GO = $f[$h{GO}];

    next if ($identity < $min_identity);
    my @ids = ($GO);

    while (@ids != 0) {
	my $id = pop(@ids);
	next if ($id eq "NA");
	defined($tree{$id}) or die;
	$tree{$id}->{genes}->{$gene} = 1;

	next if ($tree{$id}->{root});
    }
}
close(IN);

###############################################################################################
# output
###############################################################################################

my $ofn_process = $ofn_prefix."_process";
my $ofn_function = $ofn_prefix."_function";
my $ofn_component = $ofn_prefix."_component";

print STDERR "generating file: $ofn_process\n";
print STDERR "generating file: $ofn_function\n";
print STDERR "generating file: $ofn_component\n";

open(OUT_P, ">", $ofn_process) || die $ofn_process;
print OUT_P "id\tcount\tdesc\n";

open(OUT_F, ">", $ofn_function) || die $ofn_function;
print OUT_F "id\tcount\tdesc\n";

open(OUT_C, ">", $ofn_component) || die $ofn_component;
print OUT_C "id\tcount\tdesc\n";

foreach my $id (keys %tree) {
    $tree{$id}->{gene_count} = scalar keys %{$tree{$id}->{genes}};
}

foreach my $id (sort { $tree{$b}->{gene_count} <=> $tree{$a}->{gene_count} } keys %tree) {
    my $type = $tree{$id}->{type};
    my $gene_count = scalar keys %{$tree{$id}->{genes}};
    next if ($gene_count == 0);
    if ($type eq "biological_process") {
	print OUT_P $id, "\t", $gene_count, "\t", $tree{$id}->{desc}, "\n";
    }
    if ($type eq "molecular_function") {
	print OUT_F $id, "\t", $gene_count, "\t", $tree{$id}->{desc}, "\n";
    }
    if ($type eq "cellular_component") {
	print OUT_C $id, "\t", $gene_count, "\t", $tree{$id}->{desc}, "\n";
    }
}
close(OUT_P);
close(OUT_C);
close(OUT_F);

######################################################################################################
# Subroutines
######################################################################################################


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
