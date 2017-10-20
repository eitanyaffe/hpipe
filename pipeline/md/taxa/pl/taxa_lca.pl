#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);

if ($#ARGV == -1) {
    print STDERR "usage: $0 <gene group table> <group field> <taxa field> <taxa nodes> <taxa names> <taxa merged> <taxa deleted> <ofn>\n";
	exit 1;
}

my @pnames = ("ifn_genes",
	      "group_field",
	      "taxa_field",
	      "ifn_tax_nodes",
	      "ifn_tax_names",
	      "ifn_tax_merged",
	      "ifn_tax_deleted",
	      "ofn"
    );
my %p;
@p{@pnames} = @ARGV;

print "=============================================\n";
foreach my $key (keys %p) {
    defined($p{$key}) or die "parameter $key not defined (check if all parameters defined)";
    print $key, ": ", $p{$key}, "\n";
}
print "=============================================\n";

#######################################################################################
# read taxa table
#######################################################################################

# counts for group/tax_id pairs
our %taxa;

print "reading taxa names: $p{ifn_tax_names}\n";
open(IN, $p{ifn_tax_names}) || die;
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line, -1);
    my $id = $f[0];
    my $value = $f[2];
    my $type = $f[6];

    next if ($type ne "scientific name");
    !defined($taxa{$id}) or die;
    $taxa{$id} = {};
    $taxa{$id}->{name} = $value;
}
close(IN);

print "reading taxa nodes: $p{ifn_tax_nodes}\n";
open(IN, $p{ifn_tax_nodes}) || die;
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line, -1);

    my $id = $f[0];
    my $parent = $f[2];
    my $level = $f[4];

    defined($taxa{$id}) or die $id;
    $taxa{$id}->{parent} = $parent;
    $taxa{$id}->{level} = $level;
    $taxa{$id}->{children} = {};
}
close(IN);

my %merged;
print "reading taxa merged: $p{ifn_tax_merged}\n";
open(IN, $p{ifn_tax_merged}) || die;
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line, -1);

    my $id = $f[0];
    my $target = $f[2];
    !defined($merged{$id}) or die $id;
    $merged{$id} = $target;
}
close(IN);

my %deleted;
print "reading taxa deleted: $p{ifn_tax_deleted}\n";
open(IN, $p{ifn_tax_deleted}) || die;
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line, -1);

    my $id = $f[0];
    !defined($deleted{$id}) or die $id;
    $deleted{$id} = 1;
}
close(IN);

# make parent of root negative
our $root_id = 1;
$taxa{$root_id}->{parent} == 1 or die;
$taxa{$root_id}->{parent} = -1;

# set children
foreach my $id (keys %taxa) {
    defined($taxa{$id}->{parent}) or die $id;
    my $parent = $taxa{$id}->{parent};
    next if ($parent < 0);
    defined($taxa{$parent}) or die $parent;
    $taxa{$parent}->{children}->{$id} = 1;
}

print "number of taxa: ", scalar(keys %taxa)-1, "\n";

#######################################################################################
# read gene table
#######################################################################################

my %groups;

print "reading gene group table: $p{ifn_genes}\n";
open(IN, $p{ifn_genes}) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $group = $f[$h{$p{group_field}}];
    my $id = $f[$h{$p{taxa_field}}];

    $id = $merged{$id} if (defined($merged{$id}));
    next if (defined($deleted{$id}));
    
    $groups{$group} = {} if (!defined($groups{$group}));
    $groups{$group}->{$id} = 1;
}

#######################################################################################
# output lca per group
#######################################################################################

print "writing tree table: $p{ofn_tree}\n";
open(OUT, ">", $p{ofn_tree}) || die;
print OUT $p{taxa_field}, "\tcount\ttax_id\n";

foreach my $group (keys %groups) {
    my @ids = keys %{$groups{$group}};
    my $lca_id = compute_lca(@ids);
    print OUT $group, "\t", scalar(@ids), "\t", $lca_id, "\n";
}
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

sub compute_lca
{
    my @tax_ids = @_;

    # create tree
    my %lca_tree;
    for my $tax_id (@tax_ids) {
	my $node_id = $tax_id;
	while ($node_id > 0) {
	    defined($taxa{$node_id}) or die;
 	    if (!defined($lca_tree{$node_id})) {
		$lca_tree{$node_id} = {};
		$lca_tree{$node_id}->{count} = 0;
		$lca_tree{$node_id}->{children} = {};
	    }
	    $node_id = $taxa{$node_id}->{parent};
	}
    }

    # set counts
    for my $tax_id (@tax_ids) {
	my $node_id = $tax_id;
	while ($node_id > 0) {
	    defined($taxa{$node_id}) or die;
	    my $parent_id = $taxa{$node_id}->{parent};

 	    defined($lca_tree{$node_id}) or die;
	    $lca_tree{$node_id}->{count}++;

	    if ($parent_id > 0) {
		defined($lca_tree{$parent_id}) or die;
		$lca_tree{$parent_id}->{children}->{$node_id} = 1;
	    }
	    $node_id = $parent_id;
	}
    }

    my $n_tax = scalar(@tax_ids);
    $lca_tree{$root_id}->{count} == $n_tax or die;

    # move down from root until LCA is found
    my $lca_id = $root_id;
    while (1) {
	my $found = 0;
	defined($lca_tree{$lca_id}) or die;
	foreach my $child_id (keys %{$lca_tree{$lca_id}->{children}}) {
	    if ($lca_tree{$child_id}->{count} == $n_tax) {
		$found = 1;
		$lca_id = $child_id;
		last;
	    }
	}
	last if (!$found);
    }
    return ($lca_id);
}
