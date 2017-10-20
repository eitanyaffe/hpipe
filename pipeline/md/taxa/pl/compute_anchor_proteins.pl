#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);

if ($#ARGV == -1) {
    print STDERR "usage: $0 <gene anchor table> <uniref top hits> <taxa nodes> <taxa names> <taxa merged> <taxa deleted> <min identity> <min gene cover> <ofn anchor summary> <ofn taxa table> <ofn trees>\n";
	exit 1;
}

my @pnames = ("ifn_genes_anchors",
	     "ifn_uniref",
	     "uniref_tax_lookup",
	     "uniref_table",
	     "ifn_tax_nodes",
	     "ifn_tax_names",
	     "ifn_tax_merged",
	     "ifn_tax_deleted",
	     "min_identity",
	     "min_cover",
	     "mask_ids",
	     "ofn"
    );
my %p;
@p{@pnames} = @ARGV;

print "=============================================\n";
foreach my $key (keys %p) {
    defined($p{$key}) or die "parameter $key not defined (check if all parameters defined";
    print $key, ": ", $p{$key}, "\n";
}
print "=============================================\n";

my @mask_ids = split(" ", $p{mask_ids});
my %masked;
@masked{@mask_ids} = (1) x scalar(@mask_ids);
print "masking tax_ids: ", join(",", keys %masked), "\n";

#######################################################################################
# read uniref desc table
#######################################################################################

my %uniref_desc;

print "reading uniref table: $p{uniref_table}\n";
open(IN, $p{uniref_table}) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line, -1);
    my $uniref = $f[$h{id}];
    my $desc = $f[$h{desc}];
    $uniref_desc{$uniref} = $desc;
}
close(IN);

#######################################################################################
# read taxa table
#######################################################################################

# counts for anchor/tax_id pairs
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

    print "masking $id: $value\n" if (defined($masked{$id}));
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
# read uniref_tax_lookup table
#######################################################################################

my %uniref2tax;

print "reading uniref to taxa lookup table: $p{uniref_tax_lookup}\n";
open(IN, $p{uniref_tax_lookup}) || die;
$header = <IN>;
%h = parse_header($header);
my $count = 0;
my $skipped = 0;
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line, -1);
    my $uniref = $f[$h{uniref}];
    my $ids = $f[$h{tax_ids}];
    $count++;
    if (!defined($ids) or $ids eq "") {
	$skipped++;
	next;
    }
    $uniref2tax{$uniref} = $ids;

    # !!!
#    last if ($count > 10000000);
}
close(IN);
print "number of genes in uniref lookup table: $count\n";
print "skipped genes (no assigned taxa): $skipped\n";

#######################################################################################
# reading uniref hits
#######################################################################################

my %uniref;

print "reading uniref table: $p{ifn_uniref}\n";
open(IN, $p{ifn_uniref}) || die;
$header = <IN>;
%h = parse_header($header);

$count = 0;
$skipped = 0;

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    $uniref{$gene} = {};
    $uniref{$gene}->{uniref} = $f[$h{top_uniref}];
    $uniref{$gene}->{identity} = $f[$h{identity_max}];
}
close(IN);

print "number of uniref hits: ", scalar(keys %uniref), "\n";

#######################################################################################
# read ga table
#######################################################################################

my %anchors;
my %genes;

print "reading gene anchor table: $p{ifn_genes_anchors}\n";
open(IN, $p{ifn_genes_anchors}) || die;
$header = <IN>;
%h = parse_header($header);
$count = 0;
$skipped = 0;
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $anchor = $f[$h{anchor}];

    $anchors{$anchor} = {} if (!defined($anchors{$anchor}));
    $anchors{$anchor}->{$gene} = 1;

    $genes{$gene} = {} if (!defined($genes{$gene}));
    $genes{$gene}->{$anchor} = 1;
}
close(IN);

print "number of anchors: ", scalar(keys %anchors), "\n";
print "number of genes: ", scalar(keys %genes), "\n";

#######################################################################################
# process all gene/achor pairs
#######################################################################################

print "saving ofn: ", $p{ofn}, "\n";
open(OUT, ">", $p{ofn}) || die;
print OUT "anchor\tgene\tshared\tidentity\tuniref\tprotein_desc\tlca_tax_id\tlca_name\tlineage\n";

foreach my $anchor (keys %anchors) {
foreach my $gene (keys %{$anchors{$anchor} }) {

    # gene shared by multiple anchors (2 because of fake global anchor)
    my $shared = scalar(keys %{$genes{$gene}}) > 1 ? "T" : "F";

    # no hits at all
    if (!defined($uniref{$gene})) {
	next;
    }
    my $uniref = $uniref{$gene}->{uniref};
    my $identity = $uniref{$gene}->{identity};

    #################################################
    # convert uniref ids to tax ids
    #################################################

    next if (!defined($uniref2tax{$uniref}));
    my @tax_ids;
    push(@tax_ids, split(";", $uniref2tax{$uniref}));
    my %tids;
    foreach my $id (@tax_ids) {
	$id = $merged{$id} if (defined($merged{$id}));
	next if (defined($deleted{$id}));
	$tids{$id} = 1;
    }
    if (scalar(keys %tids) == 0) {
	$skipped++;
	next;
    }

    #################################################
    # removed masked and inner leaves
    #################################################

    my %stree = %{compute_spanning_tree(keys %tids)};

    my @ids;
    foreach my $id (keys %tids) {
	# masked
	next if (defined($masked{$id}));

	# inner leaf
	defined($stree{$id}) or die;
	next if ($stree{$id} > 0);
	push(@ids, $id);
    }
    if (scalar(@ids) == 0) {
	next;
    }
    my $lca_id = compute_lca(@ids);

    my $lineage = lineage($lca_id);
    my $desc = $uniref_desc{$uniref};
    print OUT $anchor, "\t", $gene, "\t", $shared, "\t", $identity, "\t", $uniref, "\t", $desc, "\t", $lca_id, "\t", $taxa{$lca_id}->{name}, "\t", $lineage, "\n";

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

sub median {
    ( sort { $a <=> $b } @_ )[ int( $#_/2 ) ];
}

sub mean {
    (@_ > 0) ? (sum(@_)/@_) : 0;
}


sub lineage
{
    my ($tax_id) = @_;

    my $node_id = $tax_id;
    my $result = $taxa{$node_id}->{name};
    while ($node_id > 0) {
	defined($taxa{$node_id}) or die $node_id;
	my $parent_id = $taxa{$node_id}->{parent};
	$result = $taxa{$node_id}->{name}." | ".$result if ($node_id != $tax_id);
	$node_id = $parent_id;
    }
    return ($result);
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

sub compute_spanning_tree
{
    my @tax_ids = @_;

    # init tree counts for LCA
    my %tree;
    for my $tax_id (@tax_ids) {
	my $node_id = $tax_id;
	while ($node_id > 0) {
	    defined($taxa{$node_id}) or die;
	    my $parent_id = $taxa{$node_id}->{parent};

 	    $tree{$node_id} = 0 if (!defined($tree{$node_id}));
 	    $tree{$parent_id} = 0 if (!defined($tree{$parent_id}));

	    $tree{$parent_id}++ if ($parent_id > 0);
	    $node_id = $parent_id;
	}
    }
    for my $tax_id (@tax_ids) {
	my $node_id = $tax_id;
	while ($node_id > 0) {
	    defined($taxa{$node_id}) or die;
 	    if (!defined($tree{$node_id})) {
		$tree{$node_id} = 0;
	    }
	    my $parent_id = $taxa{$node_id}->{parent};
	    $tree{$parent_id}++ if ($parent_id > 0);

	    $node_id = $parent_id;
	}
    }
    return \%tree;
}
