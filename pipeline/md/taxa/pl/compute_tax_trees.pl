#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);

if ($#ARGV == -1) {
    print STDERR "usage: $0 <gene set table> <set field> <uniref top hits> <taxa nodes> <taxa names> <taxa merged> <taxa deleted> <min identity> <min gene cover> <mask cag T|F> <mask ids> <ofn set summary> <ofn taxa table> <ofn trees>\n";
	exit 1;
}

my @pnames = ("ifn_genes_sets",
	      "set_field",
	      "ifn_uniref",
	      "uniref_tax_lookup",
	      "ifn_tax_nodes",
	      "ifn_tax_names",
	      "ifn_tax_merged",
	      "ifn_tax_deleted",
	      "min_identity",
	      "min_cover",
	      "mask_cag",
	      "mask_ids",
	      "ofn_summary",
	      "ofn_taxa",
	      "ofn_tree"
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
# read taxa table
#######################################################################################

# counts for set/tax_id pairs
our %taxa;

print "masking CAG nodes: $p{mask_cag}\n";
my $mask_cag = $p{mask_cag} eq "T";

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

    $masked{$id} = 1 if ($mask_cag && index($value, " CAG:") != -1);
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
my $header = <IN>;
my %h = parse_header($header);
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
    my $t_uniref = $f[$h{top_uniref}];
    my $other_unirefs = $f[$h{secondary_unirefs}];

    my @unirefs = ($t_uniref);
    push(@unirefs, split(";", $other_unirefs)) if ($other_unirefs ne "NA");

    !defined($uniref{$gene}) or die;

    $uniref{$gene} = {};
    $uniref{$gene}->{unirefs} = [@unirefs];
    $uniref{$gene}->{identity} = $f[$h{identity_max}];
    $uniref{$gene}->{coverage} = $f[$h{coverage}];
}
close(IN);

print "number of uniref hits: ", scalar(keys %uniref), "\n";

#######################################################################################
# read ga table
#######################################################################################

our %sets;
my %genes;

my @thresholds = (70,80,90,98,100);

print "reading gene set table: $p{ifn_genes_sets}\n";
open(IN, $p{ifn_genes_sets}) || die;
$header = <IN>;
%h = parse_header($header);
$count = 0;
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $set = $f[$h{$p{set_field}}];

    # skip if zero (relevant for contig_anchor field)
    # next if ($set == 0);

    if (!defined($sets{$set})) {
	$sets{$set} = {};
	$sets{$set}->{genes} = {};

	# general counts per set
	$sets{$set}->{count} = {};
	$sets{$set}->{count}->{no_hit} = 0;
	$sets{$set}->{count}->{no_taxa} = 0;
	$sets{$set}->{count}->{poor_hit} = 0;
	$sets{$set}->{count}->{masked} = 0;

	$sets{$set}->{count}->{hits} = {};
	foreach my $threshold (@thresholds) {
	    $sets{$set}->{count}->{hits}->{$threshold} = 0;
	}

	# all taxa of an set including internal nodes
	$sets{$set}->{node} = {};
    }
    $sets{$set}->{genes}->{$gene} = 1;

    if (!defined($genes{$gene})) {
	$genes{$gene} = {};
	$genes{$gene}->{sets} = {};
    }
    $genes{$gene}->{sets}->{$set} = 1;
}
close(IN);

print "number of sets: ", scalar(keys %sets), "\n";
print "number of genes: ", scalar(keys %genes), "\n";

#######################################################################################
# process all gene/achor pairs
#######################################################################################

$skipped = 0;
$count = 0;
foreach my $set (keys %sets) {
foreach my $gene (keys %{$sets{$set}->{genes}}) {
    # no hits at all
    if (!defined($uniref{$gene})) {
	$sets{$set}->{count}->{no_hit}++;
	$skipped++;
	next;
    }

    my $identity = $uniref{$gene}->{identity};
    my $coverage = $uniref{$gene}->{coverage};
    my @unirefs = @{$uniref{$gene}->{unirefs}};
    scalar(@unirefs) > 0 or die;

    #################################################
    # convert uniref ids to tax ids
    #################################################

    my @tax_ids;
    foreach my $uniref (@unirefs) {
	next if (!defined($uniref2tax{$uniref}));
	push(@tax_ids, split(";", $uniref2tax{$uniref}));
    }
    my %tids;
    foreach my $id (@tax_ids) {
	$id = $merged{$id} if (defined($merged{$id}));
	next if (defined($deleted{$id}));
	$tids{$id} = 1;
    }
    if (scalar(keys %tids) == 0) {
	$sets{$set}->{count}->{no_taxa}++;
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
	$sets{$set}->{count}->{masked}++;
	$skipped++;
	next;
    }

    #################################################
    # hit quality
    #################################################

    if ($coverage < $p{min_cover}/100.0 || $identity < $p{min_identity}) {
	$sets{$set}->{count}->{poor_hit}++;
	$skipped++;
	next;
    }

    if ($identity < 70) {
	$sets{$set}->{count}->{hits}->{70}++;
    } elsif ($identity < 80) {
	$sets{$set}->{count}->{hits}->{80}++;
    } elsif ($identity < 90) {
	$sets{$set}->{count}->{hits}->{90}++;
    } elsif ($identity < 98) {
	$sets{$set}->{count}->{hits}->{98}++;
    } else {
	$sets{$set}->{count}->{hits}->{100}++;
    }

    $count++;

    #################################################
    # update weighted gene counts and identity
    #################################################

    my $step = 1.0 / scalar(@ids);

    foreach my $id (@ids) {
	check_create_node($set, $id);
	!defined($sets{$set}->{node}->{$id}->{genes}->{$gene}) or die;
	$sets{$set}->{node}->{$id}->{genes}->{$gene} = 1;
	$sets{$set}->{node}->{$id}->{weight_count} += $step;
	push(@{$sets{$set}->{node}->{$id}->{identity}}, $identity);
	push(@{$sets{$set}->{node}->{$id}->{coverage}}, $coverage);
    }

    my $lca_id = compute_lca(@ids);

    check_create_node($set, $lca_id);
    $sets{$set}->{node}->{$lca_id}->{lca_count}++;
} }

$count > 0 or die "no genes with taxa found";

print "gene/set pairs with assigned taxa: $count\n";
print "skipped gene/set pairs: $skipped\n";

#######################################################################################
# summary table
#######################################################################################

print "writing summary table: $p{ofn_summary}\n";
open(OUT, ">", $p{ofn_summary}) || die;
print OUT "anchor\tno_hit\tno_taxa\tmasked\tpoor_hit";
foreach my $threshold (@thresholds) {
    print OUT "\t", "hit_", $threshold;
}
print OUT "\n";

foreach my $set (sort {$a <=> $b} keys %sets) {
    my $no_hit = $sets{$set}->{count}->{no_hit};
    my $no_taxa = $sets{$set}->{count}->{no_taxa};
    my $poor = $sets{$set}->{count}->{poor_hit};
    my $masked = $sets{$set}->{count}->{masked};
    print OUT $set, "\t", $no_hit, "\t", $no_taxa, "\t", $masked, "\t", $poor;
    foreach my $threshold (@thresholds) {
	print OUT "\t", $sets{$set}->{count}->{hits}->{$threshold};
    }
    print OUT "\n";
}
close(OUT);

#######################################################################################
# mark tree and output set trees
#######################################################################################

print "writing tree table: $p{ofn_tree}\n";
open(OUT, ">", $p{ofn_tree}) || die;
print OUT "anchor\ttax_id\tis_leaf\tgene_count\tweight\tlca_count\tidentity\tcoverage\n";

my %visited;

foreach my $set (sort {$a <=> $b} keys %sets) {

    my %counts;
    my $wtotal = 0;

    # print "set: ", $set, ", nodes=", scalar(keys %{$sets{$set}->{node}}), "\n";

    # counts and identities
    foreach my $leaf_id (keys %{$sets{$set}->{node}}) {
	my $wcount = $sets{$set}->{node}->{$leaf_id}->{weight_count};
	my $lcount = $sets{$set}->{node}->{$leaf_id}->{lca_count};
	my @genes = keys %{$sets{$set}->{node}->{$leaf_id}->{genes}};
	my $acount = scalar(@genes);

	$wtotal += $wcount;
	my $id = $leaf_id;
	while ($id > 0) {
	    defined($taxa{$id}) or die;
	    $visited{$id} = 1;
	    if (!defined($counts{$id})) {
		$counts{$id} = { wcount => 0, lcount => 0, genes => {}, identity => [], coverage => [] };
	    }
	    foreach my $gene (@genes) { $counts{$id}->{genes}->{$gene} = 1; }

	    $counts{$id}->{wcount} += $wcount;
	    $counts{$id}->{lcount} += $lcount;

	    # identity
	    my @leaf_iden = @{$sets{$set}->{node}->{$leaf_id}->{identity}};
	    push(@{$counts{$id}->{identity}}, @leaf_iden);

	    # coverage
	    my @leaf_cov = @{$sets{$set}->{node}->{$leaf_id}->{coverage}};
	    push(@{$counts{$id}->{coverage}}, @leaf_cov);

	    $id = $taxa{$id}->{parent};
	}
    }

    if ($wtotal == 0) {
	print "skipping set because no genes with taxa found, id: $set\n";
	next;
    }
    $counts{$root_id}->{wcount} == $wtotal or die;

    for my $id (sort {$a <=> $b} keys %counts) {
	my @genes = keys %{$counts{$id}->{genes}};
	my $acount = scalar(@genes);
	my $wcount = $counts{$id}->{wcount};
	my $lcount = $counts{$id}->{lcount};
	my $is_leaf = (scalar(keys %{$taxa{$id}->{children}}) == 0) ? "T" : "F";

	$wcount = $lcount if ($wcount < $lcount);
	$wcount = $acount if ($wcount > $acount);

	my $identity = median(@{$counts{$id}->{identity}});
	my $coverage = median(@{$counts{$id}->{coverage}});

	$acount >= $wcount && $wcount >= $lcount or die;

	print OUT $set, "\t", $id, "\t", $is_leaf, "\t", $acount, "\t", $wcount, "\t", $lcount,
			  "\t", $identity, "\t", $coverage, "\n";
    }
}
close(OUT);

#######################################################################################
# limit taxa table to visited taxa
#######################################################################################

print "Restricting taxa table to visited taxas: ", scalar(keys %visited), "\n";

print "writing restricted taxa table: $p{ofn_taxa}\n";
open(OUT, ">", $p{ofn_taxa}) || die;
print OUT "tax_id\tparent_id\tname\tlevel\tlineage\n";
foreach my $id (keys %taxa) {
    next if (!defined($visited{$id}));
    my $parent = $taxa{$id}->{parent};
    my $name = $taxa{$id}->{name};
    my $level = $taxa{$id}->{level};
    my $lineage = lineage($id);
    print OUT $id, "\t", $parent, "\t", $name, "\t", $level, "\t", $lineage, "\n";
}
close(IN);
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

sub check_create_node
{
    my ($set, $tax_id) = @_;
    if (!defined($sets{$set}->{node}->{$tax_id})) {
	$sets{$set}->{node}->{$tax_id} = {};
	$sets{$set}->{node}->{$tax_id}->{identity} = [];
	$sets{$set}->{node}->{$tax_id}->{coverage} = [];
	$sets{$set}->{node}->{$tax_id}->{weight_count} = 0;
	$sets{$set}->{node}->{$tax_id}->{lca_count} = 0;
	$sets{$set}->{node}->{$tax_id}->{genes} = {};
    }
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
