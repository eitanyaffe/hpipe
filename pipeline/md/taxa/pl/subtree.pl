#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);

if ($#ARGV == -1) {
    print STDERR "usage: $0 <input table> <taxa nodes> <taxa names> <taxa merged> <taxa deleted> <ofn tree>\n";
	exit 1;
}

my @pnames = ("ifn",
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
    defined($p{$key}) or die "parameter $key not defined (check if all parameters defined";
    print $key, ": ", $p{$key}, "\n";
}
print "=============================================\n";

#######################################################################################
# read taxa table
#######################################################################################

my $count = 0;
my %taxa;
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

#    last if ($count++ > 10000);
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
# go over input table and output trees
#######################################################################################

print "reading input: $p{ifn}\n";
open(IN, $p{ifn}) || die;

print "writing output: $p{ofn}\n";
open(OUT, ">", $p{ofn}) || die;
print OUT "tax_id\tparent_id\troot_id\tname\tlevel\tlineage\n";

my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line, -1);
    my $rid = $f[$h{id}];
    defined($taxa{$rid}) or die;
    
    print "Creating tree for id=", $rid, ", name=", $taxa{$rid}->{name}, "\n";
    my @ids = get_subtree($rid);
    foreach my $id (@ids) {
	my $parent = $taxa{$id}->{parent};
	my $name = $taxa{$id}->{name};
	my $level = $taxa{$id}->{level};
	my $lineage = lineage($id);
	print OUT $id, "\t", $parent, "\t", $rid, "\t", $name, "\t", $level, "\t", $lineage, "\n";
    }
}
close(IN);
close(OUT);

#######################################################################################
# utils
#######################################################################################

sub get_subtree
{
    my ($id) = @_;
    my @ids = ($id);
    defined($taxa{$id}) or die;
    foreach my $cid (keys %{$taxa{$id}->{children}}) {
	push(@ids, get_subtree($cid));
    }
    return (@ids);
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
