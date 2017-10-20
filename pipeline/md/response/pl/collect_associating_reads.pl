#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);

if ($#ARGV == -1) {
    print STDERR "usage: $0 <>\n";
	exit 1;
}

my @pnames = ("ifn_network", "ifn_clusters", "ifn_anchors", "paired_dir", "ofn");
my %p;
@p{@pnames} = @ARGV;

print "=============================================\n";
foreach my $key (keys %p) {
    defined($p{$key}) or die "parameter $key not defined (check if all parameters defined";
    print $key, ": ", $p{$key}, "\n";
}
print "=============================================\n";

# turn on autoflush, for progress
$|++;

#######################################################################################
# read anchors
#######################################################################################

my %anchors;
print "reading anchors: $p{ifn_anchors}\n";
open(IN, $p{ifn_anchors}) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $anchor = $f[$h{anchor}];
    $anchors{$anchor} = {} if (!defined($anchors{$anchor}));
    $anchors{$anchor}->{$contig} = 1;
}
close(IN);

#######################################################################################
# read clusters
#######################################################################################

my %clusters;
print "reading clusters: $p{ifn_clusters}\n";
open(IN, $p{ifn_clusters}) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $cluster = $f[$h{cluster}];
    $clusters{$cluster} = {} if (!defined($clusters{$cluster}));
    $clusters{$cluster}->{$contig} = 1;
}
close(IN);

#######################################################################################
# read network
#######################################################################################

# map from contig1_contig2 to (anchor,cluster) pairs
my %network;

print "reading network: $p{ifn_network}\n";
open(IN, $p{ifn_network}) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $anchor = $f[$h{anchor}];
    my $cluster = $f[$h{cluster}];

    defined($anchors{$anchor}) && defined($clusters{$cluster}) or die;

    foreach my $anchor_contig (keys %{$anchors{$anchor}}) {
    foreach my $cluster_contig (keys %{$clusters{$cluster}}) {
	my $key1 = $anchor_contig."_".$cluster_contig;
	my $key2 = $cluster_contig."_".$anchor_contig;

	$network{$key1} = {};
	$network{$key1}->{anchor_side} = 1;
	$network{$key1}->{anchor} = $anchor;
	$network{$key1}->{cluster} = $cluster;

	$network{$key2} = {};
	$network{$key2}->{anchor_side} = 2;
	$network{$key2}->{anchor} = $anchor;
	$network{$key2}->{cluster} = $cluster;
    } }
}
close(IN);
print "size of memory network hashtable: ", scalar(keys %network), "\n";

#######################################################################################
# go over all reads
#######################################################################################

print "writing output table: $p{ofn}\n";
open(OUT, ">", $p{ofn}) || die;

my $first = 1;
my @ifns = <$p{paired_dir}/*.pair>;
print "going over input files: ", scalar(@ifns), "\n";
for my $ifn (@ifns) {
    print ".";
    open(IN, $ifn) || die;
    $header = <IN>;
    %h = parse_header($header);

    if ($first) {
	$first = 0;
	chomp($header);
	print OUT "anchor\tcluster\tanchor_side\t", $header, "\n";
    }
    while (my $line = <IN>) {
	chomp($line);
	my @f = split("\t", $line);
	my $contig1 = $f[$h{contig1}];
	my $contig2 = $f[$h{contig2}];
	my $key = $contig1."_".$contig2;
	next if (!defined($network{$key}));
	print OUT $network{$key}->{anchor}, "\t", $network{$key}->{cluster}, "\t", $network{$key}->{anchor_side}, "\t";
	print OUT $line, "\n";
    }
    close(IN);
}
close(OUT);
print "\n";

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
