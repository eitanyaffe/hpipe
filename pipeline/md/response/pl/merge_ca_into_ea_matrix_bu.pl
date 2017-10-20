#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <cluster table> <contig-anchor table> <cluster-cluster table> <ofn>\n";
	exit 1;
}

my $ifn_cluster_table = $ARGV[0];
my $ifn_ca = $ARGV[1];
my $ifn_cluster_matrix = $ARGV[2];
my $ofn = $ARGV[3];

###############################################################################################
# read cluster table
###############################################################################################

my %contigs;
my %clusters;

print STDERR "reading cluster table: $ifn_cluster_table\n";
open(IN, $ifn_cluster_table) || die $ifn_cluster_table;
my $header = <IN>;
my %h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $cluster = $f[$h{cluster}];
    $contigs{$contig} = $cluster;
    $clusters{$cluster} = {} if (!defined($clusters{$cluster}));
    $clusters{$cluster}->{$contig} = 1;
}
close(IN);

###############################################################################################
# read ca matrix
###############################################################################################

my %anchors;
print STDERR "reading ca table: $ifn_ca\n";
open(IN, $ifn_ca) || die $ifn_ca;
$header = <IN>;
%h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $anchor = $f[$h{anchor}];
    defined($contigs{$contig}) or die;
    my $cluster = $contigs{$contig};
    defined($clusters{$cluster}) or die;
    $anchors{$anchor} = {} if (!defined($anchors{$anchor}));
    $anchors{$anchor}->{$cluster} = {} if (!defined($anchors{$anchor}->{$cluster}));
    $anchors{$anchor}->{$cluster}->{$contig} = 1;
}
close(IN);

###############################################################################################
# traverse cluster matrix
###############################################################################################

print STDERR "reading cluster-cluster matrix: $ifn_cluster_matrix\n";
open(IN, $ifn_cluster_matrix) || die $ifn_cluster_matrix;
$header = <IN>;
%h = parse_header($header);

print STDERR "writing anchor-cluster matrix: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "anchor\tcluster\tobserved\texpected\tn.connected.contigs\tn.contigs\n";

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $anchor = $f[$h{cluster1}];
    my $cluster = $f[$h{cluster2}];
    my $observed = $f[$h{observed}];
    my $expected = $f[$h{expected}];

    next if (!defined($anchors{$anchor}));
    next if (!defined($clusters{$cluster}));
    my $n_contigs = scalar(keys %{$clusters{$cluster}});
    my $n_connected_contigs = defined($anchors{$anchor}->{$cluster}) ? scalar(keys %{$anchors{$anchor}->{$cluster}}) : 0;

    print OUT $anchor, "\t", $cluster, "\t", $observed, "\t", $expected, "\t", $n_connected_contigs, "\t", $n_contigs, "\n";
}

close(OUT);
close(IN);

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
