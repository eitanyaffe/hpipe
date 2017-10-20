#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;
use List::Util;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig table> <cluster table> <contig-anchor matrix> <cluster-cluster table> <ofn>\n";
	exit 1;
}

my $ifn_contig_table = $ARGV[0];
my $ifn_cluster_table = $ARGV[1];
my $ifn_ca = $ARGV[2];
my $ifn_cluster_matrix = $ARGV[3];
my $ofn = $ARGV[4];

###############################################################################################
# read contig table
###############################################################################################

my %contig_length;

print STDERR "reading contig table: $ifn_contig_table\n";
open(IN, $ifn_contig_table) || die $ifn_contig_table;
my $header = <IN>;
my %h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $length = $f[$h{length}];
    $contig_length{$contig} = $length;
}
close(IN);

###############################################################################################
# read cluster table
###############################################################################################

my %contigs;
my %clusters;

print STDERR "reading cluster table: $ifn_cluster_table\n";
open(IN, $ifn_cluster_table) || die $ifn_cluster_table;
$header = <IN>;
%h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $cluster = $f[$h{cluster}];
    $contigs{$contig} = $cluster;
    $clusters{$cluster} = {} if (!defined($clusters{$cluster}));
    defined ($contig_length{$contig}) or die;
    $clusters{$cluster}->{$contig} = {};
    $clusters{$cluster}->{$contig}->{length} = $contig_length{$contig};
}
close(IN);

###############################################################################################
# set length based weights
###############################################################################################

foreach my $cluster (keys %clusters) {
    my $n_contigs = scalar(keys %{$clusters{$cluster}});
    my $sum_length = 0;
    foreach my $contig (keys %{$clusters{$cluster}}) {
	$sum_length += $clusters{$cluster}->{$contig}->{length};
    }
    foreach my $contig (keys %{$clusters{$cluster}}) {
	$clusters{$cluster}->{$contig}->{weight} = $clusters{$cluster}->{$contig}->{length} / $sum_length;
    }
}


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
    my $obs = $f[$h{contig_total_count}];
    my $exp = $f[$h{contig_expected}];

    next if (!defined($contigs{$contig}));
    my $cluster = $contigs{$contig};

    $anchors{$anchor} = {} if (!defined($anchors{$anchor}));
    $anchors{$anchor}->{$cluster} = {} if (!defined($anchors{$anchor}->{$cluster}));
    $anchors{$anchor}->{$cluster}->{$contig} = $obs > 0 ? log10($obs/$exp) : 0;
    # $anchors{$anchor}->{$cluster}->{$contig}->{observed} = $obs;
    # $anchors{$anchor}->{$cluster}->{$contig}->{expected} = $exp;
    # $anchors{$anchor}->{$cluster}->{$contig}->{score} = $obs > 0 ? log10($obs/$exp) : 0;
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
#print OUT "anchor\tcluster\tobserved\texpected\tn.connected.contigs\tn.contigs\n";
print OUT "anchor\tcluster\tobserved\texpected\tn_contigs\tmean_contig_score\n";

my $lcount = 0;
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $anchor = $f[$h{cluster1}];
    my $cluster = $f[$h{cluster2}];
    my $observed = $f[$h{observed}];
    my $expected = $f[$h{expected}];

    # cluster-cluster map contains all clusters on both sides
    next if (!defined($anchors{$anchor}));
    next if (!defined($clusters{$cluster}));

    my $n_contigs = scalar(keys %{$clusters{$cluster}});
    my $sum_score = 0;
    # print "anchor: $anchor, cluster: $cluster, n_contigs: $n_contigs\n";
    foreach my $contig (keys %{$clusters{$cluster}}) {
	# print " contig: $contig\n";
	defined($clusters{$cluster}->{$contig}->{weight}) or die;
	my $weight = $clusters{$cluster}->{$contig}->{weight};
	# print " weight: $weight\n";
	my $score = defined($anchors{$anchor}->{$cluster}) && defined($anchors{$anchor}->{$cluster}->{$contig}) ?
	    $anchors{$anchor}->{$cluster}->{$contig} : 0;
	# print " score: $score\n";
	$sum_score += $weight * $score;
    }
    $n_contigs > 0 or die;
    my $mean_score = $sum_score / $n_contigs;
    print OUT $anchor, "\t", $cluster, "\t", $observed, "\t", $expected, "\t", $n_contigs, "\t", $mean_score, "\n";
    $lcount++;
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

sub median
{
    my ($a1, $a2) = (sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ];
    ($a1 + $a2)/2;
}
