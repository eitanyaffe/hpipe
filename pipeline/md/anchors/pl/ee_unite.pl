#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <observed table> <expected table> <cluster bins> <ofn>\n";
	exit 1;
}

my $ifn_obs = $ARGV[0];
my $ifn_expected = $ARGV[1];
my $ifn_bins = $ARGV[2];
my $ofn = $ARGV[3];

#######################################################################################
# read bins
#######################################################################################

my %clusters;
print STDERR "reading bin table: $ifn_bins\n";
open(IN, $ifn_bins) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    $clusters{$f[$h{bin}]} = $f[$h{cluster}];
}
close(IN);

#######################################################################################
# load obs to memory
#######################################################################################

my %matrix;

print STDERR "reading obs into memory: $ifn_obs\n";
open(IN, $ifn_obs) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $cluster_bin1 = $f[$h{cluster_bin1}];
    my $cluster_bin2 = $f[$h{cluster_bin2}];
    my $count = $f[$h{count}];
    $matrix{$cluster_bin1} = {} if (!defined($matrix{$cluster_bin1}));
    $matrix{$cluster_bin1}->{$cluster_bin2} = $count;
}
close(IN);

#######################################################################################
# go over exp table
#######################################################################################

print STDERR "reading expected: $ifn_expected\n";
open(IN, $ifn_expected) || die;
$header = <IN>;
%h = parse_header($header);

print STDERR "creating ofn: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT "cluster1\tcluster2\tobserved\texpected\n";

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $cluster_bin1 = $f[$h{cluster_bin1}];
    my $cluster_bin2 = $f[$h{cluster_bin2}];
    my $exp = $f[$h{value}];

    (defined($clusters{$cluster_bin1}) and defined($clusters{$cluster_bin2})) or die;
    my $cluster1 = $clusters{$cluster_bin1};
    my $cluster2 = $clusters{$cluster_bin2};

    my $obs = (defined($matrix{$cluster_bin1}) && defined($matrix{$cluster_bin1}->{$cluster_bin2})) ?
	$matrix{$cluster_bin1}->{$cluster_bin2} : 0;

    print OUT $cluster1, "\t", $cluster2, "\t", $obs, "\t", $exp, "\n";
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
