#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <input fends> <input mat> <output file>\n";
	exit 1;
}

my $fends_ifn = $ARGV[0];
my $mat_ifn = $ARGV[1];
my $ofn = $ARGV[2];

#############################################################################################
# fend file
#############################################################################################

my %fends;
my %contigs;
my %anchors;

open(IN, $fends_ifn) || die;
print STDERR "Reading file $fends_ifn into hash...\n";
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $fend = $f[$h{fend}];
    my $cluster_bin = $f[$h{cluster_bin}];

    !defined($fends{$fend}) or die "non-unique fend";
    $fends{$fend} = $cluster_bin;
}
close(IN);

#############################################################################################
# mat file
#############################################################################################

my $appr_lines = apprx_lines($mat_ifn);
print STDERR "traversing file $mat_ifn, with about ".int($appr_lines/1000000)."M lines\n";

my %matrix;

my $lcount = 0;
open(IN, $mat_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>)
{
    $lcount++;
    print STDERR "line: $lcount\n" if ($lcount % 1000000 == 0);

    chomp $line;
    my @f = split("\t", $line);
    my $fend1 = $f[$h{fend1}];
    my $fend2 = $f[$h{fend2}];

    next if (!defined ($fends{$fend1}) or !defined ($fends{$fend2}));

    my $cluster_bin1 = $fends{$fend1};
    my $cluster_bin2 = $fends{$fend2};
    next if ($cluster_bin1 eq $cluster_bin2);

    $matrix{$cluster_bin1} = {} if (!defined($matrix{$cluster_bin1}));
    $matrix{$cluster_bin1}->{$cluster_bin2} = 0 if (!defined($matrix{$cluster_bin1}->{$cluster_bin2}));
    $matrix{$cluster_bin1}->{$cluster_bin2}++;

    $matrix{$cluster_bin2} = {} if (!defined($matrix{$cluster_bin2}));
    $matrix{$cluster_bin2}->{$cluster_bin1} = 0 if (!defined($matrix{$cluster_bin2}->{$cluster_bin1}));
    $matrix{$cluster_bin2}->{$cluster_bin1}++;
}
close(IN);

print STDERR "Writing output file: $ofn\n";

open(OUT, ">", $ofn) || die;
print OUT "cluster_bin1\tcluster_bin2\tcount\n";
foreach my $cluster_bin1 (keys %matrix) {
foreach my $cluster_bin2 (keys %{$matrix{$cluster_bin1}}) {
    my $count = $matrix{$cluster_bin1}->{$cluster_bin2};
    print OUT "$cluster_bin1\t$cluster_bin2\t$count\n";
} }
close(OUT);

######################################################################################################
# Subroutines
######################################################################################################

sub median {
  (sort { $a <=> $b } @_ )[ int( $#_/2 ) ];
}

sub min
{
    my ($a, $b) = @_;
    return $a < $b ? $a : $b;
}

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

