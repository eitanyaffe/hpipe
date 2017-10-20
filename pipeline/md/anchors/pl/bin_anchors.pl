#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <fend table> <output anchor bin file> <output fend>\n";
	exit 1;
}

my ($ifn, $ofn_fends, $ofn_anchor_bins) = @ARGV;

#######################################################################################
# identify and bind anchors
#######################################################################################

print STDERR "traversing fends file: $ifn\n";
open(IN, $ifn) || die $ifn;
my $header = <IN>;
chomp($header);

my %anchors;
my %anchors_back;
my %h = parse_header($header);

my $anchor_bin = 1;
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);

    my $anchor = $f[$h{anchor}];
    next if ($anchor <= 0) || defined($anchors{$anchor});

    $anchors_back{$anchor_bin} = $anchor;
    $anchors{$anchor} = $anchor_bin;
    $anchor_bin++;
}
close(IN);

#######################################################################################
# output anchor bins
#######################################################################################

print STDERR "generating anchor bins: $ofn_anchor_bins\n";
open(OUT, ">", $ofn_anchor_bins) || die $ofn_anchor_bins;
print OUT "bin\tanchor\n";
for my $anchor_bin (sort {$a <=> $b} keys %anchors_back) {
    print OUT $anchor_bin, "\t", $anchors_back{$anchor_bin}, "\n";
}
close(OUT);

#######################################################################################
# second pass - generate output file
#######################################################################################

print STDERR "generating file: $ofn_fends\n";
open(OUT, ">", $ofn_fends) || die $ofn_fends;
print OUT $header, "\tanchor_bin\n";

print STDERR "traversing again fends file: $ifn\n";
open(IN, $ifn) || die $ifn;
$header = <IN>;
%h = parse_header($header);

my $anchor_count = 0;
my $other_count = 0;
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);

    my $anchor = $f[$h{anchor}];
    if ($anchor <= 0) {
	$other_count++;
    } else {
	defined($anchors{$anchor}) or die "anchor not defined";
	print OUT $line, "\t", $anchors{$anchor}, "\n";
	$anchor_count++;
    }
}
close(IN);
close(OUT);

my $sum = $anchor_count + $other_count;
my $ap = int(1000 * $anchor_count / $sum)/10;
my $op = int(1000 * $other_count / $sum)/10;
print "number of anchor fends: ", $anchor_count, "($ap%)\n";
print "number of non-anchor fends: ", $other_count, "($op%)\n";

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
