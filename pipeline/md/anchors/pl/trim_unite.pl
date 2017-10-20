#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig bins> <contig anchor table> <anchor field in contig table> <observed> <observed uses contig bins T|F> <expected> <output file>\n";
	exit 1;
}

my $icontig = $ARGV[0];
my $ianchors = $ARGV[1];
my $anchor_field = $ARGV[2];
my $iobs = $ARGV[3];
my $obs_contig_bins = $ARGV[4] eq "T";
my $iexp = $ARGV[5];
my $ofn = $ARGV[6];

#######################################################################################
# read contig table
#######################################################################################

# contig length table
my %contigs;
my %contigs_index;

print STDERR "reading contig bin table: $icontig\n";
open(IN, $icontig) || die $icontig;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $index = $f[$h{index}];

    $contigs_index{$index} = $contig;

    $contigs{$contig} = {};
    $contigs{$contig}->{index} = $index;
    $contigs{$contig}->{anchor} = -1;
    $contigs{$contig}->{contact_count} = {};
    $contigs{$contig}->{contig_count} = {};
}
close(IN);

#######################################################################################
# read contig anchor table
#######################################################################################

my %anchors;

print STDERR "reading contig anchor table: $ianchors\n";
open(IN, $ianchors) || die $ianchors;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $anchor = $f[$h{$anchor_field}];
    $anchors{$contig} = $anchor;
}
close(IN);

#######################################################################################
# read observed table
#######################################################################################

print STDERR "reading obs: $iobs\n";
open(IN, $iobs) || die $iobs;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig_value = $f[$h{contig}];
    my $anchor = $f[$h{anchor}];

    my $contig = "";
    if ($obs_contig_bins) {
	defined($contigs_index{$contig_value}) or die;
	$contig = $contigs_index{$contig_value};
    } else {
	$contig = $contig_value;
    }

    $contigs{$contig}->{total_count}->{$anchor} = $f[$h{contig_total_count}];
    $contigs{$contig}->{median_count}->{$anchor} = $f[$h{contig_median_count}];
}
close(IN);

#######################################################################################
# go over exp table
#######################################################################################

print STDERR "reading exp: $iexp\n";
open(IN, $iexp) || die $iexp;
$header = <IN>;
%h = parse_header($header);

print STDERR "writing file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "contig\tcell\tother_cell\tobserved_contact_count\tobserved_median_count\texpected\n";

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig_i = $f[$h{contig_bin1}];
    my $oanchor = $f[$h{anchor2}];
    my $exp = $f[$h{value}];

    defined($contigs_index{$contig_i}) or die;
    my $contig = $contigs_index{$contig_i};

    my $anchor = defined($anchors{$contig}) ? $anchors{$contig} : -1;

    my $obs_total = defined($contigs{$contig}->{total_count}->{$oanchor}) ? $contigs{$contig}->{total_count}->{$oanchor} : 0;
    my $obs_median = defined($contigs{$contig}->{median_count}->{$oanchor}) ? $contigs{$contig}->{median_count}->{$oanchor} : 0;

    print OUT $contig, "\t", $anchor, "\t", $oanchor, "\t", $obs_total, "\t", $obs_median, "\t", $exp, "\n";

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
