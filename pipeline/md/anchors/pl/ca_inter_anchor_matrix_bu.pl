#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn obs> <ifn exp> <ifn exp contig bins> <ifn exp anchor bins> <ifn contig table> <ofn>";
	exit 1;
}

my $ifn_o = $ARGV[0];
my $ifn_expected = $ARGV[1];
my $ifn_contig_bins = $ARGV[2];
my $ifn_anchor_bins = $ARGV[3];
my $ifn_contig_table = $ARGV[4];
my $ofn = $ARGV[5];

#######################################################################################
# read bins
#######################################################################################

my %contig_bins;
print STDERR "reading contig bin table: $ifn_contig_bins\n";
open(IN, $ifn_contig_bins) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    $contig_bins{$f[$h{bin}]} = $f[$h{contig}];
}
close(IN);

my %anchor_bins;
print STDERR "reading anchor bin table: $ifn_anchor_bins\n";
open(IN, $ifn_anchor_bins) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    $anchor_bins{$f[$h{bin}]} = $f[$h{anchor}];
}
close(IN);

my %contigs;
print STDERR "reading contig table: $ifn_contig_table\n";
open(IN, $ifn_contig_table) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    $contigs{$f[$h{contig}]} = $f[$h{anchor}];
}
close(IN);

#######################################################################################
# go over expected
#######################################################################################

my %matrix;

print STDERR "reading contig expected: $ifn_expected\n";
open(IN, $ifn_expected) || die;
$header = <IN>;
%h = parse_header($header);

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig_bin = $f[$h{contig_bin1}];
    my $anchor_bin = $f[$h{anchor_bin2}];

    defined($contig_bins{$contig_bin}) or die;
    defined($anchor_bins{$anchor_bin}) or die;

    my $contig = $contig_bins{$contig_bin};
    my $anchor = $anchor_bins{$anchor_bin};
    my $exp = $f[$h{value}];

    next if (!defined($contigs{$contig}));
    my $sanchor = $contigs{$contig};

    $matrix{$sanchor} = {} if (!defined($matrix{$sanchor}));
    if (!defined($matrix{$sanchor}->{$anchor})) {
	$matrix{$sanchor}->{$anchor}->{exp} = 0;
	$matrix{$sanchor}->{$anchor}->{obs} = 0;
    }
    $matrix{$sanchor}->{$anchor}->{exp} += $exp;
}
close(IN);

#######################################################################################
# go over observed
#######################################################################################

print STDERR "reading obs: $ifn_o\n";
open(IN, $ifn_o) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $anchor = $f[$h{anchor}];
    my $count = $f[$h{contig_total_count}];

    next if (!defined($contigs{$contig}));
    my $sanchor = $contigs{$contig};

    defined($matrix{$sanchor}) or die;
    defined($matrix{$sanchor}->{$anchor}) or die;
    $matrix{$sanchor}->{$anchor}->{obs} += $count;
}
close(IN);

#######################################################################################
# output
#######################################################################################

print STDERR "creating ofn: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT "anchor_src\tanchor_tgt\tobserved\texpected\n";
foreach my $sanchor (sort {$a <=> $b} keys %matrix) {
foreach my $anchor (sort {$a <=> $b} keys %{$matrix{$sanchor}}) {
    print OUT
	$sanchor, "\t", $anchor, "\t",
	$matrix{$sanchor}->{$anchor}->{obs}, "\t",
	$matrix{$sanchor}->{$anchor}->{exp}, "\n";
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
