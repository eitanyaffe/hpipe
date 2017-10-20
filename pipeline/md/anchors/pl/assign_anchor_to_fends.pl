#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
    print STDERR "usage: $0 <fend table> <contig table> <anchor table> <gene table> <ofn>\n";
    exit 1;
}

my $ifn = $ARGV[0];
my $contig_table = $ARGV[1];
my $anchor_table = $ARGV[2];
my $gene_table = $ARGV[3];
my $ofn = $ARGV[4];

my $step = 10;

#######################################################################################
# read contig table
#######################################################################################

# contig length table
my %contigs;

print STDERR "reading contig table: $contig_table\n";
open(IN, $contig_table) || die $contig_table;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $length = $f[$h{length}];
    my $nbins = int($length/$step) + 1;
    $contigs{$contig} = {};
    $contigs{$contig}->{coords} = [(0) x $nbins];
    $contigs{$contig}->{anchor} = 0;
}
close(IN);

#######################################################################################
# read anchor table
#######################################################################################

print STDERR "reading anchor table: $anchor_table\n";
open(IN, $anchor_table) || die $anchor_table;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $anchor = $f[$h{anchor}];
    defined($contigs{$contig}) or die;
    $anchor > 0 or die;
    $contigs{$contig}->{anchor} = $anchor;
}
close(IN);

#######################################################################################
# read gene table
#######################################################################################

# contig length table
my %genes;
my $gene_index = 1;

print STDERR "reading gene table: $gene_table\n";
open(IN, $gene_table) || die $gene_table;
$header = <IN>;
%h = parse_header($header);
my $overlap = 0;
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $start = int($f[$h{start}] / $step);
    my $end = int($f[$h{end}] / $step);
    my $gene = $f[$h{gene}];

    !defined($genes{$gene_index}) or die;
    defined($contigs{$contig}) or die $contig;
    $genes{$gene_index} = $gene;

    for (my $i=$start; $i<$end; $i++) {
	if ($contigs{$contig}->{coords}[$i] != 0) {
	    $contigs{$contig}->{coords}[$i] = -1;
	    $overlap++;
	} else {
	    $contigs{$contig}->{coords}[$i] = $gene_index;
	}
    }
    $gene_index++;
}
close(IN);

print "total overlapping bp of genes: $overlap\n";

#######################################################################################
# traverse fends
#######################################################################################

print STDERR "traversing fends file: $ifn\n";
open(IN, $ifn) || die $ifn;
$header = <IN>;
chomp($header);
%h = parse_header($header);

print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;

print OUT $header, "\tgene\tanchor\n";

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $coord = $f[$h{coord}];
    my $bin = int($coord / $step);

    defined($contigs{$contig}) or die;
    my $anchor = $contigs{$contig}->{anchor};
    my $gene_index = $contigs{$contig}->{coords}->[$bin];
    my $gene;
    $gene = "_multi_" if ($gene_index == -1);
    $gene = "_none_" if ($gene_index == 0);
    if ($gene_index > 0) {
	defined($genes{$gene_index}) or die;
	$gene = $genes{$gene_index};
    }

    print OUT $line, "\t", $gene, "\t", $anchor, "\n";
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
