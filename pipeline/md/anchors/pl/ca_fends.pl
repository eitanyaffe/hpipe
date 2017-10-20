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
my $ofn = $ARGV[3];

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
# traverse fends
#######################################################################################

print STDERR "traversing fends file: $ifn\n";
open(IN, $ifn) || die $ifn;
$header = <IN>;
chomp($header);
%h = parse_header($header);

print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;

print OUT $header, "\tanchor\n";

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $coord = $f[$h{coord}];
    my $bin = int($coord / $step);

    defined($contigs{$contig}) or die $contig;
    my $anchor = $contigs{$contig}->{anchor};
    print OUT $line, "\t", $anchor, "\n";
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
