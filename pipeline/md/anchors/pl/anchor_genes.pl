#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig ifn> <gene ifn> <ofn>\n";
	exit 1;
}

my ($icontig, $igenes, $ofn) = @ARGV;

#######################################################################################
# read anchor table
#######################################################################################

my %contigs;
print STDERR "reading anchor table: icontig\n";
open(IN, $icontig) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $anchor = $f[$h{anchor}];
    $contigs{$contig} = {} if (!defined($contigs{$contig}));
    $contigs{$contig}->{$anchor} = $line;
}
close(IN);

#######################################################################################
# read gene table
#######################################################################################

print STDERR "generating table: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT "gene\t$header";

print STDERR "reading gene table: $igenes\n";
open(IN, $igenes) || die;
$header = <IN>;
%h = parse_header($header);

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $gene = $f[$h{gene}];

    next if (!defined($contigs{$contig}));
    foreach my $anchor (sort {$a <=> $b} keys %{$contigs{$contig}}) {
	print OUT $gene, "\t", $contigs{$contig}->{$anchor}, "\n";
    }
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
