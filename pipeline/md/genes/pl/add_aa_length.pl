#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn fasta> <ifn gene table> <ofn>\n";
	exit 1;
}

my $ifn_fasta = $ARGV[0];
my $ifn_genes = $ARGV[1];
my $ofn = $ARGV[2];

my %genes;

print STDERR "traversing file: $ifn_fasta\n";
open(IN, $ifn_fasta) || die $ifn_fasta;

my $seq = "";
my $gene = "";
while (my $line = <IN>) {
    chomp($line);
    if (substr($line, 0, 1) eq ">") {
	$genes{$gene} = length($seq) if ($gene ne "");
	$gene = substr($line, 1);
	$seq = "";
    } else {
	$seq .= $line;
    }

}
close(IN);
$genes{$gene} = length($seq) if ($gene ne "");

print STDERR "creating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;

print STDERR "traversing file: $ifn_genes\n";
open(IN, $ifn_genes) || die $ifn_genes;

my $header = <IN>;
my %h = parse_header($header);
chomp($header);
print OUT $header, "\taa_length\n";
my $in = 0;
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    next if (!defined($genes{$gene}));
    print OUT $line, "\t", $genes{$gene}, "\n";
}
close(IN);

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
