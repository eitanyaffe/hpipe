#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Scalar::Util qw(looks_like_number);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <mgm gene table> <gene_prefix> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $prefix = $ARGV[1];
my $ofn = $ARGV[2];

print STDERR "reading file: $ifn\n";
open(IN, $ifn) || die $ifn;

print STDERR "writing file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "gene\tcontig\tstart\tend\tstrand\tlength\tclass\ttruncated\n";

my $count = 0;
my $contig = "";
while (my $line = <IN>) {
    chomp($line);
    $line =~ s/^\s+//;
    my @f = split(/\s+/, $line);

    $contig = $f[3] if ($line ne "" && $f[0] eq "FASTA");
    if (scalar(@f) == 6 && ($f[1] eq "+" or $f[1] eq "-")) {
	my $gene = $f[0];
	my $strand = $f[1];
	my $start = $f[2];
	my $end = $f[3];
	my $length = $f[4];
	my $class = $f[5];

	my $truncated = index($start, "<") != -1 || index($end, ">") != -1 ? "T" : "F";

	$start =~ s/<//g;
	$end =~ s/>//g;

	$gene = $prefix."_".$gene;

	print OUT "$gene\t$contig\t$start\t$end\t$strand\t$length\t$class\t$truncated\n";
	$count++;
    }
}
close(IN);
close(OUT);

print "gene count: $count\n";
