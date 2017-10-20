#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Scalar::Util qw(looks_like_number);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ofn = $ARGV[1];

print STDERR "reading file: $ifn\n";
open(IN, $ifn) || die $ifn;

print STDERR "writing file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "uniprot\tGO\ttype\tgene_name\tgene_desc\n";

my %term;
my $meta = 1;
my $in_term = 0;
while (my $line = <IN>) {
    chomp($line);
    next if (substr($line, 0, 1) eq "!");

    # UniProtKB	A0A000	moeA5		GO:0003824	GO_REF:0000002	IEA	InterPro:IPR015421|InterPro:IPR015422	F	MoeA5	A0A000_9ACTN|moeA5	protein	taxon:35758	20161029	InterPro

    my @f = split("\t", $line);
    next if ($f[0] ne "UniProtKB");
    my $uniprot = $f[1];
    my $name = $f[2];
    my $GO = $f[4];
    my $type = $f[8];
    my $desc = $f[9];
    print OUT $uniprot, "\t", $GO, "\t", $type, "\t", $name, "\t", $desc, "\n";
}
close(IN);
close(OUT);
