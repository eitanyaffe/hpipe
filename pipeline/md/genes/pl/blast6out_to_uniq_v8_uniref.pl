#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <usearch blast6out ifn> <ofn>\n";
	print STDERR "NOTE: for file format see http://www.drive5.com/usearch/manual/blast6out.html\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ofn = $ARGV[1];

print STDERR "reading file: $ifn\n";
open(IN, $ifn) || die $ifn;

my %gene_hits = ();
my $prev_gene = "";

print STDERR "writing file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "gene\tuniref\tidentity\tmlength\tqlength\ttlength\tevalue\tbitscore\n";

# for checking we don't handle gene twice
my %genes_flag;

my $count = 1;
while (my $line = <IN>) {
    print "line: ", $count, "\n" if ($count % 10000000 == 0);
    $count++;

    chomp($line);
    my @f = split("\t", $line);

    my @f1 = split(/\|/,, $f[0]);

    my $gene = $f1[0];
    my $uniref = $f[1];

    # match details
    my $identity = $f[2];
    my $mlength = $f[3];
    my $qlength = $f[7] - $f[6] + 1;
    my $tlength = $f[9] - $f[8] + 1;
    my $evalue = $f[10];
    my $bitscore = $f[11];

    my $pline = "$gene\t$uniref\t$identity\t$mlength\t$qlength\t$tlength\t$evalue\t$bitscore\n";

    if ($prev_gene ne "" && $prev_gene ne $gene) {
	foreach my $igene (keys %gene_hits) {
	    !defined($genes_flag{$igene}) or die "gene already processed: $igene";
	    $genes_flag{$igene} = 1;
	    my( $max_key, $max_val ) = each %{$gene_hits{$igene}};
	    while (my($k,$v) = each %{$gene_hits{$igene}}) {
		$k > $max_key and ($max_key,$max_val) = ($k,$v);
	    }
	    print OUT $gene_hits{$igene}->{$max_key};
	}
	undef %gene_hits;
	%gene_hits = ();

    }

    $gene_hits{$gene} = {} if (!defined($gene_hits{$gene}));
    $gene_hits{$gene}->{$identity} = $pline;

    $prev_gene = $gene;
}
close(IN);

# last gene
foreach my $igene (keys %gene_hits) {
    !defined($genes_flag{$igene}) or die "gene already processed: $igene";
    $genes_flag{$igene} = 1;
    my( $max_key, $max_val ) = each %{$gene_hits{$igene}};
    while (my($k,$v) = each %{$gene_hits{$igene}}) {
	$k > $max_key and ($max_key,$max_val) = ($k,$v);
    }
    print OUT $gene_hits{$igene}->{$max_key};
}


close(OUT);
