#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ofn = $ARGV[1];

print STDERR "reading file: $ifn\n";
open(IN, $ifn) || die $ifn;
my $header = <IN>;
my %h = parse_header($header);

my %gene_hits = ();
my $prev_gene = "";

print STDERR "writing file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "gene\tuniref\tidentity\tcoverage\tevalue\tbitscore\n";

# for checking we don't handle gene twice
my %genes_flag;

my $count = 1;
while (my $line = <IN>) {
    print "line: ", $count, "\n" if ($count % 10000000 == 0);
    $count++;

    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{source}];
    my $uniref = $f[$h{target}];
    my $identity = $f[$h{identity}];
    my $coverage = $f[$h{coverage}];
    my $evalue = $f[$h{evalue}];
    my $bitscore = $f[$h{bitscore}];

    my $pline = "$gene\t$uniref\t$identity\t$coverage\t$evalue\t$bitscore\n";

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

