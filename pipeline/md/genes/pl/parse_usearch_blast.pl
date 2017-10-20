#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <usearch blast6out ifn> <only target header T|F> <ofn>\n";
	print STDERR "NOTE: for file format see http://www.drive5.com/usearch/manual/blast6out.html\n";
	exit 1;
}

my ($ifn, $strip_str, $ofn) = @ARGV;


my $strip = $strip_str eq "T";

my $appr_lines = apprx_lines($ifn);
print STDERR "traversing file $ifn, with about ".int($appr_lines/1000000)."M lines\n";
open(IN, $ifn) || die $ifn;

my %gene_hits;
my $prev_gene = "";

print STDERR "writing file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "source\ttarget\tidentity\tmlength\tevalue\tbitscore\n";

# for checking we don't handle gene twice
my %genes_flag;

my $count = 1;
while (my $line = <IN>) {
    print "line: ", $count/1000000, "M\n" if ($count % 100000000 == 0);
    $count++;

    chomp($line);
    my @f = split("\t", $line);

    my @f1 = split(/\|/,, $f[0]);

    my $gene = $f1[0];
    my $hit_gene = $f[1];

    if ($strip) {
	$gene = substr($gene, 0, index($gene, "|")) if (index($gene, "|") > 0);
	$hit_gene = substr($hit_gene, 0, index($hit_gene, "|")) if (index($hit_gene, "|"));
    }

    # match details
    my $identity = $f[2];
    my $mlength = $f[3];
    my $qlength = $f[7] - $f[6] + 1;
    my $tlength = $f[9] - $f[8] + 1;
    my $evalue = $f[10];
    my $bitscore = $f[11];

    print OUT $gene, "\t", $hit_gene, "\t", $identity, "\t", $mlength, "\t", $evalue, "\t", $bitscore, "\n";
}
close(IN);
close(OUT);

sub apprx_lines
{
	my ($fn) = @_;
	my $tmp = "/tmp/".$$."_apprx_lines.tmp";
	system("head -n 100000 $fn > $tmp");
	my $size_head = -s $tmp;
	my $size_all = -s $fn;
	return (int($size_all/$size_head*100000));
}
