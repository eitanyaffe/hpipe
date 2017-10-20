#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <usearch blast6out ifn> <identity ratio> <only target header T|F> <ofn>\n";
	print STDERR "NOTE: for file format see http://www.drive5.com/usearch/manual/blast6out.html\n";
	exit 1;
}

my ($ifn, $identity_ratio, $identity_difference, $strip_str, $ofn) = @ARGV;

my $strip = $strip_str eq "T";

print STDERR "ratio: $identity_ratio\n";

my $appr_lines = apprx_lines($ifn);
print STDERR "traversing file $ifn, with about ".int($appr_lines/1000000)."M lines\n";
open(IN, $ifn) || die $ifn;

my %gene_hits;
my $prev_gene = "";

print STDERR "writing file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "gene\ttop_hit\tidentity_max\tidentity_min\tcoverage\tsecondary_hits\n";

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
	$hit_gene = substr($hit_gene, 0, index($hit_gene, "|"));
    }

    # match details
    my $identity = $f[2];
    my $mlength = $f[3];
    my $qlength = $f[7] - $f[6] + 1;
    my $tlength = $f[9] - $f[8] + 1;
    my $evalue = $f[10];
    my $bitscore = $f[11];

    my $coverage = $mlength / (($tlength < $qlength) ? $tlength : $qlength);
    $coverage = 1 if ($coverage > 1);

    next if (defined($gene_hits{$hit_gene}) && $gene_hits{$hit_gene}->{coverage} > $coverage);

    $gene_hits{$hit_gene} = {};
    $gene_hits{$hit_gene}->{identity} = $identity;
    $gene_hits{$hit_gene}->{coverage} = $coverage;

    if ($prev_gene ne "" && (($prev_gene ne $gene) || eof(IN))) {
	!defined($genes_flag{$prev_gene}) or die "gene already processed: $gene";
	$genes_flag{$prev_gene} = 1;
	my $t_identity = 0;
	my ($t_coverage, $t_hit_gene);
	# first pass find max identity
	foreach my $hit_gene (keys %gene_hits) {
	    if ($gene_hits{$hit_gene}->{identity} > $t_identity) {
		$t_identity = $gene_hits{$hit_gene}->{identity};
		$t_coverage = $gene_hits{$hit_gene}->{coverage};
		$t_hit_gene = $hit_gene;
	    }
	}

	my $identity_threshold1 = 100 - (100 - $t_identity) * $identity_ratio;
	my $identity_threshold2 = $t_identity - $identity_difference;
	my $identity_threshold = $identity_threshold1 > $identity_threshold2 ? $identity_threshold1 : $identity_threshold2;

	# gather all above threshold
	my @refids;
	my $identity_min = $t_identity;
	foreach my $hit_gene (keys %gene_hits) {
	    if ($hit_gene ne $t_hit_gene && $gene_hits{$hit_gene}->{identity} > $identity_threshold) {
		push(@refids, $hit_gene);
		$identity_min = $identity_min < $gene_hits{$hit_gene}->{identity} ? $identity_min : $gene_hits{$hit_gene}->{identity};
	    }
	}
	my $refs = join(";", @refids);
	print OUT $prev_gene, "\t", $t_hit_gene, "\t", $t_identity, "\t", $identity_min, "\t", $t_coverage, "\t", $refs eq "" ? "NA" : $refs, "\n";

	undef %gene_hits;
    }

    $prev_gene = $gene;
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
