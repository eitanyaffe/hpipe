#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <usearch blast6out ifn> <identity ratio> <ofn>\n";
	print STDERR "NOTE: for file format see http://www.drive5.com/usearch/manual/blast6out.html\n";
	exit 1;
}

my ($ifn, $identity_ratio, $identity_difference, $ofn) = @ARGV;

print STDERR "ratio: $identity_ratio\n";

my $appr_lines = apprx_lines($ifn);
print STDERR "traversing file $ifn, with about ".int($appr_lines/1000000)."M lines\n";
open(IN, $ifn) || die $ifn;

my %gene_hits;
my $prev_gene = "";

print STDERR "writing file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "gene\ttop_uniref\tidentity_max\tidentity_min\tcoverage\tsecondary_unirefs\n";

# for checking we don't handle gene twice
my %genes_flag;

my $count = 1;
while (my $line = <IN>) {
    print "line: ", $count/1000000, "M\n" if ($count % 100000000 == 0);
    $count++;

    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f1[0];
    my $uniref = $f[1];

    # match details
    my $identity = $f[2];
    my $mlength = $f[3];
    my $qlength = $f[7] - $f[6] + 1;
    my $tlength = $f[9] - $f[8] + 1;
    my $evalue = $f[10];
    my $bitscore = $f[11];

    my $coverage = $mlength / (($tlength < $qlength) ? $tlength : $qlength);
    $coverage = 1 if ($coverage > 1);
    
    next if (defined($gene_hits{$uniref}) && $gene_hits{$uniref}->{coverage} > $coverage);

    $gene_hits{$uniref} = {};
    $gene_hits{$uniref}->{identity} = $identity;
    $gene_hits{$uniref}->{coverage} = $coverage;

    if ($prev_gene ne "" && (($prev_gene ne $gene) || eof(IN))) {
	!defined($genes_flag{$prev_gene}) or die "gene already processed: $gene";
	$genes_flag{$prev_gene} = 1;
	my $t_identity = 0;
	my ($t_coverage, $t_uniref);
	# first pass find max identity
	foreach my $uniref (keys %gene_hits) {
	    if ($gene_hits{$uniref}->{identity} > $t_identity) {
		$t_identity = $gene_hits{$uniref}->{identity};
		$t_coverage = $gene_hits{$uniref}->{coverage};
		$t_uniref = $uniref;
	    }
	}

	my $identity_threshold1 = 100 - (100 - $t_identity) * $identity_ratio;
	my $identity_threshold2 = $t_identity - $identity_difference;
	my $identity_threshold = $identity_threshold1 > $identity_threshold2 ? $identity_threshold1 : $identity_threshold2;

	# gather all above threshold
	my @refids;
	my $identity_min = $t_identity;
	foreach my $uniref (keys %gene_hits) {
	    if ($uniref ne $t_uniref && $gene_hits{$uniref}->{identity} > $identity_threshold) {
		push(@refids, $uniref);
		$identity_min = $identity_min < $gene_hits{$uniref}->{identity} ? $identity_min : $gene_hits{$uniref}->{identity};
	    }
	}
	my $refs = join(";", @refids);
	print OUT $prev_gene, "\t", $t_uniref, "\t", $t_identity, "\t", $identity_min, "\t", $t_coverage, "\t", $refs eq "" ? "NA" : $refs, "\n";

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
