#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <identity ratio> <identity diff> <ofn>\n";
	exit 1;
}

my ($ifn, $identity_ratio, $identity_difference, $ofn) = @ARGV;

print STDERR "max identity ratio: $identity_ratio\n";
print STDERR "max identity difference: $identity_difference\n";

my $appr_lines = apprx_lines($ifn);
print STDERR "traversing file $ifn, with about ".int($appr_lines/1000000)."M lines\n";
open(IN, $ifn) || die $ifn;
my $header = <IN>;
my %h = parse_header($header);

my %gene_hits;
my $prev_gene = "";

print STDERR "writing file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "gene\ttop_uniref\tidentity_max\tidentity_min\tcoverage\tsecondary_unirefs\tsecondary_unirefs_count\n";

# for checking we don't handle gene twice
my %genes_flag;

my $count = 1;
while (my $line = <IN>) {
    print "line: ", $count/1000000, "M\n" if ($count % 100000000 == 0);
    $count++;

    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{source}];
    my $uniref = $f[$h{target}];
    my $identity = $f[$h{identity}];
    my $coverage = $f[$h{coverage}];
    my $evalue = $f[$h{evalue}];
    my $bitscore = $f[$h{bitscore}];

    if ($prev_gene ne "" && (($prev_gene ne $gene) || eof(IN))) {
	!defined($genes_flag{$prev_gene}) or die "gene already processed: $prev_gene";
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
	my $identity_threshold = $identity_threshold1 < $identity_threshold2 ? $identity_threshold1 : $identity_threshold2;

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
	print OUT $prev_gene, "\t", $t_uniref, "\t", $t_identity, "\t", $identity_min, "\t", $t_coverage, "\t", $refs eq "" ? "NA" : $refs, "\t", scalar(@refids), "\n";

	undef %gene_hits;
    }

    next if (defined($gene_hits{$uniref}) && $gene_hits{$uniref}->{identity} > $identity);

    $gene_hits{$uniref} = {};
    $gene_hits{$uniref}->{identity} = $identity;
    $gene_hits{$uniref}->{coverage} = $coverage;

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

