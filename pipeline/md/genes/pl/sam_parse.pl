#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <query_table_ifn> <target_table_ifn> <type aa|nt> <ofn>\n";
	exit 1;
}

my ($ifn, $query_table_ifn, $target_table_ifn, $type, $ofn) = @ARGV;

my $length_field = $type eq "aa" ? "aa_length" : "length";

###################################################################################################################
# load gene length into memory
###################################################################################################################

my %query_genes;
print STDERR "reading query gene table: $query_table_ifn\n";
open(IN, $query_table_ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $length = $f[$h{$length_field}];
    $query_genes{$gene} = $length;
}
close(IN);

my %target_genes;
print STDERR "reading target gene table: $target_table_ifn\n";
open(IN, $target_table_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $length = $f[$h{$length_field}];
    $target_genes{$gene} = $length;
}
close(IN);

###################################################################################################################
# parse
###################################################################################################################

my $appr_lines = apprx_lines($ifn);
print STDERR "traversing file $ifn, with about ".int($appr_lines/1000000)."M lines\n";
open(IN, $ifn) || die $ifn;

my %gene_hits;
my $prev_gene = "";

print STDERR "writing file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "source\tsource_strand\ttarget\ttarget_strand\tidentity\tcoverage\tevalue\tbitscore\n";

# for checking we don't handle gene twice
my %genes_flag;

my $count = 1;
while (my $line = <IN>) {
    print "line: ", $count/1000000, "M\n" if ($count % 100000000 == 0);
    $count++;

    # 0:qseqid 1:sseqid 2:pident 3:length 4:mismatch 5:gapopen 6:qstart 7:qend 8:sstart 9:send 10:evalue 11:bitscore
    chomp($line);
    my @f = split("\t", $line);
    my $query_gene = $f[0];
    my $target_gene = $f[1];

    defined($query_genes{$query_gene}) or die $query_gene;
    defined($target_genes{$target_gene}) or die $target_gene;
    my $query_length = $query_genes{$query_gene};
    my $target_length = $target_genes{$target_gene};

    # match details
    my $identity = $f[2];
    my $mlength = $f[3];
    my $qlength = $f[7] - $f[6] + 1;
    my $tlength = $f[9] - $f[8] + 1;
    my $evalue = $f[10];
    my $bitscore = $f[11];

    my $strand_t = $tlength > 0 ? 1 : -1;
    $tlength = ($strand_t == 1) ? $tlength : -$tlength;

    my $strand_q = $qlength > 0 ? 1 : -1;
    $qlength = ($strand_q == 1) ? $qlength : -$qlength;

    $tlength > 0 or die $line;
    $qlength > 0 or die $line;

    my $qcoverage = $qlength / $query_length;
    my $tcoverage = $tlength / $target_length;
    my $coverage = $qcoverage < $tcoverage ? $qcoverage : $tcoverage;

    # correct identity by coverage
    $identity = int(1000 * $identity * $coverage) / 1000;

    print OUT $query_gene, "\t", $strand_q, "\t", $target_gene, "\t", $strand_t, "\t", $identity, "\t", $coverage, "\t", $evalue, "\t", $bitscore, "\n";
}
close(IN);
close(OUT);

#######################################################################################
# utils
#######################################################################################

sub apprx_lines
{
	my ($fn) = @_;
	my $tmp = "/tmp/".$$."_apprx_lines.tmp";
	system("head -n 100000 $fn > $tmp");
	my $size_head = -s $tmp;
	my $size_all = -s $fn;
	return (int($size_all/$size_head*100000));
}

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
