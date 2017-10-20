#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);
use POSIX qw( strftime );
use File::stat;
use File::Basename;

if ($#ARGV == -1) {
    print STDERR "usage: $0 <fasta file> <contig table> <gene table> <contig anchor table> <type A|U|X|nS> <gene gap> <odir>\n";
    print STDERR "   type: A=Anchor, U=Union, X=Accessory (Union-Anchor), S=only shared, nS=no shared (Union-shared_contigs)\n";
    exit 1;
}

my @pnames = (
    "fasta_ifn",
    "contig_ifn",
    "gene_ifn",
    "ca_ifn",
    "type",
    "gene_gap",
    "odir"
    );
my %p;
@p{@pnames} = @ARGV;

print "=============================================\n";
foreach my $key (keys %p) {
    defined($p{$key}) or die "parameter $key not defined (check if all parameters defined";
    print $key, ": ", $p{$key}, "\n";
}
print "=============================================\n";

my $fasta_ifn = $p{fasta_ifn};
my $contig_ifn = $p{contig_ifn};
my $genes_ifn = $p{gene_ifn};
my $type = $p{type};
my $gene_gap = $p{gene_gap};
my $ca_ifn = $p{ca_ifn};

$type eq "A" || $type eq "U" || $type eq "X" || $type eq "S" || $type eq "nS" or die "unknown type: $type";

system "rm -rf $p{odir}";
system "mkdir -p $p{odir}";

#######################################################################################
# read gene fasta into memory
#######################################################################################

my %genes;

print "reading fasta into memory: $fasta_ifn\n";
open(IN, $fasta_ifn) || die;
my $gene = "";
my $seq = "";
while (my $line = <IN>) {
    chomp($line);
    if (substr($line,0,1) eq ">") {
	$genes{$gene} = $seq if ($seq ne "");
	$gene = substr($line, 1);
	$seq = "";
    } else {
	$seq .= $line;
    }
}
$genes{$gene} = $seq if ($seq ne "");
close(IN);

#######################################################################################
# read contig table
#######################################################################################

my %contigs;

print "reading contig table: $contig_ifn\n";
open(IN, $contig_ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $length = $f[$h{length}];
    !defined($contigs{$contig}) or die;
    $contigs{$contig} = {};
    $contigs{$contig}->{length} = $length;
    $contigs{$contig}->{genes} = {};
    $contigs{$contig}->{anchors} = {};
}
close(IN);

#######################################################################################
# read gene table
#######################################################################################

my $total_genes = 0;
my $trimmed_genes = 0;

print "reading gene table: $genes_ifn\n";
open(IN, $genes_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $contig = $f[$h{contig}];
    my $start = $f[$h{start}];
    my $end = $f[$h{end}];
    next if (!defined($contigs{$contig}));

    $total_genes++;
    if (($start < $gene_gap) || $end > ($contigs{$contig}->{length} - $gene_gap)) {
	$trimmed_genes++;
	next;
    }
    $contigs{$contig}->{genes}->{$gene} = 1;
}
close(IN);
print "total genes: $total_genes\n";
print "trimmed genes: $trimmed_genes\n";

#######################################################################################
# read ca table
#######################################################################################

# first compute multiplicity
print "reading ca table: $ca_ifn\n";
open(IN, $ca_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $anchor = $f[$h{anchor}];
    my $contig = $f[$h{contig}];
    next if (!defined($contigs{$contig}));
    $contigs{$contig}->{anchors}->{$anchor} = 1;
}
close(IN);
foreach my $contig (keys %contigs) {
    $contigs{$contig}->{is_multi_anchor} = scalar(keys %{$contigs{$contig}->{anchors}}) > 1;
}

my %anchors;

print "reading ca table: $ca_ifn\n";
open(IN, $ca_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $anchor = $f[$h{anchor}];
    my $contig = $f[$h{contig}];
    my $is_anchor = $f[$h{anchor}] eq $f[$h{contig_anchor}];

    next if (!defined($contigs{$contig}));
    my $is_multi_anchor = $contigs{$contig}->{is_multi_anchor};

    next if ($type eq "A" && !$is_anchor);
    next if ($type eq "X" && $is_anchor);
    next if ($type eq "nS" && $is_multi_anchor);
    next if ($type eq "S" && !$is_multi_anchor);

    for my $gene (keys %{$contigs{$contig}->{genes}}) {
	my $seq = $genes{$gene};
	$anchors{$anchor} = {} if (!defined($anchors{$anchor}));
	$anchors{$anchor}->{$gene} = $seq;
    }
}
close(IN);

#######################################################################################
# generate anchor seq files
#######################################################################################

print "generating anchor files in directory: $p{odir}\n";
foreach my $anchor (keys %anchors) {
    my $ofn = $p{odir}."/".$anchor.".fasta";
    # print "generating anchor file: $ofn\n";
    open(OUT, ">", $ofn) || die;
    foreach my $gene (keys %{$anchors{$anchor}}) {
	my $seq = $anchors{$anchor}->{$gene};
	print OUT ">", $gene, "\n";
	print OUT $seq, "\n";
    }
    close(OUT);
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


