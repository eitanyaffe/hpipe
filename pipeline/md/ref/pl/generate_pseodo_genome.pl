#!/usr/bin/env perl

use strict;
use File::Basename;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn contig table> <ifn fasta> <ofn genome table> < ofn fasta>\n";
	exit 1;
}

my $ifn_contigs = $ARGV[0];
my $ifn_fasta = $ARGV[1];
my $ofn_table = $ARGV[2];
my $ofn_fasta = $ARGV[3];

####################################################################################################
# contig table
####################################################################################################

my %genomes;
my %contigs;

print "reading table: ", $ifn_contigs, "\n";
open(IN, $ifn_contigs) or die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $genome = $f[$h{genome}];
    my $contig = $f[$h{contig}];
    my $index = $f[$h{contig_index}];

    # contigs
    $contigs{$contig} = {};
    $contigs{$contig}->{genome} = $genome;
    $contigs{$contig}->{index} = $index;

    # genomes
    $genomes{$genome} = {} if (!defined($genomes{$genome}));
    $genomes{$genome}->{$index} = $contig;
}

####################################################################################################
# contig fasta
####################################################################################################

my $seq = "";
my $contig = "";
print "reading fasta: ", $ifn_fasta, "\n";
open(IN, $ifn_fasta) or die;
while (my $line = <IN>) {
    chomp($line);
    if (substr($line, 0, 1) ne ">") {
	$line =~ s/N//g;
	$seq .= $line;
    } else {
 	if ($seq ne "") {
	    $contigs{$contig}->{seq} = $seq;
	    $contigs{$contig}->{length} = length($seq);
	}
	my @f = split(" ", substr($line,1));
	$contig = $f[0];
	$seq = "";
    }
}
close(IN);
if ($seq ne "") {
    $contigs{$contig}->{seq} = $seq;
    $contigs{$contig}->{length} = length($seq);
}

####################################################################################################
# create output
####################################################################################################

print "generating fasta file: $ofn_fasta\n";
open(OUT_FASTA, ">", $ofn_fasta) or die;

print "generating table file: $ofn_table\n";
open(OUT_TABLE, ">", $ofn_table) or die;

print OUT_TABLE "genome\tlength\n";

foreach my $genome (keys %genomes) {
    my $length = 0;
    my $seq = "";
    foreach my $index (sort {$a <=> $b} keys %{$genomes{$genome}}) {
	my $contig = $genomes{$genome}->{$index};
	defined($contigs{$contig}) or die $contig;
	$length += $contigs{$contig}->{length};
	$seq .= $contigs{$contig}->{seq};
    }

    print OUT_TABLE "$genome\t$length\n";

    print OUT_FASTA ">$genome $length\n";
    print OUT_FASTA "$seq\n";
}
close(OUT_FASTA);
close(OUT_TABLE);

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

