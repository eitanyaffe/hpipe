#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig table> <nt_table> <link_table> <variant> <percentage cutoff> <count cutoff> <ofn>\n";
	exit 1;
}

my $ifn_contig = $ARGV[0];
my $ifn_nt = $ARGV[1];
my $ifn_link = $ARGV[2];
my $min_perc = $ARGV[3];
my $min_count = $ARGV[4];
my $ofn = $ARGV[5];

###############################################################################################
# setup contig bins
###############################################################################################

my %contigs;
print STDERR "reading contig table: $ifn_contig\n";
open(IN, $ifn_contig) || die $ifn_contig;
my $header = <IN>;
my %h = parse_header($header);

my @types = ("snp", "insert", "delete", "dangle");

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $length = $f[$h{length}];
    $contigs{$contig} = {};
    $contigs{$contig}->{length} = $length;
    $contigs{$contig}->{counts} = {};
    for my $type (@types) {
	$contigs{$contig}->{counts}->{$type} = 0;
    }
}

###############################################################################################
# traverse all reads
###############################################################################################

foreach my $ifn (($ifn_nt, $ifn_link)) {
    open(IN, $ifn) || die $ifn;
    my $header = <IN>;
    my %h = parse_header($header);
    while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $contig = $f[$h{contig}];
	my $coord = $f[$h{coord}];
	my $type = $f[$h{type}];
	my $count = $f[$h{count}];
	my $perc = $f[$h{percent}];
	next if (!defined($contigs{$contig}));
	next if ($count < $min_count || $perc < $min_perc);
	next if ($type eq "REF");

	$type = "dangle" if ($type eq "dangle_left" || $type eq "dangle_right");
	$contigs{$contig}->{counts}->{$type}++;
    }
    close(IN);
}

open(OUT, ">", $ofn) || die $ofn;
print STDERR "generating file: $ofn\n";
print OUT "contig\tlength";
for my $type (@types) {
    print OUT "\t", $type, "\t", $type, "_per_kb";
}
print OUT "\ttotal\ttotal_per_kb\n";
foreach my $contig (sort keys %contigs) {
    my $length = $contigs{$contig}->{length};
    print OUT $contig, "\t", $length;
    my $total = 0;
    for my $type (@types) {
	my $count = $contigs{$contig}->{counts}->{$type};
	$total += $count;
	my $count_per_kb = floor(100 * (1000 * $count / $length)) / 100;
	print OUT "\t", $count, "\t", $count_per_kb;
    }
    my $total_per_kb = floor(100 * (1000 * $total / $length)) / 100;
    print OUT "\t", $total, "\t", $total_per_kb, "\n";
}
close(OUT);

######################################################################################################
# Subroutines
######################################################################################################


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
