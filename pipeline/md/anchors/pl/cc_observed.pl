#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <input fends> <cluster field in fend table> <input mat> <bins per contig> <output file>\n";
	exit 1;
}

my $fends_ifn = $ARGV[0];
my $cluster_field = $ARGV[1];
my $mat_ifn = $ARGV[2];
my $bins_per_contig = $ARGV[3];
my $ofn = $ARGV[4];

#############################################################################################
# fend file
#############################################################################################

our %fends;
my %contigs;

open(IN, $fends_ifn) || die;
print STDERR "Reading file $fends_ifn into hash...\n";
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $fend = $f[$h{fend}];
    my $contig = $f[$h{contig_bin}];
    
    !defined($fends{$fend}) or die "non-unique fend";
    $fends{$fend} = {};
    $fends{$fend}->{anchor} = $f[$h{$cluster_field}];
    $fends{$fend}->{contig} = $contig;

    if (!defined($contigs{$contig})) {
	$contigs{$contig} = {};
	$contigs{$contig}->{fend_count} = 0;
    }
    
    $fends{$fend}->{contig_index} = $contigs{$contig}->{fend_count};
    $contigs{$contig}->{fend_count}++;

}
close(IN);

#############################################################################################
# mat file
#############################################################################################

my $appr_lines = apprx_lines($mat_ifn);
print STDERR "traversing file $mat_ifn, with about ".int($appr_lines/1000000)."M lines\n";

my %contig_table;

my $lcount = 0;
open(IN, $mat_ifn) || die;
$header = <IN>;
%h = parse_header($header);
my $ucount = 0;
while (my $line = <IN>) 
{
    $lcount++;
    print STDERR "line: $lcount\n" if ($lcount % 1000000 == 0);

    chomp $line;
    my @f = split("\t", $line);
    my $fend1 = $f[$h{fend1}];
    my $fend2 = $f[$h{fend2}];
    my $count = $f[$h{count}];
    
    next if (!defined ($fends{$fend1}) or !defined ($fends{$fend2}));
    
    my $contig1 = $fends{$fend1}->{contig};
    my $contig2 = $fends{$fend2}->{contig};
    next if ($contig1 eq $contig2);
    
    my $anchor1 = $fends{$fend1}->{anchor};
    my $anchor2 = $fends{$fend2}->{anchor};
	
    defined($contigs{$contig1}) or die $contig1;
    defined($contigs{$contig2}) or die;

    $contigs{$contig1}->{nbins} = min($bins_per_contig, $contigs{$contig1}->{fend_count}) if (!defined($contigs{$contig1}->{nbins}));
    $contigs{$contig2}->{nbins} = min($bins_per_contig, $contigs{$contig2}->{fend_count}) if (!defined($contigs{$contig2}->{nbins}));

    if ($anchor1 != 0) {
	my $nbins2 = $contigs{$contig2}->{nbins};
	my $contig_bin2 = int(($fends{$fend2}->{contig_index} / $contigs{$contig2}->{fend_count}) * $nbins2);
	add_contig($contig2, $contig_bin2, $nbins2, $anchor1);
    }
    if ($anchor2 != 0) {
	my $nbins1 = $contigs{$contig1}->{nbins};
	my $contig_bin1 = int(($fends{$fend1}->{contig_index} / $contigs{$contig1}->{fend_count}) * $nbins1);
	add_contig($contig1, $contig_bin1, $nbins1, $anchor2);
    }

    $ucount++;
}
close(IN);
print STDERR "Number of inter contig unique contacts: $ucount\n";

print STDERR "Writing output file: $ofn\n";

open(OUT, ">", $ofn) || die;
print OUT "contig\tanchor\tcontig_total_count\tcontig_median_count\n";
foreach my $contig (keys %contig_table) {
foreach my $anchor (sort {$a <=> $b} keys %{$contig_table{$contig}}) {

    my $nbins = $contigs{$contig}->{nbins};
    my $contig_total_count = $contig_table{$contig}->{$anchor}->{total};
    my $contig_median_count = median(@{$contig_table{$contig}->{$anchor}->{bins}}) * $nbins;
    print OUT "$contig\t$anchor\t$contig_total_count\t$contig_median_count\n";
} }
close(OUT);

######################################################################################################
# Subroutines
######################################################################################################

sub median {
  (sort { $a <=> $b } @_ )[ int( $#_/2 ) ];
}

sub min
{
    my ($a, $b) = @_;
    return $a < $b ? $a : $b;
}

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

sub add_contig
{
    my ($contig, $bin, $nbins, $anchor) = @_;
    $contig_table{$contig} = {} if (!defined($contig_table{$contig}));
    ($bin >= 0 and $bin < $nbins) or die;

    if (!defined($contig_table{$contig}->{$anchor})) {
	$contig_table{$contig}->{$anchor} = {};
	$contig_table{$contig}->{$anchor}->{bins} = [(0) x $nbins];
	$contig_table{$contig}->{$anchor}->{total} = 0;
    }
    $contig_table{$contig}->{$anchor}->{bins}[$bin]++;
    $contig_table{$contig}->{$anchor}->{total}++;
}
