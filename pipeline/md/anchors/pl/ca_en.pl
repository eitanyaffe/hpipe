#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <input fends> <input mat> <output file>\n";
	exit 1;
}

my $fends_ifn = $ARGV[0];
my $mat_ifn = $ARGV[1];
my $ofn = $ARGV[2];

#############################################################################################
# fend file
#############################################################################################

my %fends;
my %contigs;
my %anchors;

open(IN, $fends_ifn) || die;
print STDERR "Reading file $fends_ifn into hash...\n";
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $fend = $f[$h{fend}];
    my $anchor = $f[$h{anchor}];
    my $contig = $f[$h{contig}];
    my $coord = $f[$h{coord}];

    !defined($fends{$fend}) or die "non-unique fend";
    $fends{$fend} = {};
    $fends{$fend}->{anchor} = $anchor;
    $fends{$fend}->{contig} = $contig;

    # save all anchors
    $anchors{$anchor} = 0;

    if (!defined($contigs{$contig})) {
	$contigs{$contig} = {};
	$contigs{$contig}->{anchor} = $anchor;
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

our %contig_table;

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

    defined($contigs{$contig1}) or die;
    defined($contigs{$contig2}) or die;

    add_contig($contig1, $anchor2) if ($anchor2 != 0);
    add_contig($contig2, $anchor1) if ($anchor1 != 0);
    $ucount++;
}

close(IN);
print STDERR "Number of inter contig unique contacts: $ucount\n";

print STDERR "Writing output file: $ofn\n";

open(OUT, ">", $ofn) || die;
print OUT "contig\tcontig_anchor\tanchor\tcount\tmarginal_count\n";
foreach my $contig (keys %contigs) {
foreach my $anchor (sort {$a <=> $b} keys %anchors) {
    next if ($anchor == 0);
    my $contig_anchor = $contigs{$contig}->{anchor};

    # skip if there are no contacts between anchor and all of contig
    next if (!defined($contig_table{$contig}->{$anchor}));

    my $count = $contig_table{$contig}->{$anchor};
    my $marginal = $contig_table{$contig}->{marginal};
    print OUT "$contig\t$contig_anchor\t$anchor\t$count\t$marginal\n";
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
    my ($contig, $anchor) = @_;
    if (!defined($contig_table{$contig})) {
	$contig_table{$contig} = {};
	$contig_table{$contig}->{marginal} = 0;
    }
    $contig_table{$contig}->{marginal}++;
    
    $contig_table{$contig}->{$anchor} = 0 if (!defined($contig_table{$contig}->{$anchor}));
    $contig_table{$contig}->{$anchor}++;
}
