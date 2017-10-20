#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ga_ifn> <cgene_table_ifn> <map_ifn> <ofn>\n";
	exit 1;
}

my ($ga_ifn, $cgene_table_ifn, $map_ifn, $function, $ofn) = @ARGV;

###################################################################################################################
# read ga
###################################################################################################################

my %anchors;

print STDERR "reading ga: $ga_ifn\n";
open(IN, $ga_ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $cgene = $f[$h{cgene}];
    my $anchor = $f[$h{anchor}];

    $anchors{$anchor} = {} if (!defined($anchors{$anchor}));
    $anchors{$anchor}->{$cgene} = 0;
}
close(IN);


###################################################################################################################
# read gene table into memory
###################################################################################################################

my %cgenes;
print STDERR "reading gene table: $cgene_table_ifn\n";
open(IN, $cgene_table_ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $cgene = $f[$h{cgene}];
    if (!defined($cgenes{$cgene}) {
	$cgenes{$cgene} = {};
    $genes{$gene}->{cgene} = $cgene;
}
close(IN);

###################################################################################################################
# go over map
###################################################################################################################

# from 1 to 2
our %gmap;
read_map($map_ifn, "source", "target");

###################################################################################################################
# output graph
###################################################################################################################

print STDERR "generating output graph: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT "cgene1\tcgene2\tidentity\tcoverage\tcount\n";

foreach my $cgene1 (keys %gmap) {
    foreach my $cgene2 (keys %{$gmap{$cgene1}}) {
	next if ($cgene1 eq $cgene2);

	my ($identity_v, $coverage_v, $count, $first) = (0,0,0,1);
	foreach my $key (keys %{$gmap{$cgene1}->{$cgene2}}) {
	    $count++;
	    my $tidentity = $gmap{$cgene1}->{$cgene2}->{$key}->{identity};
	    my $tcoverage = $gmap{$cgene1}->{$cgene2}->{$key}->{coverage};
	    if ($function eq "min") {
		$identity_v = $first ? $tidentity : ($tidentity < $identity_v ? $tidentity : $identity_v);
		$coverage_v = $first ? $tcoverage : ($tcoverage < $coverage_v ? $tcoverage : $coverage_v);
	    } elsif ($function eq "max") {
		$identity_v = $first ? $tidentity : ($tidentity > $identity_v ? $tidentity : $identity_v);
		$coverage_v = $first ? $tcoverage : ($tcoverage > $coverage_v ? $tcoverage : $coverage_v);
	    }  elsif ($function eq "mean") {
		$identity_v += $tidentity;
		$coverage_v += $tcoverage;
	    } else {
		die "unknown collapse function: $function";
	    }

	}
	if ($function eq "mean") {
	    $identity_v = $identity_v / $count;
	    $coverage_v = $coverage_v / $count;
	}
	print OUT "$cgene1\t$cgene2\t$identity_v\t$coverage_v\t$count\n";
    }
}
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

sub read_map
{
    my ($ifn, $field1, $field2) = @_;
    print STDERR "reading map: $ifn\n";
    open(IN, $ifn) || die;
    my $header = <IN>;
    my %h = parse_header($header);
    while (my $line = <IN>) {
	chomp($line);
	my @f = split("\t", $line);
	my $identity = $f[$h{identity}];
	my $coverage = $f[$h{coverage}];
	my $evalue = $f[$h{evalue}];
	my $gene1 = $f[$h{$field1}];
	my $gene2 = $f[$h{$field2}];

	defined($genes{$gene1}) or die;
	defined($genes{$gene2}) or die;
	my $cgene1 = $genes{$gene1};
	my $cgene2 = $genes{$gene2};
	$gmap{$cgene1} = {} if (!defined($gmap{$cgene1}));
	$gmap{$cgene1}->{$cgene2} = {} if (!defined($gmap{$cgene1}->{$cgene2}));

	my $key = $gene1."_".$gene2;
	next if (defined($gmap{$cgene1}->{$cgene2}->{$key}) && $gmap{$cgene1}->{$cgene2}->{$key}->{evalue} < $evalue);

	$gmap{$cgene1}->{$cgene2}->{$key} = {};
	$gmap{$cgene1}->{$cgene2}->{$key}->{coverage} = $coverage;
	$gmap{$cgene1}->{$cgene2}->{$key}->{identity} = $identity;
	$gmap{$cgene1}->{$cgene2}->{$key}->{evalue} = $evalue;
    }
    close(IN);
}

