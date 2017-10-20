#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <field> <max|min> <ofn> <map_ifn1, map_ifn...>\n";
	exit 1;
}

my ($max_field, $function, $ofn) = @ARGV[0..2];
shift; shift; shift;
my @ifns = @ARGV;

# print "merging maps: ", join(" \n", @ifns), "\n";

###################################################################################################################
# go over map
###################################################################################################################

# from 1 to 2
our %gmap;
for my $ifn (@ifns) {
    read_map($ifn, "gene", "target", 1);
}

###################################################################################################################
# output graph
###################################################################################################################

# print STDERR "generating output graph: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "gene\ttarget\tidentity\tcoverage\tfwd_count\tbck_count\n";

foreach my $gene (keys %gmap) {
    my $coverage = $gmap{$gene}->{coverage};
    my $identity = $gmap{$gene}->{identity};
    my $fwd_count = $gmap{$gene}->{fwd_count};
    my $bck_count = $gmap{$gene}->{bck_count};
    my $best = $gmap{$gene}->{best_match};

    print OUT "$gene\t$best\t$identity\t$coverage\t$fwd_count\t$bck_count\n";
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
    my ($ifn, $field1, $field2, $fwd) = @_;
    # print STDERR "reading map: $ifn\n";
    open(IN, $ifn) || die;
    my $header = <IN>;
    my %h = parse_header($header);
    while (my $line = <IN>) {
	chomp($line);
	my @f = split("\t", $line);
	my $identity = $f[$h{identity}];
	my $coverage = $f[$h{coverage}];
	my $gene1 = $f[$h{$field1}];
	my $gene2 = $f[$h{$field2}];
	my $value = $f[$h{$max_field}];

	if ($function eq "min") {
	    $value *= -1;
	}

	if (!defined($gmap{$gene1})) {
	    $gmap{$gene1} = {};
	    $gmap{$gene1}->{fwd_count} = 0;
	    $gmap{$gene1}->{bck_count} = 0;
	    $gmap{$gene1}->{value} = $value;
	    $gmap{$gene1}->{coverage} = $coverage;
	    $gmap{$gene1}->{identity} = $identity;
	    $gmap{$gene1}->{best_match} = $gene2;
	}
	if ($value > $gmap{$gene1}->{value}) {
	    $gmap{$gene1}->{coverage} = $coverage;
	    $gmap{$gene1}->{identity} = $identity;
	    $gmap{$gene1}->{best_match} = $gene2;
	}

	$gmap{$gene1}->{fwd_count}++ if ($fwd);
	$gmap{$gene1}->{bck_count}++ if (!$fwd);
    }
    close(IN);
}

