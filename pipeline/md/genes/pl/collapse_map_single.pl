#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <map_ifn> <field> <max|min> <ofn>\n";
	exit 1;
}

my ($map_ifn, $max_field, $function, $ofn) = @ARGV;

###################################################################################################################
# go over map
###################################################################################################################

# from 1 to 2
our %gmap;
read_map($map_ifn, "source", "target", 1);

###################################################################################################################
# output graph
###################################################################################################################

print STDERR "generating output graph: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT "gene1\tgene2\tidentity\tcoverage\tfwd_count\tbck_count\n";

foreach my $gene1 (keys %gmap) {
    foreach my $gene2 (keys %{$gmap{$gene1}}) {
	my $coverage = $gmap{$gene1}->{$gene2}->{coverage};
	my $identity = $gmap{$gene1}->{$gene2}->{identity};
 	my $fwd_count = $gmap{$gene1}->{$gene2}->{fwd_count};
 	my $bck_count = $gmap{$gene1}->{$gene2}->{bck_count};

	print OUT "$gene1\t$gene2\t$identity\t$coverage\t$fwd_count\t$bck_count\n";
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
    my ($ifn, $field1, $field2, $fwd) = @_;
    print STDERR "reading map: $ifn\n";
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

	$gmap{$gene1} = {} if (!defined($gmap{$gene1}));
	if (!defined($gmap{$gene1}->{$gene2})) {
	    $gmap{$gene1}->{$gene2} = {};
	    $gmap{$gene1}->{$gene2}->{fwd_count} = 0;
	    $gmap{$gene1}->{$gene2}->{bck_count} = 0;
	    $gmap{$gene1}->{$gene2}->{value} = $value;
	    $gmap{$gene1}->{$gene2}->{coverage} = $coverage;
	    $gmap{$gene1}->{$gene2}->{identity} = $identity;
	}
	if ($value > $gmap{$gene1}->{$gene2}->{value}) {
	    $gmap{$gene1}->{$gene2}->{coverage} = $coverage;
	    $gmap{$gene1}->{$gene2}->{identity} = $identity;
	}

	$gmap{$gene1}->{$gene2}->{fwd_count}++ if ($fwd);
	$gmap{$gene1}->{$gene2}->{bck_count}++ if (!$fwd);
    }
    close(IN);
}

