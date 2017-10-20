#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
    print STDERR "usage: $0 <map> <self T|F> <ofn1> <ofn2>\n";
    exit 1;
}

my ($map_ifn, $self, $ofn1, $ofn2) = @ARGV;

$self = $self eq "T";

############################################################################################################
# reading map
############################################################################################################

my %map1;
my %map2;

print STDERR "reading file: $map_ifn\n";
open(IN, $map_ifn) or die;
my $header = <IN>;
my %h = parse_header($header);

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $gene1 = $f[$h{gene1}];
    my $gene2 = $f[$h{gene2}];

    next if ($self && $gene1 eq $gene2);

    my $iden = $f[$h{identity}];
    my $coverage = $f[$h{coverage}];

    if (!defined($map1{$gene1})) {
	$map1{$gene1} = {};
	$map1{$gene1}->{identity} = 0;
    }
    if ($map1{$gene1}->{identity} < $iden) {
	$map1{$gene1}->{coverage} = $coverage;
	$map1{$gene1}->{identity} = $iden;
	$map1{$gene1}->{target} = $gene2;
    }

    if (!defined($map2{$gene2})) {
	$map2{$gene2} = {};
	$map2{$gene2}->{identity} = 0;
    }
    if ($map2{$gene2}->{identity} < $iden) {
	$map2{$gene2}->{coverage} = $coverage;
	$map2{$gene2}->{identity} = $iden;
	$map2{$gene2}->{target} = $gene1;
    }
}
close(IN);

############################################################################################################
# saving ouput
############################################################################################################

print "generating table1: $ofn1\n";
open(OUT, ">", $ofn1) || die;
print OUT "gene\tidentity\tcoverage\ttarget\n";
foreach my $gene (keys %map1) {
    print OUT $gene, "\t", $map1{$gene}->{identity}, "\t", $map1{$gene}->{coverage}, "\t", $map1{$gene}->{target}, "\n";
}
close(OUT);

print "generating table2: $ofn2\n";
open(OUT, ">", $ofn2) || die;
print OUT "gene\tidentity\tcoverage\ttarget\n";
foreach my $gene (keys %map2) {
    print OUT $gene, "\t", $map2{$gene}->{identity}, "\t", $map2{$gene}->{coverage}, "\t", $map2{$gene}->{target}, "\n";
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


