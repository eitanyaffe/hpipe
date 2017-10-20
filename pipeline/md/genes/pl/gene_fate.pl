#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util;

if ($#ARGV == -1) {
    print STDERR "usage: $0 <cgene_table_ifn1> <cgene_table_ifn2> <unique_ifn1> <unique_ifn2> <identity_threshold> <coverage_threshold> <ofn1> <ofn2>\n";
	exit 1;
}

my ($cgene_table_ifn1, $cgene_table_ifn2, $umap_ifn1, $umap_ifn2, $identity_threshold, $coverage_threshold, $ofn1, $ofn2) = @ARGV;

print "identity threshold: $identity_threshold\n";
print "coverage threshold: $coverage_threshold\n";

print "ofn1: $ofn1\n";
print "ofn2: $ofn2\n";

my $self = $cgene_table_ifn1 eq $cgene_table_ifn2;

process($cgene_table_ifn1, $umap_ifn1, $identity_threshold, $coverage_threshold, "lost", $ofn1);
process($cgene_table_ifn2, $umap_ifn2, $identity_threshold, $coverage_threshold, "gained", $ofn2);

sub process
{
    my ($cgene_table_ifn, $umap_ifn, $identity_threshold, $coverage_threshold, $default, $ofn) = @_;

    my %cgenes;
    print STDERR "reading gene table: $cgene_table_ifn\n";
    open(IN, $cgene_table_ifn) || die;
    my $header = <IN>;
    my %h = parse_header($header);
    while (my $line = <IN>) {
	chomp($line);
	my @f = split("\t", $line);
	my $cgene = $f[$h{cgene}];

	$cgenes{$cgene} = {};
	$cgenes{$cgene}->{fate} = $default;
	$cgenes{$cgene}->{identity} = 0;
    }
    close(IN);

    print STDERR "reading unique map: $umap_ifn\n";
    open(IN, $umap_ifn) || die;
    $header = <IN>;
    %h = parse_header($header);
    while (my $line = <IN>) {
	chomp($line);
	my @f = split("\t", $line);
	my $identity = $f[$h{identity}];
	my $coverage = $f[$h{coverage}];
	my $cgene = $f[$h{cgene}];
	defined($cgenes{$cgene}) or die;
	if ($identity < $identity_threshold || $coverage < $coverage_threshold) {
 	    $cgenes{$cgene}->{fate} = "weak";
 	    $cgenes{$cgene}->{identity} = $identity;
	} else {
	    $cgenes{$cgene}->{fate} = "found";
 	    $cgenes{$cgene}->{identity} = $identity;
	}
    }
    close(IN);

    print STDERR "generating output graph: $ofn\n";
    open(OUT, ">", $ofn) || die;
    print OUT "cgene\tfate\tidentity\n";
    foreach my $cgene (keys %cgenes) {
	print OUT $cgene, "\t", $cgenes{$cgene}->{fate}, "\t", $cgenes{$cgene}->{identity}, "\n";
    }
    close(OUT);
}

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
    my ($ifn) = @_;
}

