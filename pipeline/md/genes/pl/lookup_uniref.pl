#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <uniref table> <usearch table> <ofn>\n";
	exit 1;
}

my $ifn_uniref = $ARGV[0];
my $ifn_genes = $ARGV[1];
my $ofn = $ARGV[2];

#######################################################################################
# read uniref table
#######################################################################################

# contig length table
my %uniref;

print STDERR "reading uniref table: $ifn_uniref\n";
open(IN, $ifn_uniref) || die $ifn_uniref;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $id = $f[$h{id}];
    $uniref{$id} = {};
    $uniref{$id}->{desc} = $f[$h{desc}];
    $uniref{$id}->{tax} = $f[$h{tax}];
    $uniref{$id}->{count} = $f[$h{count}];
}
close(IN);

print STDERR "number of uniref items: ", scalar(keys %uniref), "\n";

#######################################################################################
# traverse reads
#######################################################################################

print STDERR "going over table: $ifn_genes\n";
open(IN, $ifn_genes) || die $ifn_genes;

print STDERR "creating output: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
$header = <IN>;
chomp($header);
%h = parse_header($header);

print OUT $header, "\tprot_desc\ttax\tuniref_count\n";
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $id = $f[$h{uniref}];
    if (!defined($uniref{$id})) {
	print "Warning: id $id not found in uniref table. line=$line\n";
	next;
    }
    print OUT $line, "\t", $uniref{$id}->{desc}, "\t", $uniref{$id}->{tax}, "\t", $uniref{$id}->{count}, "\n";
}
close(IN);
close(OUT);

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
