#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
    print STDERR "usage: $0 <fend table> <contig cluster table> <ofn>\n";
    exit 1;
}

my $ifn = $ARGV[0];
my $cluster_table = $ARGV[1];
my $ofn = $ARGV[2];

#######################################################################################
# read cluster table
#######################################################################################

my %contigs;

print STDERR "reading cluster table: $cluster_table\n";
open(IN, $cluster_table) || die $cluster_table;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $cluster = $f[$h{cluster}];
    $contigs{$contig} = $cluster;
}
close(IN);

#######################################################################################
# traverse fends
#######################################################################################

print STDERR "traversing fends file: $ifn\n";
open(IN, $ifn) || die $ifn;
$header = <IN>;
chomp($header);
%h = parse_header($header);

print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;

print OUT $header, "\tcluster\n";

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    next if (!defined($contigs{$contig}));
    my $cluster = $contigs{$contig};
    print OUT $line, "\t", $cluster, "\n";
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
