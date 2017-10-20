#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig table> <min relative abundance> <fends table> <output file>\n";
	exit 1;
}

my $icontig = $ARGV[0];
my $min_abundance = $ARGV[1];
my $ifn = $ARGV[2];
my $ofn = $ARGV[3];

print "minimal relative abundance: ", $min_abundance, "\n";

#######################################################################################
# read contig table
#######################################################################################

# contig length table
my %contigs;

print STDERR "reading contig table: $icontig\n";
open(IN, $icontig) || die $icontig;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $coverage = $f[$h{"abundance.enrichment"}];

    $contigs{$contig} = $coverage;
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
print OUT $header, "\tabundance\n";

my $skip = 0;
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $coverage = $contigs{$contig};


    if ($coverage < $min_abundance) {
	$skip++;
	next
    }

    print OUT $line, "\t", $coverage, "\n";
}
close(IN);
close(OUT);

print "number of skipped fends: ", $skip, "\n";

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
