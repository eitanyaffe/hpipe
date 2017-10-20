#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig table> <key field> <value field> <fend table> <output file>\n";
	exit 1;
}

my ($icontig, $ikey, $ivalue, $ifn, $ofn) = @ARGV;

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

    my $contig = $f[$h{$ikey}];
    $contigs{$contig} = $f[$h{$ivalue}];
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
print OUT $header, "\tanchor\n";

my $anchor_count = 0;
my $other_count = 0;
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $anchor;
    if (!defined($contigs{$contig})) {
	$anchor = 0;
	$other_count++;
    } else {
	$anchor = $contigs{$contig};
	$anchor_count++;
    }
    $anchor = 0 if ($anchor <= 0);
    print OUT $line, "\t", $anchor, "\n";
}
close(IN);
close(OUT);

my $sum = $anchor_count + $other_count;
my $ap = int(1000 * $anchor_count / $sum)/10;
my $op = int(1000 * $other_count / $sum)/10;
print "number of anchor-assigned fends (including background anchor): ", $anchor_count, "($ap%)\n";
print "number of other fends: ", $other_count, "($op%)\n";

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
