#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <fend table> <anchor field> <output file>\n";
	exit 1;
}

my ($ifn, $anchor_field, $ofn) = @ARGV;

#######################################################################################
# identify and bind anchors
#######################################################################################

print STDERR "traversing fends file: $ifn\n";
open(IN, $ifn) || die $ifn;
my $header = <IN>;
chomp($header);
my %h = parse_header($header);

#######################################################################################
# second pass - generate output file
#######################################################################################

print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT $header, "\n";

my $anchor_count = 0;
my $other_count = 0;
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);

    my $anchor = $f[$h{$anchor_field}];
    if ($anchor <= 0) {
	$other_count++;
    } else {
	print OUT $line, "\n";
	$anchor_count++;
    }
}
close(IN);
close(OUT);

my $sum = $anchor_count + $other_count;
my $ap = int(1000 * $anchor_count / $sum)/10;
my $op = int(1000 * $other_count / $sum)/10;
print "number of anchor fends: ", $anchor_count, "($ap%)\n";
print "number of non-anchor fends: ", $other_count, "($op%)\n";

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
