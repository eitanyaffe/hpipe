#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

my $ofn = $ARGV[0];

my $qcontig = "";
my $qstrand = "+";
open(OUT, ">", $ofn);
print OUT "contig1\tstart1\tend1\tcontig2\tstart2\tend2\tstrand\tlength\n";

while (my $line = <STDIN>) {
    chomp($line);
    my @f = split(/\s+/, $line);
    if (substr($line, 0, 1) eq ">") {
	$qcontig = $f[1];
	
	if (scalar(@f) == 2) {
	    $qstrand = "+";
	} else {
	    (scalar(@f) == 3 and $f[2] eq "Reverse") or die;
	    $qstrand = "-";
	}

    } else {
	$qcontig ne "" or die;
	my $length = $f[4];

	my $contig1 = $qcontig;
	my $start1 = ($qstrand eq "+") ? $f[3] : ($f[3] - $length + 1);
	my $end1 =$start1 + $length;

	my $contig2 = $f[1];
	my $start2 = $f[2];
	my $end2 =$start2 + $length;

	next if (($contig1 eq $contig2) &&  $qstrand eq "+");

	if (($contig1 cmp $contig2) == 1) {
	    print OUT "$contig1\t$start1\t$end1\t$contig2\t$start2\t$end2\t";
	} else {
	    print OUT "$contig2\t$start2\t$end2\t$contig1\t$start1\t$end1\t";
	}
	print OUT "$qstrand\t$length\n";
    }
}
close(OUT);
