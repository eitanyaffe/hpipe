#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ofn = $ARGV[1];

print STDERR "traversing file: $ifn\n";
open(IN, $ifn) || die $ifn;

print STDERR "creating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "id\tdesc\tcount\ttax\trepid\n";
while (my $line = <IN>) {
    chomp($line);
    if (substr($line, 0, 1) eq ">") {
	my @f = split("n=|Tax=|RepID=", substr($line,1));
	my @f1 = split(/\s+/, $f[0]);
	my $id = $f1[0];
	shift @f1;
	my $desc = join(" ", @f1);

	my $count = $f[1];
	my $tax = $f[2];
	my $repid = $f[3];

	print OUT "$id\t$desc\t$count\t$tax\t$repid\n";
    }
}
close(IN);
