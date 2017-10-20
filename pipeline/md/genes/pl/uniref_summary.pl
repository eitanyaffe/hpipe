#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

print "gene\taa_length\n";

my $seq = "";
my $contig = "";
while (my $line = <STDIN>) {
    chomp($line);
    if (substr($line, 0, 1) ne ">") {
	$seq .= $line;
    } else {
	print $contig, "\t", length($seq), "\n" if ($contig ne "");
	my @f = split(" ", substr($line,1));
	$contig = $f[0];
	$seq = "";
    }
    }
print $contig, "\t", length($seq), "\n" if ($contig ne "");



