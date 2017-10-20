#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

while (my $line = <STDIN>) {
    if (substr($line, 0, 1) eq ">") {
	chomp($line);
	my @f = split(" ", substr($line,1));
	my $contig = $f[0];
	$contig = substr($contig, 0, index($contig, "_"));
	print ">c", $contig, "\n";
    } else {
	print $line;
    }
}



