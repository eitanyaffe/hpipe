#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
        print STDERR "usage: $0 <is_container file>\n";
        exit 1;
}

my $is_container = -e $ARGV[0];
while (my $line = <STDIN>) {
    chomp $line;
    next if ($line eq "" || substr($line,0,1) eq "#");
    
    if ($is_container) {
	my @f = split("=", $line);
	print $f[0], "=/links/", $f[0], "\n"; 
    } else {
	print $line, "\n";
    }
}
