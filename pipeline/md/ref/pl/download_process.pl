#!/usr/bin/env perl

use strict;
use File::Basename;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <genome dir> <genome filename> <truncate T|F> <lines> <odir>\n";
	exit 1;
}

my $idir = $ARGV[0];
my $pattern = $ARGV[1];
my $truncate = $ARGV[2];
my $truncate_lines = $ARGV[3];
my $odir = $ARGV[4];

my @files = <$idir/*>;

($truncate eq "T") ? print "using only the first $truncate_lines lines of each genome\n" : print "using all downloaded genomes\n";
foreach my $dir (@files) {
    my $ifn = $dir."/".$pattern;
    my $genome_id = basename($dir);
    my $ofn = $odir."/".$genome_id;
    if ($truncate eq "T") {
 	system("head -n $truncate_lines $ifn > $ofn");
    } else {
 	system("cp $ifn $ofn");
    }
}
