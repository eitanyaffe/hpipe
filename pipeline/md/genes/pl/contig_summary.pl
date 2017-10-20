#!/usr/bin/env perl

use strict;
use File::Basename;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <genome dir> <genome filename> <ofn>\n";
	exit 1;
}

my $idir = $ARGV[0];
my $pattern = $ARGV[1];
my $ofn = $ARGV[2];

my @files = <$idir/*>;

print "generating file: $ofn\n";
open(OUT, ">", $ofn) or die $ofn;
print OUT "genome\tgenome_id\tcontig_index\tcontig\tlength\n";

my $genome = 1;

foreach my $dir (@files) {
    my $seq = "";
    my $contig = "";
    my $ifn = $dir."/".$pattern;
    my $genome_id = basename($dir);
    
    next if ($genome_id eq "united");
    open(IN, $ifn) or die $ifn;
    my $index = 1;
    while (my $line = <IN>) {
	chomp($line);
	if (substr($line, 0, 1) ne ">") {
	    $seq .= $line;
	} else {
	    print OUT  "R".$genome, "\t", $genome_id, "\t", $index++, "\t", $contig, "\t", length($seq), "\n" if ($contig ne "");
	    my @f = split(" ", substr($line,1));
	    $contig = $f[0];
	    $seq = "";
	}
    }
    close(IN);
    print OUT "R".$genome, "\t", $genome_id, "\t", $index, "\t", $contig, "\t", length($seq), "\n" if ($contig ne "");
    $genome++;
}

close(OUT);
