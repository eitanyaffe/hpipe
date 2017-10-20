#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <binsize> <ofn bins> <ofn table>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $binsize = $ARGV[1];
my $ofn_bins = $ARGV[2];
my $ofn_table = $ARGV[3];

print STDERR "traversing file: $ifn\n";
open(IN, $ifn) || die $ifn;

print STDERR "creating file: $ofn_table\n";
open(OUT_TABLE, ">", $ofn_table) || die $ofn_table;
print OUT_TABLE "contig\tgc\n";

print STDERR "creating file: $ofn_bins\n";
open(OUT_BINS, ">", $ofn_bins) || die $ofn_bins;
print OUT_BINS "contig\tcoord\tgc\n";

my $contig = "";
my $seq = "";
while (my $line = <IN>) {
    chomp($line);
    if (substr($line, 0, 1) eq ">") {

	if ($contig ne "") {
	    print OUT_TABLE $contig, "\t", compute_gc($seq), "\n";
	    my $nbins = int(length($seq)/$binsize);
	    for (my $i=0; $i<$nbins; $i++) {
		my $sstr = substr($seq,$i*$binsize,$binsize);
		print OUT_BINS $contig, "\t", $i*$binsize, "\t", compute_gc($sstr), "\n";

	    }
	}
	my @f = split(" ", substr($line,1));
	$contig = $f[0];
	$seq = "";
    } else {
	$seq .= $line;
    }
}
print OUT_TABLE $contig, "\t", compute_gc($seq), "\n";
my $nbins = int(length($seq)/$binsize);
for (my $i=0; $i<$nbins; $i++) {
    my $sstr = substr($seq,$i*$binsize,$binsize);
    print OUT_BINS $contig, "\t", $i*$binsize, "\t", length($sstr), "\t", compute_gc($sstr), "\n";
}

close(IN);
close(OUT_TABLE);
close(OUT_BINS);

sub compute_gc {
    my ($seq) = @_;
    my $gc = 0;
    my $at = 0;
    for(my $i = 0; $i < length($seq); $i++) {
	my $c = uc(substr($seq, $i, 1));
	if($c eq "C" || $c eq "G") {
	    $gc++;
	} elsif($c eq "A" || $c eq "T") {
	    $at++;
	}
    }
    return $gc/($gc+$at);
}
