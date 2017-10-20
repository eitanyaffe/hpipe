#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print "usage: $0 <output dir> <max number of reads> <should trim> <trim read offset> <trim read length> <input fastq files>\n";
	exit 1;
}

my $odir = $ARGV[0];
my $max_reads = $ARGV[1];
my $trim = $ARGV[2] eq "T";
my $offset = $ARGV[3];
my $rlen = $ARGV[4];
shift; shift; shift; shift; shift;
my @ifns = @ARGV;

my %pfiles;
for my $ifn (@ifns) {
    my $side = (index($ifn, "R1") != -1) ? "R1" : "R2";
    my $fkey = $ifn;
    $fkey =~ s/R[12]//;
    $pfiles{$fkey}->{$side} = $ifn;
}

my $oindex = 1;

my $ofn1 = $odir."/R1_".$oindex.".fastq";
my $ofn2 = $odir."/R2_".$oindex.".fastq";

open(OUT1, ">", $ofn1) or die;
open(OUT2, ">", $ofn2) or die;
print "output file: $ofn1\n";
print "output file: $ofn2\n";

my $read_count = 0;
foreach my $fkey (keys %pfiles) {
    defined($pfiles{$fkey}->{R1}) && defined($pfiles{$fkey}->{R2}) or die "prefix doesn't have two sides: $fkey";
    my $ifn1 = $pfiles{$fkey}->{R1};
    my $ifn2 = $pfiles{$fkey}->{R2};
    open(IN1, $ifn1) or die;
    open(IN2, $ifn2) or die;

    my $iindex = 0;
    while () {
	my $line1 = <IN1>;
	my $line2 = <IN2>;
	
	last if (!defined($line1) || !defined($line2));

	chomp($line1);
	chomp($line2);

	$line1 = substr($line1, $offset, $rlen) if ($trim && $iindex % 2 == 1);
	$line2 = substr($line2, $offset, $rlen) if ($trim && $iindex % 2 == 1);
 	$iindex++;

	print OUT1 $line1, "\n";
	print OUT2 $line2, "\n";
	
	if ($iindex % 4 == 0) {
	    $read_count++;
	
	    if ($read_count >= $max_reads) {
		$read_count = 0;
		close(OUT1);
		close(OUT2);
		$oindex++;
		
		$ofn1 = $odir."/R1_".$oindex.".fastq";
		$ofn2 = $odir."/R2_".$oindex.".fastq";
		
		open(OUT1, ">", $ofn1) or die;
		open(OUT2, ">", $ofn2) or die;
		print "output file: $ofn1\n";
		print "output file: $ofn2\n";
	    }
	}
    }
}

close(OUT1);
close(OUT2);
