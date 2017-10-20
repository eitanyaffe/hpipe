#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <binary> <output file R1> <output file R2> <complexity ofn> <stats ofn> <input dir1 dir2 ...>\n";
	exit 1;
}

my $bin = $ARGV[0];
my $ofn1 = $ARGV[1];
my $ofn2 = $ARGV[2];
my $ofn_complex = $ARGV[3];
my $ofn_stats = $ARGV[4];
shift; shift; shift; shift; shift;
my @idirs = @ARGV;

my $command = sprintf("%s -ofn1 %s -ofn2 %s -mfn %s -sfn %s", $bin, $ofn1, $ofn2, $ofn_complex, $ofn_stats);
print STDERR "Input dirs: ", join(" ", @idirs), "\n";
foreach my $idir (@idirs) {
    print STDERR "looking for files in dir: $idir\n";
    my @ifns = <$idir/*.fastq>;

    @ifns > 0 or die "no files in $idir";

    my %pfiles;
    for my $ifn (@ifns) {
	my $side = (index($ifn, "R1") != -1) ? "R1" : "R2";
	my $fkey = $ifn;
	$fkey =~ s/R[12]//;
	$pfiles{$fkey}->{$side} = $ifn;
    }

    foreach my $fkey (keys %pfiles) {
	defined($pfiles{$fkey}->{R1}) && defined($pfiles{$fkey}->{R2}) or die "prefix doesn't have two sides: $fkey";
	my $ifn1 = $pfiles{$fkey}->{R1};
	my $ifn2 = $pfiles{$fkey}->{R2};
	$command .= " -ifn1 ".$ifn1." -ifn2 ".$ifn2;
    }
}
print "command: $command\n";
system($command) == 0 or die;
