#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print "usage: $0 <input dir> <input prefix pattern> <input suffix pattern> <output title> <output file>\n";
	exit 1;
}

my $idir = $ARGV[0];
my $pre_pattern = $ARGV[1];
my $suf_pattern = $ARGV[2];
my $title = $ARGV[3];
my $ofn = $ARGV[4];

print "computing stats in $title directory: $idir\n";
print "prefix pattern: $pre_pattern\n";
print "suffix pattern: $suf_pattern\n";

my $p1 = $pre_pattern."R1".$suf_pattern;
my $p2 = $pre_pattern."R2".$suf_pattern;

my @ifns1 = <$idir/$p1>;
my @ifns2 = <$idir/$p2>;

scalar(@ifns1) > 0 && scalar(@ifns2) > 0 or die "no files found in $idir";

print "side1 files: ", join(",", @ifns1), "\n";
print "side2 files: ", join(",", @ifns2), "\n";

my ($read_count1, $bp_count1) = parse_files(\@ifns1);
my ($read_count2, $bp_count2) = parse_files(\@ifns2);

open(OUT, ">", $ofn) or die;
print OUT $title, "\t", "R1", "\t", $read_count1, "\t", $bp_count1, "\n";
print OUT $title, "\t", "R2", "\t", $read_count2, "\t", $bp_count2, "\n";
close(OUT);

sub parse_files
{
    my ($ref) = @_;
    my @ifns = @{ $ref };
    my ($read_count, $bp_count) = (0,0);
    foreach my $ifn (@ifns) {
	next if (index($ifn, "~") != -1);
	open(IN, $ifn) or die;
	my $l_count = 0;
	while (my $line = <IN>) {
	    if ($l_count % 4 == 1) {
		chomp($line);
		$bp_count += length($line);
		$read_count++;
	    }
	    $l_count++;
	}
    }
    return ($read_count, $bp_count);
}
