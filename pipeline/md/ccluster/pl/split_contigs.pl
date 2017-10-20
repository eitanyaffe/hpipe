#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
    print "usage: $0 <contig table> <contig file> <n files> <odir>\n";
    exit (0);
}

my $itable = $ARGV[0];
my $ifn = $ARGV[1];
my $ncores = $ARGV[2];
my $odir = $ARGV[3];

system("rm -rf $odir") == 0 or die;
system("mkdir -p $odir") == 0 or die;

my $tot = 0;
open(IN, $itable) or die $itable;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    $tot += $f[$h{length}];
}
close(IN);

my $len_step =  1 + int($tot / $ncores);
my $len = 0;
my $ifile = 0;

print "total contig length: $tot\n";
print "number of split files: $ncores\n";
print "length of contigs per file: $len_step\n";
print "reading contig file: $ifn\n";
print "generating output in: $odir\n";

open(IN, $ifn) or die $ifn;
while (my $line = <IN>) {
    chomp($line);
    $len += length($line) if (substr($line, 0, 1) ne ">");
    if (($ifile == 0) || ( (substr($line, 0, 1) eq ">") && ($len > $len_step))) {
	close(OUT) if ($ifile > 0);
	$ifile++;
	my $ofn = $odir."/".$ifile;
	open(OUT, ">", $ofn);
	$len = 0;
    }
    print OUT $line, "\n";
}
close(OUT);
close(IN);

$ifile <= $ncores  or die;


sub parse_header
{
	my ($header) = @_;
	chomp($header);
	my @f = split("\t", $header);
	my %result;
	for (my $i = 0; $i <= $#f; $i++) {
		$result{$f[$i]} = $i;
	}
	return %result;
}
