#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig table> <input dir> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $idir = $ARGV[1];
my $odir = $ARGV[2];

my %contigs;

###############################################################################################
# setup contig bins
###############################################################################################

print STDERR "reading contig table: $ifn\n";
open(IN, $ifn) || die $ifn;
my $header = <IN>;
my %h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    $contigs{$contig} = {};
}

###############################################################################################
# traverse all reads
###############################################################################################

print STDERR "Input dir: $idir\n";
my @ifns = <$idir/*.pair>;

print STDERR "number of input files: ", scalar(@ifns), "\n";

foreach my $ifn (@ifns) {
    # print STDERR "traversing file: $ifn\n";
    print STDERR ".";
    open(IN, $ifn) || die $ifn;
    my $header = <IN>;
    my %h = parse_header($header);
    while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $contig1 = $f[$h{contig1}];
	my $coord1 = $f[$h{coord1}];
	
	my $contig2 = $f[$h{contig2}];
	my $coord2 = $f[$h{coord2}];

	next if ($contig1 ne $contig2);
	my $contig = $contig1;
	next if (!defined($contigs{$contig}));

	my ($coord_l, $coord_r) = $coord1 < $coord2 ? ($coord1,$coord2) : ($coord2,$coord1);

	my $key = $coord_l."_".$coord_r;

	if (!defined($contigs{$contig}->{$key})) {
	    $contigs{$contig}->{$key} = {};
	    $contigs{$contig}->{$key}->{coord1} = $coord_l;
	    $contigs{$contig}->{$key}->{coord2} = $coord_r;
	    $contigs{$contig}->{$key}->{count}= 0;
	}
	$contigs{$contig}->{$key}->{count}++;
    }
    close(IN);
}
print STDERR "\n";

print STDERR "generating files in: $odir\n";
foreach my $contig (keys %contigs) {
    my $ofn = $odir."/".$contig;
    open(OUT, ">", $ofn) || die $ofn;
    print OUT "coord1\tcoord2\tcount\n";
    foreach my $key (keys %{$contigs{$contig}}) {
	
    print OUT 
	$contigs{$contig}->{$key}->{coord1}, "\t", 
	$contigs{$contig}->{$key}->{coord2}, "\t", 
	$contigs{$contig}->{$key}->{count}, "\n"; 
    }
    close(OUT);
}

######################################################################################################
# Subroutines
######################################################################################################


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
