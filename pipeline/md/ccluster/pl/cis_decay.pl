#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig table> <input dir> <bin start> <bin end> <bin step> <gap> <max reads> <ofn decay> <ofn summary>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $idir = $ARGV[1];
my $bstart = $ARGV[2];
my $bend = $ARGV[3];
my $bstep = $ARGV[4];
my $gap = $ARGV[5];
my $mreads = $ARGV[6];
my $ofn_decay = $ARGV[7];
my $ofn_summary = $ARGV[8];

print "bin log start: $bstart\n";
print "bin log end: $bend\n";
print "bin log step: $bstep\n";
print "edge gap: $gap\n";
my %contigs;
my %bins;

###############################################################################################
# bins
###############################################################################################

my $nbins = int(($bend - $bstart) / $bstep);
for (my $i = 0; $i <= $nbins; $i++)
{
    my $log_bin = $bstart + $bstep * $i;
    my $bin = int(10**$log_bin);
    $bins{$bin} = {};
    $bins{$bin}->{log_start} = $log_bin;
    $bins{$bin}->{log_end} = $log_bin+$bstep;
    $bins{$bin}->{start} = $bin;
    $bins{$bin}->{end} = int(10**($log_bin+$bstep))-1;
    $bins{$bin}->{potential} = 0;
    $bins{$bin}->{face} = 0;
    $bins{$bin}->{back} = 0;
    $bins{$bin}->{pos} = 0;
    $bins{$bin}->{neg} = 0;
    $bins{$bin}->{inter} = 0;
    $bins{$bin}->{contigs} = {};
}

my @sorted_bin_keys = sort {$a <=> $b} keys %bins;

###############################################################################################
# contigs
###############################################################################################

print STDERR "reading contig table: $ifn\n";
open(IN, $ifn) || die $ifn;
my $header = <IN>;
my %h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $length = $f[$h{length}];

    # trimmed length
    my $tlength = $length - 2 * $gap;
    next if ($tlength <= 0);

    $contigs{$contig} = {};
    $contigs{$contig}->{length} = $length;

    foreach my $bin (keys %bins) {
	my $start = $bins{$bin}->{start};
	my $end = $bins{$bin}->{end};
	$end  = ($end < $tlength) ? $end : $tlength;
	my $nn = $end - $start;
	if ($start < $end) {
	    $bins{$bin}->{potential} += $nn * ($tlength - ($start+$end)/2) / 1000000;
	    $bins{$bin}->{contigs}->{$contig} = 1;
	}
    }
}

###############################################################################################
# traverse all reads
###############################################################################################

print STDERR "Input dir: $idir\n";
my @ifns = <$idir/P*>;

print STDERR "number of input files: ", scalar(@ifns), "\n";
scalar(@ifns) > 0 or die;

my $count = 0;
my %skipped;
$skipped{short_contig} = 0;
$skipped{gap} = 0;
$skipped{out_of_range} = 0;

foreach my $ifn (@ifns) {
    print STDERR ".";
    open(IN, $ifn) || die $ifn;
    my $header = <IN>;
    my %h = parse_header($header);
    while (my $line = <IN>) {
	chomp $line;
	$count++;
	my @f = split("\t", $line);
	my $contig1 = $f[$h{contig1}];
	my $coord1 = $f[$h{back_coord1}];
	my $strand1 = $f[$h{strand1}];
	my $contig2 = $f[$h{contig2}];
	my $coord2 = $f[$h{back_coord2}];
	my $strand2 = $f[$h{strand2}];

	# skip undefined (short contigs)
	if (!defined($contigs{$contig1}) || !defined($contigs{$contig2})) {
	    $skipped{short_contig}++;
	    next;
	}

	my $length1 = $contigs{$contig1}->{length};
	my $length2 = $contigs{$contig2}->{length};

	# skip if near edges
	if ($coord1 < $gap || $coord2 < $gap || $coord1 > $length1 - $gap || $coord2 > $length2 - $gap) {
	    $skipped{gap}++;
	    next;
	}

	# handle trans
	if ($contig1 ne $contig2) {
	    my $dist1 = $coord1 < $length1/2 ? $coord1 : $length1 - $coord1;
	    my $dist2 = $coord2 < $length2/2 ? $coord2 : $length2 - $coord2;
	    my $tdist = $dist1 + $dist2;
	    my $tindex = binary_search(\@sorted_bin_keys, $tdist, 1);
	    if ($tindex < 0 || $tindex >= scalar(@sorted_bin_keys)) {
		$skipped{out_of_range}++;
		next;
	    }
	    my $tbin = $sorted_bin_keys[$tindex];
	    $bins{$tbin}->{inter}++;
	    next;
	}

	my $contig = $contig1;

	# find bin
	my $dist = abs($coord1-$coord2);
	my $index = binary_search(\@sorted_bin_keys, $dist, 1);
	if ($index < 0 || $index >= scalar(@sorted_bin_keys)) {
	    $skipped{out_of_range}++;
	    next;
	}
	my $bin = $sorted_bin_keys[$index];

	my $lcoord = $coord1 < $coord2 ? $coord1 : $coord2;
	my $rcoord = $coord1 < $coord2 ? $coord2 : $coord1;

	my $lstrand = $coord1 < $coord2 ? $strand1 : $strand2;
	my $rstrand = $coord1 < $coord2 ? $strand2 : $strand1;

	my $type;
	if ($lstrand == 1 && $rstrand == -1) {
	    $type = "face";
	} elsif ($lstrand == -1 && $rstrand == 1) {
	    $type = "back";
	} elsif ($lstrand == 1 && $rstrand == 1) {
	    $type = "pos";
	} else {
	    $type = "neg";
	}
	$bins{$bin}->{$type}++;

	last if ($mreads != -1 && $count >= $mreads);
    }
    close(IN);
    last if ($mreads != -1 && $count >= $mreads);
}
print STDERR "\n";

print STDERR "generating file: $ofn_decay\n";
open(OUT, ">", $ofn_decay) || die $ofn_decay;
print OUT "log_start\tlog_end\tface\tback\tpos\tneg\tinter\tcontig_count\tpotential\n";
foreach my $bin (sort {$a <=> $b} keys %bins) {
    print OUT
	$bins{$bin}->{log_start}, "\t",
	$bins{$bin}->{log_end}, "\t",
	$bins{$bin}->{face}, "\t",
	$bins{$bin}->{back}, "\t",
	$bins{$bin}->{pos}, "\t",
	$bins{$bin}->{neg}, "\t",
	$bins{$bin}->{inter}, "\t",
	scalar(keys %{$bins{$bin}->{contigs}}), "\t",
	$bins{$bin}->{potential}, "\n";
}
close(OUT);


print STDERR "generating file: $ofn_summary\n";
open(OUT, ">", $ofn_summary) || die $ofn_summary;
print OUT "total";
foreach my $key (keys %skipped) {
    print OUT "\t", $key;
}
print OUT "\n";
print OUT "$count";
foreach my $key (keys %skipped) {
    print OUT "\t", $skipped{$key};
}
print OUT "\n";
close(OUT);

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

sub binary_search
{
    my $arr = shift;
    my $value = shift;
    my $above = shift;

    my $left = 0;
    my $right = $#$arr;

    while ($left <= $right) {
        my $mid = ($right + $left) >> 1;
        my $c = $arr->[$mid] <=> $value;
        return $mid if ($c == 0);
        if ($c > 0) {
            $right = $mid - 1;
        } else {
            $left  = $mid + 1;
        }
    }
    $left = -1 if ($left > $#$arr);
    $right = -1 if ($right < 0);
    return $above ? $left : $right;
}

