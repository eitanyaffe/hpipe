#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ofn>\n";
	exit 1;
}

my $ofn = $ARGV[0];

open(OUT, ">", $ofn) or die $ofn;

# skip lines
<STDIN>; <STDIN>; <STDIN>; <STDIN>;

# [S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[FRM]	[TAGS]
# start1  end1	   start2  end2	   len1	   len2	   iden	   strand1  strand2 contig1	   contig2
# 0	  1	   2	   3	   4	   5	   6	   7	    8	    9	           10
# 1	  131576   1	   131576  131576  131576  100.00  1	    1	    contig-100_0   contig-100_0

print OUT "contig1\tstart1\tend1\tstrand1\tcontig2\tstart2\tend2\tstrand2\tidentity\n";
while (my $line = <STDIN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig1 = $f[9];
    my $start1 = ($f[0] < $f[1]) ? $f[0] : $f[1];
    my $end1 = ($f[0] < $f[1]) ? $f[1] : $f[0];
    my $strand1 = $f[7];

    my $contig2 = $f[10];
    my $start2 = ($f[2] < $f[3]) ? $f[2] : $f[3];
    my $end2 = ($f[2] < $f[3]) ? $f[3] : $f[2];
    my $strand2 = $f[8];

    my $identity = $f[6];

    next if ($contig1 eq $contig2);

    print OUT "$contig1\t$start1\t$end1\t$strand1\t$contig2\t$start2\t$end2\t$strand2\t$identity\n";
}
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
