#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <input pattern> <fold> <ofn>\n";
	exit 1;
}

my $fold = $ARGV[0];
my $ofn = $ARGV[1];
my $pattern = $ARGV[2];
shift; shift; shift;
my @idirs = @ARGV;

# from log to real factor
my $factor = int(2 ** $fold);
print "subsample x-fold: $factor\n";
my @ifns;
foreach my $dir (@idirs) {
    my @ix = <$dir/$pattern>;
    push(@ifns, @ix);
}
print "files found: ".scalar(@ifns)."\n";
print "generating file: $ofn\n";
open(OUT, ">", $ofn) || die "$ofn";

my $count = 0;
foreach my $ifn (@ifns) {
    open(IN, $ifn) || die $ifn;

    my $selected = 1;
    while (my $line = <IN>) {
	if (substr($line,0,1) eq "@") {
	    $count++;
	    $selected = ($count % $factor == 0);
	}
	print OUT $line if ($selected);
    }
    close(IN);
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
