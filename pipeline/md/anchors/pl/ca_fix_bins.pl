#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ofn = $ARGV[1];

#############################################################################################
# fend file
#############################################################################################


open(IN, $ifn) || die;
print STDERR "Reading file $ifn\n";
my $header = <IN>;
my %h = parse_header($header);

open(OUT, ">", $ofn) || die;
print STDERR "writing $ofn\n";
print OUT $header;

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $anchor = $f[$h{anchor}];
    $f[$h{anchor_bin}] = 0 if ($anchor == 0);
    print OUT join("\t", @f), "\n";
}
close(IN);
close(OUT);

######################################################################################################
# Subroutines
######################################################################################################

sub apprx_lines
{
	my ($fn) = @_;
	my $tmp = "/tmp/".$$."_apprx_lines.tmp";
	system("head -n 100000 $fn > $tmp");
	my $size_head = -s $tmp;
	my $size_all = -s $fn;
	return (int($size_all/$size_head*100000));
}

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
