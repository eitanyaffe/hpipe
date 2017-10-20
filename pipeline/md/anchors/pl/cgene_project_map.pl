#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util;

if ($#ARGV == -1) {
    print STDERR "usage: $0 <map_ifn> <identity_threshold> <coverage_threshold> <ofn>\n";
	exit 1;
}

my ($ifn, $identity_threshold, $coverage_threshold, $ofn) = @ARGV;

print "identity threshold: $identity_threshold\n";
print "coverage threshold: $coverage_threshold\n";

print STDERR "reading map: $ifn\n";
open(IN, $ifn) || die;
my $header = <IN>;
my %h = parse_header($header);

print STDERR "writing table: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT $header;

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $identity = $f[$h{identity}];
    my $coverage = $f[$h{coverage}];

    next if ($identity < $identity_threshold || 100*$coverage < $coverage_threshold);
    print OUT $line, "\n";
}
close(IN);
close(OUT);

#######################################################################################
# utils
#######################################################################################

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

sub read_map
{
    my ($ifn) = @_;
}

