#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util;

if ($#ARGV == -1) {
    print STDERR "usage: $0 <ifn> <ofn>\n";
	exit 1;
}

my ($ifn, $ofn) = @ARGV;

print STDERR "reading gene table: $ifn\n";
open(IN, $ifn) || die;
my $header = <IN>;
my %h = parse_header($header);

print STDERR "writing table: $ofn\n";
open(OUT, ">", $ofn) || die;
chomp($header);
print OUT $header, "\tgenome\n";

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $contig = $f[$h{contig}];
    my $genome = $contig;
    print OUT $line, "\t", $genome, "\n";
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

