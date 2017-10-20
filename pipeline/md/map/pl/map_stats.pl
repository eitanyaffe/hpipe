#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
    print STDERR "usage: $0 <mapped dir> <parsed dir> <filtered dir> <paired dir> <ofn>\n";
    exit 1;
}

my $mapped_dir = $ARGV[0];
my $parsed_dir = $ARGV[1];
my $filter_dir = $ARGV[2];
my $pair_dir = $ARGV[3];
my $ofn = $ARGV[4];

print "generating file: $ofn\n";
open(OUT, ">", $ofn);
print OUT "field\tcount\n";

line_count($mapped_dir."/R1*", "mapped_r1");
line_count($mapped_dir."/R2*", "mapped_r2");

line_count($parsed_dir."/R1*", "parsed_r1");
line_count($parsed_dir."/R2*", "parsed_r2");


my $command1 = 'cut -f 15 | grep 1';
line_count2($filter_dir."/R1*", "filter_r1", $command1);
line_count2($filter_dir."/R2*", "filter_r2", $command1);

my $command2 = 'cut -f 14,27 | grep -v 0';
line_count2($pair_dir."/P*", "pairs", $command2);

close(OUT);

sub line_count
{
    my ($ifn, $title) = @_;
    print "counting lines for files matching pattern $ifn...\n";
    my @lines = `wc -l $ifn`;
    my $line = $lines[-1];
    my @f = split(/\s+/, $line);
    my $lc = $f[scalar(@lines) == 1 ? 0 : 1];
    print OUT $title, "\t", $lc, "\n";
    print $title, ": ", $lc, "\n";
}

sub line_count2
{
    my ($ifn, $title, $command) = @_;
    my $xcommand = "cat $ifn | $command | wc -l";
    print "counting lines for files matching pattern: $xcommand...\n";
    my @lines = `$xcommand`;
    my $line = $lines[-1];
    my @f = split(/\s+/, $line);
    my $lc = $f[scalar(@lines) == 1 ? 0 : 1];
    print OUT $title, "\t", $lc, "\n";
    print $title, ": ", $lc, "\n";
}

sub get_field_index
{
	my ($ifn, $field) = @_;
	open(IN, $ifn) || die $ifn;
	my $header = <IN>;
	close(IN);
	chomp($header);
	my @f = split("\t", $header);
	my $result = -1;
	for (my $i = 0; $i <= $#f; $i++) {
		$result = $i if $field eq $f[$i];
	}
	return $result;
}
