#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Data::Dumper;
use List::Util qw(first);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <base input file> <appended input file> <output file> <value field> <default value> <symmetric keys T|F> <mfn>";
	print STDERR "      appends a single column of a table to another table according to a key defined by multiple columns\n";
	exit 1;
}

my $mfn = $ARGV[6];
my $wd = $ARGV[7];

print STDERR "Reading model table file: $mfn\n";
open(IN, $mfn) || die;
my $header = <IN>;
my %h = parse_header($header);

my %model;
my @fields;
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $rfield = $f[$h{raw_field}];

    !defined($model{$rfield}) or die;
    $model{$rfield} = {};
    $model{$rfield}->{field} = $f[$h{field}];
    $model{$rfield}->{size} = $f[$h{size}];
    $model{$rfield}->{type} = $f[$h{type}];
    push(@fields, $f[$h{field}]."1");
    push(@fields, $f[$h{field}]."2");
}

my $command = $wd."/pl/append_table.pl ".join(" ", @ARGV[0..$#ARGV-2])." ".join(" ", @fields);
print "running command: ", $command, "\n";
system($command) == 0 or die;

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

sub apprx_lines
{
	my ($fn) = @_;
	my $tmp = "/tmp/_apprx_lines.".int(rand(1000000));
	system("head -n 100000 $fn > $tmp");
	my $size_head = -s $tmp;
	my $size_all = -s $fn;
	return (int($size_all/$size_head*100000));
}

