#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Data::Dumper;
use List::Util qw(first);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <base input file> <appended input file> <output file> <value field> <default value> <symmetric keys T|F> <key1> [key2 ...] \n";
	print STDERR "      appends a single column of a table to another table according to a key defined by multiple columns\n";
	exit 1;
}

my $bfn = $ARGV[0];
my $afn = $ARGV[1];
my $ofn = $ARGV[2];
my $value_field = $ARGV[3];
my $default_value = $ARGV[4];
my $sym_key = $ARGV[5];
shift; shift; shift; shift; shift; shift;

($sym_key eq "T") or ($sym_key eq "F") or die "'symmetric key' must be T or F: $sym_key";

my @keys;
while ($#ARGV >= 0)
{
	push(@keys, $ARGV[0]); shift;
}
print STDERR "keys: ", join(",", @keys), "\n";
# all fields
my %fields;

# appended values
my %table;


# traverse appended file
open(IN, $afn) || die $afn;

my $appr_lines = apprx_lines($afn);
print STDERR "traversing file $afn, with about ".int($appr_lines/1000000)."M lines\n";

my $header = <IN>;
my %h = parse_header($header);

my $count = 1;
while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $key = "";
	$key .= ":".$f[$h{$_}] foreach (@keys);
	$table{$key} = $f[$h{$value_field}];

	if ($sym_key eq "T")
	{
		$key = "";
		$key = ":".$f[$h{$_}].$key foreach (@keys);
		$table{$key} = $f[$h{$value_field}];
	}

	print STDERR "line: $count\n" if ($count % 1000000 == 0);
	$count++;
}
close(IN);

open(IN, $bfn) || die $bfn;
$appr_lines = apprx_lines($bfn);
print STDERR "traversing file $bfn, with about ".int($appr_lines/1000000)."M lines\n";

# write output
print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die;

$header = <IN>;
%h = parse_header($header);
chomp($header);
print OUT $header, "\t", $value_field, "\n";
$count = 1;
while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $key = "";
	$key .= ":".$f[$h{$_}] foreach (@keys);
	my $value = defined($table{$key}) ? $table{$key} : $default_value;
	print OUT $line, "\t", $value, "\n";

	print STDERR "line: $count\n" if ($count % 1000000 == 0);
	$count++;
}
close(IN);
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

sub apprx_lines
{
	my ($fn) = @_;
	my $tmp = "/tmp/_apprx_lines.".int(rand(1000000));
	system("head -n 100000 $fn > $tmp");
	my $size_head = -s $tmp;
	my $size_all = -s $fn;
	return (int($size_all/$size_head*100000));
}

