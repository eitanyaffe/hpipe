#!/usr/bin/perl

use strict;
use Data::Dumper;
use warnings FATAL => qw( all );

if ($#ARGV == -1) {
	print STDERR "usage: $0 <file prefix> <scope> <symmetric T|F> <fend suffix> <mfn> <work dir>\n";
	print STDERR "  otag can be set to 'none'\n";
	exit 1;
}

my $mfn = $ARGV[4];
my $wd = $ARGV[5];


open(IN, $mfn) || die;
print STDERR "Reading model table file: $mfn\n";
my $header = <IN>;
my %h = parse_header($header);

my %fields;
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $field = $f[$h{field}];
    !defined($fields{$field}) or die;
    $fields{$field} = 1;
}

my $command = $wd."/pl/compute_n.pl ".join(" ", @ARGV[0..$#ARGV-2])." ".join(" ", sort keys %fields);
print "running command: ", $command, "\n";
system($command) == 0 or die;

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
