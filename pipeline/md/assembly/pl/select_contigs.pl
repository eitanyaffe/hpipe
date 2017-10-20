#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <table> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $table = $ARGV[1];
my $ofn = $ARGV[2];

my %contigs;

# read table
print STDERR "reading table: $table\n";
open(IN, $table) || die $table;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    $contigs{$f[$h{contig}]} = 1;
}
close(IN);

print STDERR "traversing file: $ifn\n";
open(IN, $ifn) || die $ifn;
print STDERR "generating output file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;

my $selected = 0;
while (my $line = <IN>) {
    chomp($line);
    if (substr($line, 0, 1) eq ">") {
	my @f = split(" ", substr($line,1));
	my $contig = $f[0];
	$selected = defined($contigs{$contig});
    }
    print OUT $line, "\n" if ($selected);
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
