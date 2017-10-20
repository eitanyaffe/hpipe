#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <input fends> <contig field in fend table> <output fends> <output bins>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $field = $ARGV[1];
my $ofn = $ARGV[2];
my $ofn_bins = $ARGV[3];

#############################################################################################
# fend file
#############################################################################################



print STDERR "Reading file $ifn into hash...\n";
open(IN, $ifn) || die $ifn;
my $header = <IN>;
my %h = parse_header($header);

print STDERR "Writing output file: $ofn\n";
open(OUT, ">", $ofn) || die;
chomp($header);
print OUT $header, "\tcontig_bin\n";

my $index = 1;
my %contigs;
my %bins;
while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $contig = $f[$h{$field}];

	if (!defined($contigs{$contig})) {
	    $contigs{$contig} = $index;
	    $bins{$index} = $contig;
	    $index++;
	}
	
	print OUT $line, "\t", $contigs{$contig}, "\n";
}
close(IN);
close(OUT);

print STDERR "Writing output bin file: $ofn_bins\n";
open(OUT, ">", $ofn_bins) || die;
print OUT "index\tcontig\n";
foreach my $index (sort {$a <=> $b} keys %bins) {
    my $contig = $bins{$index};
    print OUT "$index\t$contig\n";
}
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
