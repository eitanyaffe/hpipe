#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Data::Dumper;
use Switch;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <fends table> <output mat prefix> <input dir> <type>\n";
	exit 1;
}

my $in_fends_fn = $ARGV[0];
my $mat_fn_prefix = $ARGV[1];
my $idir =  $ARGV[2];
my $type = $ARGV[3];
my @ifns = <$idir/*_$type.mat>;

#print STDERR "Input files: ", join(",", @ifns), "\n";
print STDERR "Output prefix: $mat_fn_prefix\n";

##########################################################################################
# read fends file
##########################################################################################

my %fends;

open(IN, $in_fends_fn) || die $in_fends_fn;
print STDERR "Reading input file $in_fends_fn into hash...\n";
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $fend = $f[$h{fend}];
	my $contig = $f[$h{contig}];
	my $coord = $f[$h{coord}];
	my $strand = $f[$h{strand}];
	my $frag = $f[$h{frag}];
	my $frag_len = $f[$h{frag_len}];

	!defined($fends{$fend}) or die "non-unique fend";
	$fends{$fend} = {};
	$fends{$fend}->{contig} = $contig;
	$fends{$fend}->{coord} = $coord;
}
close(IN);

##########################################################################################
# parse pair files
##########################################################################################

my %fend_matrix;
my %fends_covered;

my $line_count = 0;
my $file_count = 0;
my $total_count = @ifns;
foreach my $ifn (@ifns) {
	$file_count++;
	print STDERR "reading file: $ifn ($file_count/$total_count)\n";
	open(IN, $ifn) || die $ifn;

	my $header = <IN>;
	my %h = parse_header($header);
	while (my $line = <IN>) {
		$line_count++;
		# print STDERR "line: $line_count\n" if ($line_count % 100000 == 0);
		chomp $line;
		my @f = split("\t", $line);
		my $fend1 = $f[$h{fend1}];
		my $fend2 = $f[$h{fend2}];

		# track number of covered fends
		$fends_covered{$fend1} = 1;
		$fends_covered{$fend2} = 1;

		my $fend_small = ($fend1 <= $fend2) ? $fend1 : $fend2;
		my $fend_large = ($fend1 <= $fend2) ? $fend2 : $fend1;
		$fend_matrix{$fend_small} = {} if !defined($fend_matrix{$fend_small});
		$fend_matrix{$fend_small}->{$fend_large} = 0 if !defined($fend_matrix{$fend_small}->{$fend_large});
		$fend_matrix{$fend_small}->{$fend_large}++;
	}
	close(IN);
}

######################################################################################################
# write mat file
######################################################################################################

my %fend_stats;
$fend_stats{count} = 0;
$fend_stats{cis_0_1k} = 0;
$fend_stats{cis_1k_10k} = 0;
$fend_stats{cis_10k_100k} = 0;
$fend_stats{cis_100k_1m} = 0;
$fend_stats{cis_1m_max} = 0;
$fend_stats{trans} = 0;

my $mat_fn = $mat_fn_prefix.".mat";
print STDERR "writing file: $mat_fn\n";
open(OUT, ">", $mat_fn) || die;
print OUT "fend1\tfend2\tcount\n";
foreach my $fend1 (sort { $a <=> $b } keys %fend_matrix)
{
	foreach my $fend2 (sort { $a <=> $b } keys %{$fend_matrix{$fend1}})
	{
		$fend_stats{count}++;
		my $count = $fend_matrix{$fend1}->{$fend2};
		print OUT $fend1, "\t" ,$fend2, "\t", $count, "\n";
		
		# trans
		if ($fends{$fend1}->{contig} ne $fends{$fend2}->{contig}) {
			$fend_stats{trans}++;
			next;
		}
		
		# cis
		my $dist = abs($fends{$fend1}->{coord} - $fends{$fend2}->{coord});
		if ($dist < 1000) {
			$fend_stats{cis_0_1k}++;
		} elsif ($dist < 10000) {
			$fend_stats{cis_1k_10k}++;
		} elsif ($dist < 100000) {
			$fend_stats{cis_10k_100k}++;
		} elsif ($dist < 1000000) {
			$fend_stats{cis_100k_1m}++;
		} else {
			$fend_stats{cis_1m_max}++;
		}
	}
}
close(OUT);

######################################################################################################
# Write fend_stats
######################################################################################################

my $fend_stats_fn = $mat_fn_prefix.".fend.stats";
open(OUT, ">", $fend_stats_fn) || die;
print OUT "total_covered\ttotal_pairs\tcis_0_1k\tcis_1k_10k\tcis_10k_100k\tcis_100k_1m\tcis_1m_max\ttrans\n";
print OUT scalar keys %fends_covered, "\t", $fend_stats{count}, "\t"; 
print OUT $fend_stats{cis_0_1k}, "\t", $fend_stats{cis_1k_10k}, "\t", $fend_stats{cis_10k_100k}, "\t", $fend_stats{cis_100k_1m}, "\t", $fend_stats{cis_1m_max}, "\t", $fend_stats{trans}, "\n";
close(OUT);

######################################################################################################
# Subroutines
######################################################################################################

# check if C is between and A,B
sub between
{
	my ($A, $B, $C) = @_;
	$B = $B - $A;
	$C = $C - $A;
	return (($B>$C && $C>0) || ($B<$C && $C<0));
}

sub apprx_lines
{
	my ($fn) = @_;
	my $tmp = "/tmp/".getlogin()."_apprx_lines.tmp";
	system("head -n 100000 $fn > $tmp");
	my $size_head = -s $tmp;
	my $size_all = -s $fn;
	$size_head > 0 or die;
	return (int($size_all/$size_head*100000));
}

sub perc_str
{
	my ($n, $total) = @_;
    return ($n." (".(int(1000 * $n / $total)/10)."%)");
}

sub perc_str2
{
	my ($n, $total) = @_;
    return ((int(1000 * $n / $total)/10)."%");
}

# returns first element above/below value in sorted array
sub binary_search 
{
	my $arr = shift;
	my $value = shift;
	my $above = shift; 

	my $left = 0;
	my $right = $#$arr;
	
	while ($left <= $right) {
		my $mid = ($right + $left) >> 1;
		my $c = $arr->[$mid] <=> $value;
		return $mid if ($c == 0);
		if ($c > 0) {
			$right = $mid - 1;
		} else {
			$left  = $mid + 1;
		}
	}
	$left = -1 if ($left > $#$arr);
	$right = -1 if ($right < 0);
	return (($above eq "+") ? $left : $right);
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
