#!/usr/bin/perl

use strict;
use Data::Dumper;
use warnings FATAL => qw( all );

if ($#ARGV == -1) {
	print STDERR "usage: $0 <file prefix> <scope> <symmetric T|F> <fend suffix> <field1> [field2 ...]\n";
	exit 1;
}

my $iprefix = $ARGV[0];
my $scope = $ARGV[1];
my $sym_key = $ARGV[2];
my $fend_suffix = $ARGV[3];

($sym_key eq "T") or ($sym_key eq "F") or die "'symmetric key' must be T or F, found $sym_key";

($scope eq "inter_anchor") or ($scope eq "intra_anchor") or die "unknown scope: $scope";

print STDERR "scope: $scope\n";

my @fields;
shift;shift;shift;shift;
while ($#ARGV >= 0) {
	push(@fields, $ARGV[0]); shift;
}
print STDERR "Computing bin counts for fields: ", join(",", @fields), "\n";

my $site_fn = $iprefix.".".$fend_suffix;
my $mat_fn = $iprefix.".mat";
my $out_fn = $iprefix.".counts";
my $blacklist_fn = $iprefix.".blacklist";

# init fend blacklist if one exists
our %blacklist;
if (-e $blacklist_fn)
{
	open(IN, $blacklist_fn) || die;
	<IN>;
	while (my $line = <IN>) {
		chomp $line;
		$blacklist{$line} = 1;
	}
	close(IN);
	print STDERR "Blacklisted fends: ", scalar keys %blacklist, "\n";
}

our %fends;

open(IN, $site_fn) || die;
print STDERR "Reading file $site_fn into hash...\n";
my $header = <IN>;
my %hhash = parse_header($header);

while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $fend = $f[$hhash{fend}];
	!defined($fends{$fend}) or die "non-unique fend";
	next if defined($blacklist{$fend});
	$fends{$fend} = {};
	$fends{$fend}->{fend} = $f[$hhash{fend}];
	$fends{$fend}->{contig} = $f[$hhash{contig}];
	$fends{$fend}->{anchor} = $f[$hhash{anchor}];
	$fends{$fend}->{coord} = $f[$hhash{coord}];
	foreach my $field (@fields) {
		$fends{$fend}->{$field} = $f[$hhash{$field}];
	}
}
close(IN);

my $appr_lines = apprx_lines($mat_fn);
print STDERR "traversing file $mat_fn, with about ".int($appr_lines/1000000)."M lines\n";

our %bin_counters;

my $missing_count=0;
my $used_count=0;

my $count = 0;
open(IN, $mat_fn) || die;
$header = <IN>;
%hhash = parse_header($header);
while (my $line = <IN>)
{
	$count++;
	print STDERR "line: $count\n" if ($count % 1000000 == 0);

	chomp $line;
	my @f = split("\t", $line);

	my $fend1 = $f[$hhash{fend1}];
	my $fend2 = $f[$hhash{fend2}];

	if (!defined($fends{$fend1}) or !defined($fends{$fend2}))
	{
		$missing_count++;
		next;
	}

	my $coord1 = $fends{$fend1}->{coord};
	my $coord2 = $fends{$fend2}->{coord};

	my $contig1 = $fends{$fend1}->{contig};
	my $contig2 = $fends{$fend2}->{contig};

	my $anchor1 = $fends{$fend1}->{anchor};
	my $anchor2 = $fends{$fend2}->{anchor};

	next if ($anchor1 == 0 || $anchor2 == 0 || ($contig1 eq $contig2));

	if ($scope eq "inter_anchor") {
	    next if ($anchor1 eq $anchor2);
	} elsif ($scope eq "intra_anchor") {
	    next if ($anchor1 ne $anchor2);
	} else {
	    die "unknown scope: $scope";
	}

	$used_count++;
	my @bins1;
	my @bins2;
	foreach my $field (@fields) {
		defined($fends{$fend1}->{$field}) and defined($fends{$fend2}->{$field}) or die;
		push(@bins1, $fends{$fend1}->{$field});
		push(@bins2, $fends{$fend2}->{$field});
	}
	my $abin = make_bin(join("_", @bins1), join("_", @bins2), 0);
	$bin_counters{$abin} = 0 if !defined($bin_counters{$abin});
	$bin_counters{$abin}++;

	# handle symmetric
	if ($sym_key eq "T") {
		my $sbin = make_bin(join("_", @bins1), join("_", @bins2), 1);
		$bin_counters{$sbin} = 0 if !defined($bin_counters{$sbin});
		$bin_counters{$sbin}++ if ($sbin ne $abin);
	}
}
close(IN);
print STDERR "number of used pairs: $used_count\n";
print STDERR "number of pairs discarded some where along the pipeline: $missing_count\n";
$used_count > 0 or die "no pairs found";

open(OUT, ">", $out_fn) || die;
print STDERR "Creating file $out_fn\n";

print OUT $_."1\t" foreach (@fields);
print OUT $_."2\t" foreach (@fields);
print OUT "count\n";

foreach my $bin (sort sort_bins keys %bin_counters)
{
	my @vs = split("_", $bin);
	my $counts = $bin_counters{$bin};

	print OUT join("\t", @vs);
	print OUT "\t$counts\n";
}
close(OUT);


######################################################################################################
# Subroutines
######################################################################################################

sub sort_bins
{
	my @as = split("_", $a);
	my @bs = split("_", $b);
	my $result = $as[0] <=> $bs[0];
	for (my $i = 1; $i < @as; $i++) {
		$result = $result || ($as[$i] <=> $bs[$i]);
	}
	return $result;
}

sub rsort_bins
{
	my @as = split("_", $a);
	my @bs = split("_", $b);
	my $result = $bs[0] <=> $as[0];
	for (my $i = 1; $i < @as; $i++) {
		$result = $result || ($bs[$i] <=> $as[$i]);
	}
	return $result;
}

sub make_bin
{
	my ($bin1, $bin2, $reverse) = @_;


	return ($reverse) ? join("_", sort (rsort_bins ($bin1, $bin2))) :
                        join("_", sort (sort_bins ($bin1, $bin2)));
}

sub apprx_lines
{
	my ($fn) = @_;
	my $tmp = "/tmp/_apprx_lines.".int(rand(1000000));
	system("head -n 100000 $fn > $tmp") == 0 or die;
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
