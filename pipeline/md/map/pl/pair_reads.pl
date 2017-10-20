#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <read1 ifn> <read2 ifn> <ofn> <ofn stats>\n";
	exit 1;
}

my $ifn1 = $ARGV[0];
my $ifn2 = $ARGV[1];
my $ofn = $ARGV[2];
my $ofn_stats = $ARGV[3];

# stats
my $count = 0;
my %reads;

# we discard the sequence field at this time
my @fields = ("contig", "coord", "back_coord", "strand", "edit_dist", "score", "match_length", "cigar", "substitute", "insert", "delete", "clip");

######################################################################################################
# load read1
######################################################################################################

my %stats;
$stats{ok} = 0;
$stats{no_pair} = 0;

foreach my $side ("R1", "R2") {
    my $ifn = $side eq "R1" ? $ifn1 : $ifn2;
    print STDERR "traversing read file: $ifn\n";
    open(IN, $ifn) || die $ifn;
    my $header = <IN>;
    my %h = parse_header($header);

    while (my $line = <IN>) {
	chomp($line);

	$count++;
	print $count, "\n" if ($count % 1000000 == 0);
	my @f = split("\t", $line);

	my $id = $f[$h{id}];
	my $score = $f[$h{score}];
	my $length = $f[$h{match_length}];

	if (!defined($reads{$id})) {
	    $reads{$id} = {};
	    $reads{$id}->{R1} = {};
	    $reads{$id}->{R1}->{max_length} = 0;
	    $reads{$id}->{R2} = {};
	    $reads{$id}->{R2}->{max_length} = 0;
	}

	next if ($length < $reads{$id}->{$side}->{max_length});
	$reads{$id}->{$side}->{max_length} = $length;

	$reads{$id}->{$side}->{fields} = {};
	foreach my $field (@fields) {
	    $reads{$id}->{$side}->{fields}->{$field} = $f[$h{$field}];
	}
    }
    close(IN);
}

######################################################################################################
# output
######################################################################################################

my $paired = 0;

print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "id\t";
foreach my $side (1,2) {
    foreach my $field (@fields) {
	print OUT $field.$side."\t";
    }
}
print OUT "\n";

foreach my $id (keys %reads) {
    if ($reads{$id}->{R1}->{max_length} == 0 || $reads{$id}->{R2}->{max_length} == 0) {
	$stats{no_pair}++;
	next;
    }
    $stats{ok}++;

    print OUT $id;
    foreach my $side (1,2) {
	foreach my $field (@fields) {
	    my $side_str = "R".$side;
	    print OUT "\t".$reads{$id}->{$side_str}->{fields}->{$field};
	}
    }
    print OUT "\n";
    $paired++;
}

print STDERR sprintf("total number of reads: %d\n", $count);
print STDERR sprintf("paired reads: 2x%d=%d (%.1f%%)\n", $paired, $paired*2, 100*$paired/($count/2));

print_hash($ofn_stats, %stats);

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

sub print_hash
{
    my ($ofn, %h) = @_;

    print STDERR "generating file: $ofn\n";
    open (OUT, ">", $ofn) || die $ofn;

    my $first = 1;
    foreach my $key (keys %h) {
	if ($first) {
	    print OUT $key;
	    $first = 0;
	} else {
	    print OUT "\t", $key;
	}
    }
    print OUT "\n";
    $first = 1;
    foreach my $key (keys %h) {
	if ($first) {
	    print OUT $h{$key};
	    $first = 0;
	} else {
	    print OUT "\t", $h{$key};
	}
    }
    print OUT "\n";
    close(OUT);
}

