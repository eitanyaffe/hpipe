#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);

if ($#ARGV == -1) {
    print STDERR "usage: $0 <ifn table> <ref dir> <anchor dir> <summary script>\n";
	exit 1;
}

my @pnames = ("ifn", "ref_dir", "anchor_dir", "script");
my %p = parse_parameters(\@pnames, \@ARGV);

#######################################################################################
# read genome table
#######################################################################################

my %refs;
my %anchors;

print "reading table: $p{ifn}\n";
open(IN, $p{ifn}) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line, -1);
    my $id = $f[$h{accession}];
    my $anchor = $f[$h{anchor}];
    $refs{$id} = 1;
    $anchors{$anchor} = 1;
}
close(IN);

#######################################################################################
# generate fasta summary files
#######################################################################################

foreach my $anchor (keys %anchors) {
    my $ifn = $p{anchor_dir}."/".$anchor.".fasta";
    my $ofn = $p{anchor_dir}."/".$anchor.".contig_table";
    my $command ="cat ".$ifn." | ".$p{script}." > ".$ofn;
    # print "command: $command\n";
    system($command) == 0 or die;
}

foreach my $ref (keys %refs) {
    my $ifn = $p{ref_dir}."/".$ref."/genome.fasta";
    my $ofn = $p{ref_dir}."/".$ref."/contig_table";

    # skip if file is missing
    (-e $ifn) or next;

    my $command ="cat ".$ifn." | ".$p{script}." > ".$ofn;
    # print "command: $command\n";
    system($command) == 0 or die;
}

#######################################################################################
# utils
#######################################################################################

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

sub parse_parameters
{
    my ($p_pnames, $p_argv) = @_;
    my @pnames = @{$p_pnames};
    my @argv = @{$p_argv};

    my %p;
    @p{@pnames} = @argv;
    my %po;
    my $N = scalar(@argv);
    my @pi = (1..$N);
    @po{@pi} = @pnames;

    print "=============================================\n";
    foreach my $index (sort {$a <=> $b} keys %po) {
	my $key = $po{$index};
	defined($p{$key}) or die "parameter $key not defined (check if all parameters defined)";
	print $key, ": ", $p{$key}, "\n";
    }
    print "=============================================\n";

    return (%p);
}
