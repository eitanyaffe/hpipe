#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <map_ifn> <min identity> <max identity> <min coverage> <ofn>\n";
	exit 1;
}

my ($ifn, $min_identity, $max_identity, $min_coverage, $ofn) = @ARGV;

###################################################################################################################
# output graph
###################################################################################################################

my %genes;
print STDERR "reading map: $ifn\n";
open(IN, $ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $identity = $f[$h{identity}];
    my $coverage = $f[$h{coverage}];
    my $gene = $f[$h{gene1}];
    my $ogene = $f[$h{gene2}];
    next if (100*$coverage < $min_coverage || $gene eq $ogene || $identity > $max_identity || $identity < $min_identity);
    if (!defined($genes{$gene})) {
	$genes{$gene} = {};
	$genes{$gene}->{count} = 0;
	$genes{$gene}->{min_identity} = 100;
	$genes{$gene}->{max_identity} = 0;
    }
    $genes{$gene}->{count}++;
    $genes{$gene}->{min_identity} = $identity if ($identity < $genes{$gene}->{min_identity});
    $genes{$gene}->{max_identity} = $identity if ($identity > $genes{$gene}->{max_identity});
}
close(IN);

###################################################################################################################
# output graph
###################################################################################################################

print STDERR "generating output: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT "gene\tcount\tmin_identity\tmax_identity\n";
foreach my $gene (keys %genes) {
    print OUT $gene, "\t", $genes{$gene}->{count}, "\t", $genes{$gene}->{min_identity}, "\t", $genes{$gene}->{max_identity}, "\n";
}
close(OUT);

#######################################################################################
# utils
#######################################################################################

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
