#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util;

if ($#ARGV == -1) {
    print STDERR "usage: $0 <gene_ifn> <set_field> <map_ifn> <source_field> <target_field> <identity_threshold> <coverage_threshold> <source_ofn> <target_ofn>\n";
	exit 1;
}

my ($gene_ifn, $set_field,
    $map_ifn, $source_field, $target_field,
    $identity_threshold, $coverage_threshold, $source_ofn, $target_ofn) = @ARGV;

print "identity threshold: $identity_threshold\n";
print "coverage threshold: $coverage_threshold\n";

my %genes;
print STDERR "reading gene table: $gene_ifn\n";
open(IN, $gene_ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $set = $f[$h{$set_field}];
    if (!defined($genes{$gene})) {
	$genes{$gene} = {};
	$genes{$gene}->{type} = "missing";
	$genes{$gene}->{found} = 0;
	$genes{$gene}->{sets} = {};
    }
    $genes{$gene}->{sets}->{$set} = 0;
}
close(IN);

print STDERR "reading map: $map_ifn\n";
open(IN, $map_ifn) || die;
$header = <IN>;
%h = parse_header($header);

my %out;

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $source = $f[$h{$source_field}];
    my $target = $f[$h{$target_field}];
    my $identity = $f[$h{identity}];
    my $coverage = $f[$h{coverage}];

    next if (!defined($genes{$source}));

    if ($identity < $identity_threshold || 100*$coverage < $coverage_threshold) {
	$genes{$source}->{type} = "weak" if ($genes{$source}->{type} ne "found");
	next;
    }

    # mark source
    $genes{$source}->{type} = "found";

    # aggregate all targets to report each gene once
    for my $set (keys %{$genes{$source}->{sets}}) {
	$out{$set} = {} if (!defined($out{$set}));
	$out{$set}->{$target} = 0 if (!defined($out{$set}->{$target}));
	$out{$set}->{$target}++;
    }
}
close(IN);

# write source table
print STDERR "writing table: $source_ofn\n";
open(OUT, ">", $source_ofn) || die;
print OUT "set\tgene\ttype\n";
for my $gene (keys %genes) {
    for my $set (keys %{$genes{$gene}->{sets}}) {
	my $type = $genes{$gene}->{type};
	print OUT "$set\t$gene\t$type\n";
    }
}
close(OUT);

# write target table
print STDERR "writing table: $target_ofn\n";
open(OUT, ">", $target_ofn) || die;
print OUT "set\tgene\tcount\n";
for my $set (keys %out) {
for my $gene (keys %{$out{$set}}) {
    my $count = $out{$set}->{$gene};
    print OUT "$set\t$gene\t$count\n";
} }
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

sub read_map
{
    my ($ifn) = @_;
}

