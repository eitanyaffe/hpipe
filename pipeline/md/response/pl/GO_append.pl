#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <gene/element table> <Gene2GO table> <GO tree> <ofn>\n";
	exit 1;
}

my $ifn_genes = $ARGV[0];
my $ifn_gene2go_table = $ARGV[1];
my $ifn_go_tree = $ARGV[2];
my $ofn = $ARGV[3];

###############################################################################################
# go tree
###############################################################################################

my %tree;

print STDERR "reading table: $ifn_go_tree\n";
open(IN, $ifn_go_tree) || die $ifn_go_tree;
my $header = <IN>;
my %h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $id = $f[$h{id}];
    $tree{$id} = {};
    $tree{$id}->{root} = $f[$h{root}] eq "T";
    $tree{$id}->{desc} = $f[$h{desc}];
    $tree{$id}->{type} = $f[$h{type}];
}
close(IN);
print "GO tree size: ", scalar(keys %tree), "\n";

###############################################################################################
# Gene2GO table
###############################################################################################

my %genes;

print STDERR "reading table: $ifn_gene2go_table\n";
open(IN, $ifn_gene2go_table) || die $ifn_gene2go_table;
$header = <IN>;
%h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $identity = $f[$h{identity}];
    $genes{$gene} = $f[$h{GO}];
}
close(IN);

###############################################################################################
# genes
###############################################################################################

print STDERR "reading table: $ifn_genes\n";
open(IN, $ifn_genes) || die $ifn_genes;
$header = <IN>;
%h = parse_header($header);

chomp($header);
print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "gene\tuniref\tidentity\tprot_desc\tGO\tGO_type\tGO_desc\n";

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $uniref = $f[$h{uniref}];
    next if ($uniref eq "NA");
    defined($genes{$gene}) or die $gene;
    my @GOs = split(";", $genes{$gene});
    for my $GO (@GOs) {
	next if ($GO eq "NA");
	defined($tree{$GO}) or die;
	print OUT
	    $f[$h{gene}], "\t",
	    $f[$h{uniref}], "\t",
	    $f[$h{identity}], "\t",
	    $f[$h{prot_desc}], "\t",
	    $GO, "\t",
	    $tree{$GO}->{type}, "\t",
	    $tree{$GO}->{desc}, "\n";
    }
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
