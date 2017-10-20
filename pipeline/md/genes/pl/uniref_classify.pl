#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn genes> <class table> <ofn>\n";
	exit 1;
}

my $ifn_genes = $ARGV[0];
my $ifn_class = $ARGV[1];
my $ofn = $ARGV[2];

#######################################################################################
# read class table
#######################################################################################

# contig length table
my %class_ht;

print STDERR "reading class table: $ifn_class\n";
open(IN, $ifn_class) || die $ifn_class;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    $class_ht{$f[$h{regex}]} = $f[$h{class}];
}
close(IN);

#######################################################################################
# traverse table
#######################################################################################

print STDERR "going over table: $ifn_genes\n";
open(IN, $ifn_genes) || die $ifn_genes;

print STDERR "creating output: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
$header = <IN>;
chomp($header);
%h = parse_header($header);

print OUT $header, "\tclass\n";
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $class = "none";

    # search by description
    my $desc = $f[$h{prot_desc}];
    foreach my $regex (keys %class_ht) {
	if ($desc =~ /$regex/i) {
	    $class = $class_ht{$regex};
	    last;
	}
    }

    # if no hit use count as well
    $class = "conserved" if ($class eq "none" && $f[$h{uniref_count}] > 1);
    
    print OUT $line, "\t", $class, "\n";
}
close(IN);
close(OUT);

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
