#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);
use POSIX qw( strftime );
use File::stat;
use File::Basename;

if ($#ARGV == -1) {
    print STDERR "usage: $0 <fasta file> <contig anchor table> <only anchor T|F> <odir>\n";
	exit 1;
}

my @pnames = (
    "fasta_ifn",
    "ca_ifn",
    "only_anchor",
    "odir"
    );
my %p;
@p{@pnames} = @ARGV;

print "=============================================\n";
foreach my $key (keys %p) {
    defined($p{$key}) or die "parameter $key not defined (check if all parameters defined";
    print $key, ": ", $p{$key}, "\n";
}
print "=============================================\n";

my $fasta_ifn = $p{fasta_ifn};
my $output_only_anchor = $p{only_anchor} eq "T";
my $ca_ifn = $p{ca_ifn};

#######################################################################################
# read fasta
#######################################################################################

my %contigs;

print "reading fasta into memory: $fasta_ifn\n";
open(IN, $fasta_ifn) || die;
my $contig = "";
my $seq = "";
while (my $line = <IN>) {
    chomp($line);
    if (substr($line,0,1) eq ">") {
	$contigs{$contig} = $seq if ($seq ne "");
	$contig = substr($line, 1);
	$seq = "";
    } else {
	$seq .= $line;
    }
}
$contigs{$contig} = $seq if ($seq ne "");
close(IN);

#######################################################################################
# read ca table
#######################################################################################

my %anchors;

print "reading ca table: $ca_ifn\n";
open(IN, $ca_ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $anchor = $f[$h{anchor}];
    my $contig = $f[$h{contig}];
    my $is_anchor = $f[$h{anchor}] eq $f[$h{contig_anchor}];
    defined($contigs{$contig}) or die "anchor: $anchor, contig: $contig";

    next if ($output_only_anchor && !$is_anchor);
    
    my $seq = $contigs{$contig};

    $anchors{$anchor} = {} if (!defined($anchors{$anchor}));
    $anchors{$anchor}->{$contig} = $seq;
}
close(IN);

#######################################################################################
# generate anchor seq files
#######################################################################################

print "generating anchor files in directory: $p{odir}\n";
foreach my $anchor (keys %anchors) {
    my $ofn = $p{odir}."/".$anchor.".fasta";
    # print "generating anchor file: $ofn\n";
    open(OUT, ">", $ofn) || die;
    foreach my $contig (keys %{$anchors{$anchor}}) {
	my $seq = $anchors{$anchor}->{$contig};
	print OUT ">", $contig, "\n";
	print OUT $seq, "\n";
    }
    close(OUT);
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


