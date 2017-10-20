#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
    print STDERR "usage: $0 <gene anchor table> <gene matrix> <odir>\n";
    exit 1;
}

my $ga_ifn = $ARGV[0];
my $gene_mat_ifn = $ARGV[1];
my $odir = $ARGV[2];

############################################################################################################
# gene table
############################################################################################################

my %genes;

print STDERR "reading gene file: $ga_ifn\n";
open(IN, $ga_ifn) || die $ga_ifn;
my $header = <IN>;
my %h = parse_header($header);

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{cgene}];
    my $anchor = $f[$h{anchor}];
    $genes{$gene} = {} if (!defined($genes{$gene}));
    $genes{$gene}->{$anchor} = 1;
}
close(IN);

############################################################################################################
# gene mat table
############################################################################################################

print STDERR "reading file: $gene_mat_ifn\n";
open(IN, $gene_mat_ifn) || die $gene_mat_ifn;
$header = <IN>;
%h = parse_header($header);

my %anchors;

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $cgene1 = $f[$h{cgene1}];
    my $cgene2 = $f[$h{cgene2}];
    my $iden = $f[$h{identity}];
    my $coverage = $f[$h{coverage}];

    next if (!defined($genes{$cgene1}));
    next if (!defined($genes{$cgene2}));

    my @qs = keys %{$genes{$cgene1}};
    my @ts = keys %{$genes{$cgene2}};
    foreach my $anchor_q (@qs) {
    foreach my $anchor_t (@ts) {
	$anchors{$anchor_q} = {} if (!defined($anchors{$anchor_q}));
	$anchors{$anchor_q}->{$anchor_t} = {} if (!defined($anchors{$anchor_q}->{$anchor_t}));

	my $rline = $cgene1."\t".$cgene2."\t".$iden."\t".$coverage;

	my $key = $cgene1."_".$cgene2."_".scalar (keys %{$anchors{$anchor_q}->{$anchor_t}});
	!defined($anchors{$anchor_q}->{$anchor_t}->{$key}) or die $key;
	$anchors{$anchor_q}->{$anchor_t}->{$key} = $rline;
    } }
}
close(IN);

print STDERR "generating files in: $odir\n";
foreach my $anchor_q (keys %anchors) {
foreach my $anchor_t (keys %{$anchors{$anchor_q} }) {
    my $ofn = $odir."/".$anchor_q."_".$anchor_t;
    # print STDERR "generating file: $ofn\n";
    open(OUT, ">", $ofn) || die $ofn;
    print OUT "cgene1\tcgene2\tidentity\tcoverage\n";
    foreach my $key (keys %{$anchors{$anchor_q}->{$anchor_t} } ) {
	print OUT $anchors{$anchor_q}->{$anchor_t}->{$key}, "\n";
    }
    close(OUT);
} }

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


