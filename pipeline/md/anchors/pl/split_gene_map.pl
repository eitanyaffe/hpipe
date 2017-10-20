#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
    print STDERR "usage: $0 <gene table1> <set field1> <set null value1> <gene table2> <set field1> <set null value2> <map> <odir>\n";
    exit 1;
}

my ($gene_ifn1, $field1, $null_value1, $gene_ifn2, $field2, $null_value2, $map_ifn, $base_odir) = @ARGV;

my $id_field = "gene";

############################################################################################################
# gene table
############################################################################################################

my %genes1;
print STDERR "reading gene file: $gene_ifn1\n";
open(IN, $gene_ifn1) or die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{$id_field}];
    my $set = $f[$h{$field1}];
    next if ($set eq $null_value1);
    $genes1{$gene} = {} if (!defined($genes1{$gene}));
    $genes1{$gene}->{$set} = 0;
}
close(IN);

my %genes2;
print STDERR "reading gene file: $gene_ifn2\n";
open(IN, $gene_ifn2) or die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{$id_field}];
    my $set = $f[$h{$field2}];
    next if ($set eq $null_value2);
    $genes2{$gene} = {} if (!defined($genes2{$gene}));
    $genes2{$gene}->{$set} = 0;
}
close(IN);

############################################################################################################
# gene mat table
############################################################################################################

print STDERR "reading file: $map_ifn\n";
open(IN, $map_ifn) or die;
$header = <IN>;
%h = parse_header($header);

my %map;
my %map_orphan_q;
my %map_orphan_t;

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $gene1 = $f[$h{gene1}];
    my $gene2 = $f[$h{gene2}];

    my $iden = $f[$h{identity}];
    my $coverage = $f[$h{coverage}];

    next if (!defined($genes1{$gene1}) && !defined($genes2{$gene2}));

    my @qs = defined($genes1{$gene1}) ? keys %{$genes1{$gene1}} : ();
    my @ts = defined($genes2{$gene2}) ? keys %{$genes2{$gene2}} : ();

    if (!defined($genes2{$gene2})) {
	foreach my $set_q (@qs) {
	    $map_orphan_q{$set_q} = {} if (!defined($map_orphan_q{$set_q}));
	    my $rline = $gene1."\t".$gene2."\t".$iden."\t".$coverage;
	    my $key = $gene1."_".$gene2."_".scalar (keys %{$map_orphan_q{$set_q}});
	    !defined($map_orphan_q{$set_q}->{$key}) or die $key;
	    $map_orphan_q{$set_q}->{$key} = $rline;
	}
	next;
    }

    if (!defined($genes1{$gene1})) {
	foreach my $set_t (@ts) {
	    $map_orphan_t{$set_t} = {} if (!defined($map_orphan_t{$set_t}));
	    my $rline = $gene1."\t".$gene2."\t".$iden."\t".$coverage;
	    my $key = $gene1."_".$gene2."_".scalar (keys %{$map_orphan_t{$set_t}});
	    !defined($map_orphan_t{$set_t}->{$key}) or die $key;
	    $map_orphan_t{$set_t}->{$key} = $rline;
	}
	next;
    }

    foreach my $set_q (@qs) {
    foreach my $set_t (@ts) {
	$map{$set_q} = {} if (!defined($map{$set_q}));
	$map{$set_q}->{$set_t} = {} if (!defined($map{$set_q}->{$set_t}));

	my $rline = $gene1."\t".$gene2."\t".$iden."\t".$coverage;

	my $key = $gene1."_".$gene2."_".scalar (keys %{$map{$set_q}->{$set_t}});
	!defined($map{$set_q}->{$set_t}->{$key}) or die $key;
	$map{$set_q}->{$set_t}->{$key} = $rline;
    } }
}
close(IN);

print STDERR "generating files in: $base_odir\n";
foreach my $set_q (keys %map) {
foreach my $set_t (keys %{$map{$set_q} }) {
    my $odir = $base_odir."/".$set_q;
    system(sprintf("mkdir -p %s", $odir)) == 0 or die;
    my $ofn = $odir."/".$set_t;
    open(OUT, ">", $ofn) || die $ofn;
    print OUT $id_field."1\t".$id_field."2\tidentity\tcoverage\n";
    foreach my $key (keys %{$map{$set_q}->{$set_t} } ) {
	print OUT $map{$set_q}->{$set_t}->{$key}, "\n";
    }
    close(OUT);
} }

foreach my $set_q (keys %map_orphan_q) {
    my $odir = $base_odir."/".$set_q;
    system(sprintf("mkdir -p %s", $odir)) == 0 or die;
    my $ofn = $odir."/NONE";
    open(OUT, ">", $ofn) || die $ofn;
    print OUT $id_field."1\t".$id_field."2\tidentity\tcoverage\n";
    foreach my $key (keys %{$map_orphan_q{$set_q}} ) {
	print OUT $map_orphan_q{$set_q}->{$key}, "\n";
    }
    close(OUT);
}

foreach my $set_t (keys %map_orphan_t) {
    my $odir = $base_odir."/NONE";
    system(sprintf("mkdir -p %s", $odir)) == 0 or die;
    my $ofn = $odir."/".$set_t;
    open(OUT, ">", $ofn) || die $ofn;
    print OUT $id_field."1\t".$id_field."2\tidentity\tcoverage\n";
    foreach my $key (keys %{$map_orphan_t{$set_t}} ) {
	print OUT $map_orphan_t{$set_t}->{$key}, "\n";
    }
    close(OUT);
}

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


