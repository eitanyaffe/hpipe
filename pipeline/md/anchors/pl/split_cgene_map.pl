#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
    print STDERR "usage: $0 <cgene table1> <set field1> <cgene table2> <set field1> <map> <odir>\n";
    exit 1;
}

my ($cgene_ifn1, $field1, $cgene_ifn2, $field2, $map_ifn, $odir) = @ARGV;

my $id_field = "gene";

############################################################################################################
# gene table
############################################################################################################

my %genes1;
print STDERR "reading gene file: $cgene_ifn1\n";
open(IN, $cgene_ifn1) or die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{$id_field}];
    my $set = $f[$h{$field1}];
    $genes1{$gene} = {} if (!defined($genes1{$gene}));
    $genes1{$gene}->{$set} = 0;
}
close(IN);

my %genes2;
print STDERR "reading gene file: $cgene_ifn2\n";
open(IN, $cgene_ifn2) or die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{$id_field}];
    my $set = $f[$h{$field2}];
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

    my $cgene1 = $f[$h{cgene1}];
    my $cgene2 = $f[$h{cgene2}];

    my $iden = $f[$h{identity}];
    my $coverage = $f[$h{coverage}];

    next if (!defined($genes1{$cgene1}) && !defined($genes2{$cgene2}));

    my @qs = defined($genes1{$cgene1}) ? keys %{$genes1{$cgene1}} : ();
    my @ts = defined($genes2{$cgene2}) ? keys %{$genes2{$cgene2}} : ();

    if (!defined($genes2{$cgene2})) {
	foreach my $set_q (@qs) {
	    $map_orphan_q{$set_q} = {} if (!defined($map_orphan_q{$set_q}));
	    my $rline = $cgene1."\t".$cgene2."\t".$iden."\t".$coverage;
	    my $key = $cgene1."_".$cgene2."_".scalar (keys %{$map_orphan_q{$set_q}});
	    !defined($map_orphan_q{$set_q}->{$key}) or die $key;
	    $map_orphan_q{$set_q}->{$key} = $rline;
	}
	next;
    }

    if (!defined($genes1{$cgene1})) {
	foreach my $set_t (@ts) {
	    $map_orphan_t{$set_t} = {} if (!defined($map_orphan_t{$set_t}));
	    my $rline = $cgene1."\t".$cgene2."\t".$iden."\t".$coverage;
	    my $key = $cgene1."_".$cgene2."_".scalar (keys %{$map_orphan_t{$set_t}});
	    !defined($map_orphan_t{$set_t}->{$key}) or die $key;
	    $map_orphan_t{$set_t}->{$key} = $rline;
	}
	next;
    }

    foreach my $set_q (@qs) {
    foreach my $set_t (@ts) {
	$map{$set_q} = {} if (!defined($map{$set_q}));
	$map{$set_q}->{$set_t} = {} if (!defined($map{$set_q}->{$set_t}));

	my $rline = $cgene1."\t".$cgene2."\t".$iden."\t".$coverage;

	my $key = $cgene1."_".$cgene2."_".scalar (keys %{$map{$set_q}->{$set_t}});
	!defined($map{$set_q}->{$set_t}->{$key}) or die $key;
	$map{$set_q}->{$set_t}->{$key} = $rline;
    } }
}
close(IN);

print STDERR "generating files in: $odir\n";

foreach my $set_q (keys %map) {
foreach my $set_t (keys %{$map{$set_q} }) {
    my $ofn = $odir."/".$set_q."_".$set_t;
    open(OUT, ">", $ofn) || die $ofn;
    print OUT $id_field."1\t".$id_field."2\tidentity\tcoverage\n";
    foreach my $key (keys %{$map{$set_q}->{$set_t} } ) {
	print OUT $map{$set_q}->{$set_t}->{$key}, "\n";
    }
    close(OUT);
} }

foreach my $set_q (keys %map_orphan_q) {
    my $ofn = $odir."/".$set_q."_NONE";
    open(OUT, ">", $ofn) || die $ofn;
    print OUT $id_field."1\t".$id_field."2\tidentity\tcoverage\n";
    foreach my $key (keys %{$map_orphan_q{$set_q}} ) {
	print OUT $map_orphan_q{$set_q}->{$key}, "\n";
    }
    close(OUT);
}

foreach my $set_t (keys %map_orphan_t) {
    my $ofn = $odir."/NONE_".$set_t;
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


