#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <uniref table> <uniref table> <uniprot2go table> <uniparc2go table> <use uniparc> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ifn_uniref = $ARGV[1];
my $ifn_uniprot2go = $ARGV[2];
my $ifn_uniparc2uniprot = $ARGV[3];
my $use_uniparc = $ARGV[4];
my $ofn = $ARGV[5];

###############################################################################################
# read genes to memory
###############################################################################################

my %genes;
my %uniref2genes;

print STDERR "reading table: $ifn\n";
open(IN, $ifn) || die $ifn;
my $header = <IN>;
my %h = parse_header($header);
chomp($header);
my $out_header = $header."\tGO";

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $uniref = $f[$h{uniref}];

    $genes{$gene} = {};
    $genes{$gene}->{line} = $line;
    $genes{$gene}->{uniref} = $uniref;
    $genes{$gene}->{uniprots} = {};
    $genes{$gene}->{uniparcs} = {};
    $genes{$gene}->{GOs} = {};

    $uniref2genes{$uniref} = {} if (!defined($uniref2genes{$uniref}));
    $uniref2genes{$uniref}->{$gene} = 1;
}
close(IN);

###############################################################################################
# uniref to uniprot/uniparc
###############################################################################################

print STDERR "reading table: $ifn_uniref\n";
open(IN, $ifn_uniref) || die $ifn_uniref;
$header = <IN>;
%h = parse_header($header);

my %uniprot2genes;
my %uniparc2genes;

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $uniref = $f[$h{uniref}];
    my $uniprots = $f[$h{uniprot}];
    my $uniparcs = $f[$h{uniparc}];

    next if (!defined($uniref2genes{$uniref}));

    for my $gene (keys %{$uniref2genes{$uniref}}) {
	for my $uniprot (split(";", $uniprots)) {
	    $genes{$gene}->{uniprots}->{$uniprot} = 1;
	    $uniprot2genes{$uniprot} = {} if (!defined($uniprot2genes{$uniprot}));
	    $uniprot2genes{$uniprot}->{$gene} = 1;
	}
	for my $uniparc (split(";", $uniparcs)) {
	    $genes{$gene}->{uniparcs}->{$uniparc} = 1;
	    $uniparc2genes{$uniparc} = {} if (!defined($uniparc2genes{$uniparc}));
	    $uniparc2genes{$uniparc}->{$gene} = 1;
	}
    }
}
close(IN);

###############################################################################################
# uniparc2uniprot lookup
###############################################################################################

if ($use_uniparc eq "T") {
    print STDERR "reading table: $ifn_uniparc2uniprot\n";
    open(IN, $ifn_uniparc2uniprot) || die $ifn_uniparc2uniprot;
    $header = <IN>;
    %h = parse_header($header);

    while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $uniparc = $f[$h{uniparc}];
	my $uniprots = $f[$h{uniprots}];
	next if (!defined($uniparc2genes{$uniparc}) || $uniprots eq "NA");
	for my $gene (keys %{$uniparc2genes{$uniparc}}) {
	    for my $uniprot (split(";", $uniprots)) {
		$genes{$gene}->{uniprots}->{$uniprot} = 1;
		$uniprot2genes{$uniprot} = {} if (!defined($uniprot2genes{$uniprot}));
		$uniprot2genes{$uniprot}->{$gene} = 1;
	    }
	}
    }
    close(IN);
}

###############################################################################################
# uniprot2GO lookup
###############################################################################################

print STDERR "reading table: $ifn_uniprot2go\n";
open(IN, $ifn_uniprot2go) || die $ifn_uniprot2go;
$header = <IN>;
%h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $uniprot = $f[$h{uniprot}];
    my $GOs = $f[$h{GO}];

    if (defined($uniprot2genes{$uniprot}) && $GOs ne "NA") {
	for my $GO (split(";", $GOs)) {
	    for my $gene (keys %{$uniprot2genes{$uniprot}}) {
		$genes{$gene}->{GOs}->{$GO} = 1;
	    }
	}
    }
}
close(IN);

###############################################################################################
# output genes
###############################################################################################

print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT $out_header, "\n";
for my $gene (keys %genes) {
    my $GOs = scalar(keys %{$genes{$gene}->{GOs}}) > 0 ? join(";", keys %{$genes{$gene}->{GOs}}) : "NA";
    print OUT $genes{$gene}->{line}, "\t", $GOs, "\n";
}
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
