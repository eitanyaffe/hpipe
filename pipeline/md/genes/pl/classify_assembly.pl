#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <contig_ifn> <ref_gene_ifn> <pre_gene_ifn> <map_ifn> <id_threshold_exact> <id_threshold_weak> <ofn ref genes> <ofn genes> <ofn contigs>\n";
	exit 1;
}

my ($contig_ifn, $ref_gene_ifn, $pre_gene_ifn, $map_ifn, $id_threshold_strong, $id_threshold_weak, $ofn_ref_gene, $ofn_gene, $ofn_contig) = @ARGV;

###################################################################################################################
# ref genes
###################################################################################################################

my %ref_genes;

print STDERR "load ref genes: $ref_gene_ifn\n";
open(IN, $ref_gene_ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $ref = $f[$h{contig}];

    $ref_genes{$gene} = $ref;
}

###################################################################################################################
# contigs
###################################################################################################################

my @fields = ("total", "match", "match.weak", "spurious", "no_ref", "switch");

my %contigs;

print STDERR "load contigs: $contig_ifn\n";
open(IN, $contig_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    $contigs{$contig} = {};
    $contigs{$contig}->{genes} = {};
    $contigs{$contig}->{ref} = "none";

    $contigs{$contig}->{counts} = {};
    foreach my $field (@fields) {
	$contigs{$contig}->{counts}->{$field} = 0;
    }
}

###################################################################################################################
# predicted genes
###################################################################################################################

my %pre_genes;

print STDERR "load pre genes: $pre_gene_ifn\n";
open(IN, $pre_gene_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $coord = $f[$h{start}];
    my $contig = $f[$h{contig}];

    $pre_genes{$gene} = {};
    $pre_genes{$gene}->{identity} = 0;
    $pre_genes{$gene}->{coord} = $coord;

    defined($contigs{$contig}) or die;

    $pre_genes{$gene}->{contig} = $contig;
    $contigs{$contig}->{genes}->{$coord} = $gene;

}

###################################################################################################################
# ref to pre map
###################################################################################################################

my %map_fwd;
my %map_bck;

print STDERR "load ref2pre map: $map_ifn\n";
open(IN, $map_ifn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $identity = $f[$h{identity}];
    my $rgene = $f[$h{gene1}];
    my $pgene = $f[$h{gene2}];

    defined($ref_genes{$rgene}) or die;

    $map_fwd{$pgene} = {} if (!defined($map_fwd{$rgene}));
    $map_fwd{$rgene}->{$pgene} = $identity;

    $map_bck{$pgene} = {} if (!defined($map_bck{$pgene}));
    $map_bck{$pgene}->{$rgene} = $identity;
}

###################################################################################################################
# classify contigs
###################################################################################################################

my $debug_contig = "c15537";
foreach my $contig (keys %contigs)
{
    my %refs;
    my $size = scalar(keys %{$contigs{$contig}->{genes}});

    foreach my $coord (sort {$a <=> $b} keys %{$contigs{$contig}->{genes}}) {
	my $pgene = $contigs{$contig}->{genes}->{$coord};
	defined($pre_genes{$pgene}) or die;
	next if (!defined($map_bck{$pgene}));

	my @max_rgenes;
	my $max_identity = 0;
	foreach my $rgene (keys %{$map_bck{$pgene}}) {
	    my $identity = $map_bck{$pgene}->{$rgene};
	    if ($max_identity < $identity) {
		$max_identity = $identity;
		@max_rgenes = ();
	    }
	    push(@max_rgenes, $rgene) if ($max_identity <= $identity);
	}
	next if ($max_identity < $id_threshold_weak);
	my $weight = 1 / scalar(@max_rgenes);
	foreach my $max_rgene (@max_rgenes) {
	    my $ref = $ref_genes{$max_rgene};
	    print $max_rgene,  " : ", $ref, ", ", $weight, ",", $max_identity, "\n" if ($contig eq $debug_contig);
	    $refs{$ref} = 0 if (!defined($refs{$ref}));
	    $refs{$ref} += $weight;
	}
    }
    my $max_ref = "";
    my $max_count = 0;
    foreach my $ref (keys %refs) {
	my $count = $refs{$ref};
	print $ref, " : ", $count, "\n" if ($contig eq $debug_contig);
	if ($max_count < $count) {
	    $max_count = $count;
	    $max_ref = $ref;
	}
    }

    # mark contig if popular ref supported by >1 gene or a singlton gene on a contig
    $contigs{$contig}->{ref} = $max_ref if ($max_count > 1 || ($size == 1 && $max_count == 1));

    print $contig, " : ", $contigs{$contig}->{ref}, "\n" if ($contig eq $debug_contig);
}

###################################################################################################################
# classify/output pre genes
###################################################################################################################

print STDERR "generating gene output table: $ofn_gene\n";
open(OUT, ">", $ofn_gene) || die;
print OUT "gene\tcontig\tref\ttype\tbest_ref\tbest_rgene\tbest_identity\n";

foreach my $pgene (keys %pre_genes) {
    my $contig = $pre_genes{$pgene}->{contig};
    defined($contigs{$contig}) or die;

    $contigs{$contig}->{counts}->{total}++;
    my $ref = $contigs{$contig}->{ref};

    # has a good match on selected ref?
    my $ref_identity = 0;
    my $max_identity = 0;
    my $max_ref = "none";
    my $max_rgene = "none";
    foreach my $rgene (keys %{$map_bck{$pgene}}) {
	my $identity = $map_bck{$pgene}->{$rgene};
	my $oref = $ref_genes{$rgene};
	$ref_identity = $identity if (($ref eq $oref) && ($ref_identity < $identity));
	if ($max_identity < $identity || ($max_identity == $identity && $oref eq $ref)) {
	    $max_identity = $identity;
	    $max_ref = $oref;
	    $max_rgene = $rgene;
	}
    }

    my $type;
    if ($ref eq "none") {
	$type = "no_ref";
    } elsif ($ref_identity > $id_threshold_strong) {
	$type = "match";
    } elsif ($ref_identity > $id_threshold_weak) {
	$type = "match.weak";
    } elsif ($max_identity < $id_threshold_strong) {
	$type = "spurious"
    } else {
	$type = "switch";
    }

    $contigs{$contig}->{counts}->{$type}++;
    print OUT $pgene, "\t", $contig, "\t", $ref, "\t", $type, "\t", $max_ref, "\t", $max_rgene, "\t", $max_identity, "\n";
}
close(OUT);

###################################################################################################################
# output contigs
###################################################################################################################

print STDERR "generating contig output table: $ofn_contig\n";
open(OUT, ">", $ofn_contig) || die;
print OUT "contig\tref";
foreach my $field (@fields) {
    print OUT "\t$field";
}
print OUT "\ttype\n";

foreach my $contig (keys %contigs) {
    my $ref = $contigs{$contig}->{ref};

    next if ($contigs{$contig}->{counts}->{total} == 0);

    print OUT "$contig", "\t", $ref;
    foreach my $field (@fields) {
	print OUT "\t", $contigs{$contig}->{counts}->{$field};
    }
    print OUT "\t", ($ref eq "none") ? "no_ref" : (($contigs{$contig}->{counts}->{switch} > 0) ? "chimeric" : "normal"), "\n";
}

###################################################################################################################
# output classified ref genes
###################################################################################################################

print STDERR "generating ref gene table: $ofn_ref_gene\n";
open(OUT, ">", $ofn_ref_gene) || die;
print OUT "gene\tref\ttype\n";

foreach my $rgene (keys %ref_genes) {
    my $ref = $ref_genes{$rgene};

    my $max_identity = 0;
    my $type;

    if (defined($map_fwd{$rgene})) {
	# single best hit
	foreach my $pgene (keys %{$map_fwd{$rgene}}) {
	    my $identity = $map_fwd{$rgene}->{$pgene};
	    if ($max_identity < $identity) {
		$max_identity = $identity;
	    }
	}
	if ($max_identity < $id_threshold_weak) {
	    $type = "lost";
	} elsif ($max_identity < $id_threshold_strong) {
	    $type = "weak";
	} else {
	    $type = "found";
	}
    } else {
	$type = "lost";
    }
    $type ne "" or die;
    print OUT $rgene, "\t", $ref, "\t", $type, "\n";
}

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

