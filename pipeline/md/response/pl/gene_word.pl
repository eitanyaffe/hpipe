#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <gene table> <gene table desc field> <background gene table> <background gene table desc field> <ofn> <ofn backtable> <whole words> <poor annotation descriptions>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $field = $ARGV[1];
my $ifn_bg = $ARGV[2];
my $field_bg = $ARGV[3];
my $ofn = $ARGV[4];
my $ofn_back = $ARGV[5];
my $whole_list = $ARGV[6];
my $poor_list = $ARGV[7];

our @poor_descs = split(" ", $poor_list);
for (my $i=0; $i<@poor_descs; $i++) {
    $poor_descs[$i] =~ s/_/ /g;
    $poor_descs[$i] = lc($poor_descs[$i]);
}
print "Poor annotations: \n ", join("\n ", @poor_descs), "\n";

our @whole_words = split(" ", $whole_list);
for (my $i=0; $i<@whole_words; $i++) {
    $whole_words[$i] =~ s/_/ /g;
}
print "Whole words: \n ", join("\n ", @whole_words), "\n";

###############################################################################################
# selected genes
###############################################################################################

my %words;
my %genes;

print STDERR "reading table: $ifn\n";
open(IN, $ifn) || die $ifn;
my $header = <IN>;
my %h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $desc = $f[$h{$field}];
    next if (is_poor($desc));

    $genes{$gene} = 1;
    my @lwords = get_words($desc);
    foreach my $word (@lwords) {
	$word = ucfirst($word);
	if (!defined($words{$word})) {
	    $words{$word} = {};
	    $words{$word}->{genes} = {};
	    $words{$word}->{genes_bg} = {};
	}
	$words{$word}->{genes}->{$gene} = $desc;
    }
}
close(IN);

###############################################################################################
# background genes
###############################################################################################

print STDERR "reading table: $ifn_bg\n";
open(IN, $ifn_bg) || die $ifn_bg;
$header = <IN>;
%h = parse_header($header);

my %genes_bg;

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    my $desc = $f[$h{$field}];
    next if (is_poor($desc));
    $genes_bg{$gene} = 1;
    my @lwords = get_words($desc);
    foreach my $word (@lwords) {
	$word = ucfirst($word);
	next if (!defined($words{$word}));
	$words{$word}->{genes_bg}->{$gene} = 1;
    }
}
close(IN);

###############################################################################################
# compute ratios
###############################################################################################

my $total_gene_count = scalar keys %genes;
my $total_gene_count_bg = scalar keys %genes_bg;
foreach my $word (keys %words) {
    my $word_gene_count = scalar keys %{$words{$word}->{genes}};
    my $word_gene_count_bg = scalar keys %{$words{$word}->{genes_bg}};
    $words{$word}->{gene_count} = $word_gene_count;
    $words{$word}->{gene_percent} = 100 * $word_gene_count / $total_gene_count;
    $words{$word}->{gene_count_bg} = $word_gene_count_bg;
    $words{$word}->{gene_percent_bg} = 100 * $word_gene_count_bg / $total_gene_count_bg;
    $words{$word}->{enrichment} = $words{$word}->{gene_percent} / $words{$word}->{gene_percent_bg};
}


###############################################################################################
# output
###############################################################################################

print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "word\tenrichment\tgene_count\tgene_percent\tgene_count_bg\tgene_percent_bg\n";

foreach my $word (sort { $words{$b}->{enrichment} <=> $words{$a}->{enrichment} } keys %words) {
    print OUT $word, "\t", round($words{$word}->{enrichment},4), "\t", $words{$word}->{gene_count}, "\t", round($words{$word}->{gene_percent},8), "\t";
    print OUT $words{$word}->{gene_count_bg}, "\t", round($words{$word}->{gene_percent_bg},8), "\n";
}

close(OUT);

###############################################################################################
# output back table
###############################################################################################

print STDERR "generating back table file: $ofn_back\n";
open(OUT, ">", $ofn_back) || die $ofn_back;
print OUT "word\tenrichment\tgene_count\tgene\tdesc\n";

foreach my $word (sort { $words{$b}->{gene_count} <=> $words{$a}->{gene_count} } keys %words) {
    foreach my $gene (sort keys %{$words{$word}->{genes}}) {
	my $desc = $words{$word}->{genes}->{$gene};
	print OUT $word, "\t", round($words{$word}->{enrichment},4), "\t", $words{$word}->{gene_count}, "\t", $gene, "\t", $desc, "\n";
    }
}
close(OUT);

######################################################################################################
# Subroutines
######################################################################################################

sub is_poor
{
    my ($desc) = @_;
    my $result = $desc eq "NA";
    $desc = lc($desc);
    for (my $i=0; $i<@poor_descs; $i++) {
	$result = 1 if ($desc eq $poor_descs[$i]);
    }
    return ($result);
}

sub fix_desc
{
    my ($desc) = @_;
    $desc =~ s/,/ /g;
    $desc =~ s/:/ /g;
    $desc =~ s/\(/ /g;
    $desc =~ s/\// /g;
    $desc =~ s/\)/ /g;
    $desc =~ s/-/ /g;
    return ($desc);
}

sub round
{
    my ($value,$digits) = @_;
    return (floor($value*(10**$digits))/(10**$digits));
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

sub get_words
{
    my ($desc) = @_;
    my @lwords;
    foreach my $word (@whole_words) {
	my $index = index(lc($desc), lc($word));
	if ($index != -1) {
	    my $length = length($word);
	    push(@lwords, $word);
	    # print "found: ", $word, "\n";
	    # print "prior: ", $desc, "\n";
	    substr($desc, $index, $length) = " " x $length;
	    # print "post: ", $desc, "\n";
	}
    }
    $desc = fix_desc($desc);

    push(@lwords, split(" ", $desc));
    return (@lwords);
}
