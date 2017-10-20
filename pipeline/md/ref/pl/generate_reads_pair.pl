#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
    print "usage: $0 <coord table> <build fasta> <ofn>\n";
    exit 1;
}

my $cifn = $ARGV[0];
my $bifn = $ARGV[1];
my $ofn1 = $ARGV[2];
my $ofn2 = $ARGV[3];

my %nts = (
    "A" => 1,
    "C" => 1,
    "G" => 1,
    "T" => 1);

my %rc = (
    "A" => "T",
    "C" => "G",
    "G" => "C",
    "T" => "A",
    "N"=> "N");

####################################################################################################
# build
####################################################################################################

print "reading build file: ", $bifn, "\n";

my %contigs;
my $seq = "";
my $contig = "";
open(IN, $bifn) or die $bifn;
while (my $line = <IN>) {
    chomp($line);
    if (substr($line, 0, 1) ne ">") {
	$seq .= $line;
    } else {
 	if ($seq ne "") {
	    $contigs{$contig} = {};
	    $contigs{$contig}->{seq} = $seq;
	    $contigs{$contig}->{length} = length($seq);
	}
	my @f = split(" ", substr($line,1));
	$contig = $f[0];
	$seq = "";
    }
}
close(IN);
if ($seq ne "") {
    $contigs{$contig} = {};
    $contigs{$contig}->{seq} = $seq;
    $contigs{$contig}->{length} = length($seq);
}

foreach my $contig (keys %contigs) {
    my $seq = $contigs{$contig}->{seq};
    $seq =~ s/[^ACGT]/N/g;
    my $length = length($seq);
    my $f = 4 + int(2000 / $length);
    $contigs{$contig}->{ext_seq} = $seq x $f;
}

####################################################################################################
# coords
####################################################################################################

print "reading coords table: ", $cifn, "\n";
open(IN, $cifn) or die $cifn;
my $header = <IN>;
my %h = parse_header($header);

print "generating file: ", $ofn1, "\n";
open(OUT1, ">", $ofn1) or die $ofn1;

print "generating file: ", $ofn2, "\n";
open(OUT2, ">", $ofn2) or die $ofn2;

my $count = 1;
while (my $line = <IN>) {
    chomp($line);

    print "line: $count\n" if ($count % 1000000 == 0);
    $count++;

    my @f = split("\t", $line);
    my $contig1 = $f[$h{contig1}];
    my $start1 = $f[$h{start1}];
    my $end1 = $f[$h{end1}];
    my $strand1 = $f[$h{strand1}];

    my $contig2 = $f[$h{contig2}];
    my $start2 = $f[$h{start2}];
    my $end2 = $f[$h{end2}];
    my $strand2 = $f[$h{strand2}];

    my $id = $contig1."_".$start1."_".$end1."_".$strand1."_".$contig2."_".$start2."_".$end2."_".$strand2;

    defined($contigs{$contig1}) and defined($contigs{$contig2}) or die;
    my $clength1 = $contigs{$contig1}->{length};
    my $clength2 = $contigs{$contig2}->{length};

    $start1 += $clength1;
    $start2 += $clength2;
    $end1 += $clength1;
    $end2 += $clength2;

    my $length1 = $end1 - $start1;
    my $length2 = $end2 - $start2;
    $length1 > 0 and $length2 > 0 or die;

    my $seq1 = substr($contigs{$contig1}->{ext_seq}, $start1, $length1);
    my $seq2 = substr($contigs{$contig2}->{ext_seq}, $start2, $length2);
    $seq1 = get_rc($seq1) if ($strand1 eq -1);
    $seq2 = get_rc($seq2) if ($strand2 eq -1);

    print OUT1 "@", $id, "\n";
    print OUT1 $seq1, "\n";
    print OUT1 "+", "\n";
    print OUT1 ("I") x $length1, "\n";

    print OUT2 "@", $id, "\n";
    print OUT2 $seq2, "\n";
    print OUT2 "+", "\n";
    print OUT2 ("I") x $length2, "\n";

    $id++;
}
close(IN);
close(OUT1);
close(OUT2);

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

# returns reverse-complement of the argument
sub get_rc {
    my ($str) = @_;
    my $len = length($str);
    my $rcstr = "";
    for (my($i) = $len - 1; $i >= 0; $i--) {
	$rcstr .= $rc{substr($str, $i, 1)};
    }
    return($rcstr);
}
