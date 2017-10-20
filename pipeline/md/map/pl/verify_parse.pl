#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <ref fasta>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ref = $ARGV[1];


######################################################################################################
# to validate we start with the ref genome and transform it according to subs/inserts/deletes
# and then verify it equals the  read sequence
######################################################################################################


######################################################################################################
# read contig fastq into memory
######################################################################################################

print STDERR "reading ref into memory: $ref\n";
open(IN, $ref) || die $ref;


my %href;
my $seq = "";
my $contig = "";
while (my $line = <IN>) {
    chomp($line);
    if (substr($line, 0, 1) ne ">") {
	$seq .= $line;
    } else {
	$href{$contig} = $seq;
	my @f = split(" ", substr($line,1));
	$contig = $f[0];
	$seq = "";
    }
}
$href{$contig} = $seq if ($contig ne "");

######################################################################################################
# read contig fastq into memory
######################################################################################################

open(IN, $ifn) || die;
my $header = <IN>;
my %h = parse_header($header);

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
#    print $line, "\n";
    my $id = $f[$h{id}];
    my $contig = $f[$h{contig}];
    my $strand = $f[$h{strand}];
    my $start_coord = ($strand == 1 ? $f[$h{back_coord}] : $f[$h{coord}]);
    my $end_coord = ($strand == 1 ? $f[$h{coord}] : $f[$h{back_coord}])+1;
    my $length = $end_coord - $start_coord + 1;
    my $seq = $f[$h{sequence}];

    my ($clip_start, $clip_end) = split(";", $f[$h{clip}]);

    defined($href{$contig}) or die;
    my $ref = $href{$contig};

    my $v = '.' x length($seq);

    # apply subs
    my @ff = split(";", substr($f[$h{substitute}],1));
    foreach my $sfield (@ff) {
	my ($ref_coord, $seq_coord, $ref_nt, $seq_nt) = split(",", $sfield);

	substr($seq, $seq_coord-1, 1) eq $seq_nt or die "seq=$seq coord=$seq_coord nt=".substr($seq, $seq_coord-1, 1)." expected=$seq_nt";
	substr($ref, $ref_coord-1, 1) eq $ref_nt or die "ref=$ref coord=$ref_coord nt=".substr($ref, $ref_coord-1, 1)." expected=$ref_nt";

	substr($ref, $ref_coord-1, 1) = $seq_nt;
	substr($v, $seq_coord-1, 1) = "S";
    }

    # apply delete
    @ff = split(";", substr($f[$h{delete}],1));
    foreach my $sfield (@ff) {
	my ($ref_coord, $length) = split(",", $sfield);
	my $nts = '*' x $length;
	substr($ref, $ref_coord-1, $length) = $nts;
    }

    # collect inserts into hash
    my %hinsert;
    @ff = split(";", substr($f[$h{insert}],1));
    foreach my $sfield (@ff) {
	my ($ref_coord, $seq_coord, $nts) = split(",", $sfield);
	my $length = length($nts);
	substr($seq, $seq_coord-1, $length) eq $nts or die "nts=$nts found=".substr($seq, $seq_coord-1, $length);
	substr($v, $seq_coord-1, $length) = 'I' x $length;

	$hinsert{$ref_coord} = $nts;
    }

    # add insert sequences
    my $result = "";
    my @insert_coords = sort {$a <=> $b} keys %hinsert;
    if (scalar(@insert_coords) > 0) {
	my $pcoord = $start_coord;
	foreach my $icoord (@insert_coords) {
	    $result .= substr($ref, $pcoord-1, $icoord-$pcoord).$hinsert{$icoord};
	    $pcoord = $icoord;
	}
	$result .= substr($ref, $pcoord-1, $end_coord-$pcoord);
    } else {
	$result = substr($ref, $start_coord-1, $end_coord-$start_coord);
    }

    # clean deleted sequences
    $result =~ s/\*//g;

    $seq = substr($seq, $clip_start-1, $clip_end-$clip_start+1);
    $v = substr($v, $clip_start-1, $clip_end-$clip_start+1);
    $result eq $seq or die "$line\nread:".$seq."\n    :".$v."\n ref:".$result;
    #print "read:".$seq."\n    :".$v."\n ref:".$result."\n\n";
}
close(IN);

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
