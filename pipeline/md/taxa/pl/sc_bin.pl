#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);
use POSIX qw( strftime );
use File::stat;
use File::Basename;

if ($#ARGV == -1) {
    print STDERR "usage: $0 <type dir|file> <input path> <source extension (for dir only)> <format fasta|fastq> <reads per split file> <length> <binsize> <odir>\n";
	exit 1;
}

my @pnames = (
    "type",
    "input",
    "ext",
    "format",
    "length",
    "binsize",
    "reads_per_split_file",
    "odir"
    );
my %p;
@p{@pnames} = @ARGV;

print "=============================================\n";
foreach my $key (keys %p) {
    defined($p{$key}) or die "parameter $key not defined (check if all parameters defined)";
    print $key, ": ", $p{$key}, "\n";
}
print "=============================================\n";

$p{type} eq "dir" or $p{type} eq "file" or die "unknown type: ".$p{input};
$p{format} eq "fastq" or $p{format} eq "fasta" or die "unknown format: ".$p{format};

#######################################################################################
# prepare intermediate directories
#######################################################################################

my $odir = $p{odir};
my @args = ("rm", "-rf", $odir);
system(@args) == 0 or die "system @args failed: $?";

@args = ("mkdir", "-p", $odir);
system(@args) == 0 or die "system @args failed: $?";

#######################################################################################
# split input sequence
#######################################################################################

my $read_count = 0;

# index of split file
my $split_index = 1;

# index of output read segment
my $split_count = 0;

my @split_fns;
my $split_fn = $odir."/".$split_index;
push(@split_fns, $split_fn);
print "generating file: $split_fn\n";
open(OUT, ">", $split_fn) || die;

my @ifns;

if ($p{type} eq "file") {
    push(@ifns, $p{input});
} else {
    push(@ifns, <$p{input}/*.$p{ext}>);
}
print "number of input files: ", scalar(@ifns), "\n";
for my $ifn (@ifns) {

    my $line_count = 0;
    my $seq = "";
    my $label = "";
    print "proccessing input file: $ifn\n";
    open(IN, $ifn) || die;
    while (my $line = <IN>) {
	chomp($line);
	if ($p{format} eq "fastq") {
	    if ($line_count % 4 == 0) {
		$label = "I";
	    } elsif ($line_count % 4 == 1) {
		$seq = $line;
	    }
	} else {
	    if ($line_count % 2 == 0) {
		$label = substr($line, 1);
	    } elsif ($line_count % 2 == 1) {
		$seq = $line;
	    }
	}

	if ($seq ne "") {
	    $read_count++;
	    for (my $i = 0; $i<length($seq); $i+=$p{binsize}) {
		my $sseq = substr($seq, $i, $p{length});
		next if (length($sseq) != $p{length});
		$split_count++;
		if ($split_count % $p{reads_per_split_file} == 0) {
		    close(OUT);
		    $split_index++;
		    $split_fn = $odir."/".$split_index;
		    push(@split_fns, $split_fn);
		    print "generating file: $split_fn\n";
		    open(OUT, ">", $split_fn) || die;
		}
		my $slabel = $label." read_index:".$read_count." read_offset:".$i;
		$slabel =~ s/ /_/g;
		print OUT ">".$slabel."\n";
		print OUT $sseq, "\n";
	    }
	    $seq = "";
	    $label = "";
	}
	$line_count++;
    }
    close(IN);
}
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


