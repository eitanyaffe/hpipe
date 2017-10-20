#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);
use POSIX qw( strftime );
use File::stat;
use File::Basename;

if ($#ARGV == -1) {
    print STDERR "usage: $0 <type dir|file> <input path> <source extension (for dir only)> <format fasta|fastq> <reads per split file> <length> <step> <max reads> <odir>\n";
	exit 1;
}

my @pnames = (
    "type",
    "input",
    "ext",
    "format",
    "length",
    "step",
    "reads_per_split_file",
    "max_reads",
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
print "clearing output directory: $odir\n";

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
    my $contig = "";
    print "proccessing input file: $ifn\n";
    open(IN, $ifn) || die;
    while (my $line = <IN>) {
	chomp($line);
	if ($p{format} eq "fastq") {
	    save_seq($line, "read:".$read_count++) if ($line_count % 4 == 1);
	} else {
	    if (substr($line, 0, 1) eq ">") {
		save_seq($seq, $contig) if ($seq ne "");
		$contig = substr($line, 1);
		my $ix = index($contig," ");
		$contig = substr($contig,0,$ix) if ($ix != -1);
		$seq = "";
	    } else {
		$seq .= $line;
	    }
	}
	$line_count++;

	last if ($p{max_reads} && $read_count > $p{max_reads});
    }
    close(IN);
    if ($p{format} eq "fastq") {
	save_seq($seq, $read_count++);
    } else {
	save_seq($seq, $contig) if ($seq ne "");
    }

    last if ($p{max_reads} && $read_count > $p{max_reads});
}
close(OUT);

#######################################################################################
# utils
#######################################################################################

sub save_seq
{
    my ($seq, $prefix) = @_;
    for (my $i = 0; $i<length($seq); $i+=$p{step}) {
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
	my $label;
	$label = "@".$prefix."_start:".($i+1)."_end:".($i+$p{length});
	print OUT $label."\n";
	print OUT $sseq, "\n";
    }
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


