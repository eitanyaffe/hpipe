#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);
use POSIX qw( strftime );
use File::stat;
use File::Basename;

if ($#ARGV == -1) {
    print STDERR "usage: $0 <idir> <input target fasta> <bwa binary> <parse bwa script> <profile binary> <source contig table> <target contig table> <read length> >parse source> <nthreads> <bwa index prefix> <odir>\n";
	exit 1;
}

my @pnames = (
    "idir",
    "target_fasta",
    "bwa_binary",
    "parse_bwa",
    "profile_binary",
    "src_contig_table",
    "tgt_contig_table",
    "read_length",
    "parse_source",
    "nthreads",
    "bwa_index_prefix",
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

#######################################################################################
# prepare mapped directories
#######################################################################################

my $odir = $p{odir};

my @args = ("rm", "-rf", $odir);
system(@args) == 0 or die "system @args failed: $?";

my $raw_dir = $odir."/raw";
my $parsed_dir = $odir."/parsed";
my $stat_dir = $odir."/stats";

@args = ("mkdir", "-p", $raw_dir, $parsed_dir, $stat_dir);
system(@args) == 0 or die "system @args failed: $?";

#######################################################################################
# map with bwa
#######################################################################################

my @split_fns = <$p{idir}/*>;
my @raw_fns;

print "number of input files: ", scalar(@split_fns), "\n";
foreach my $sfn (@split_fns) {
    my $mfn = $raw_dir."/".basename($sfn);
    push(@raw_fns, $mfn);
    print "=============================================\n";
    @args = ($p{bwa_binary}, "mem", "-T", 0, "-t", $p{nthreads}, $p{bwa_index_prefix}, $sfn, ">", $mfn);
    my $command = join(" ", @args);
    print "running command: $command\n";
    system($command) == 0 or die;
    print "=============================================\n";
}

#######################################################################################
# parse results
#######################################################################################

foreach my $fn (@raw_fns) {
    my $pfn = $parsed_dir."/".basename($fn);
    my $stat_fn = $stat_dir."/".basename($fn);
    @args = ($p{parse_bwa}, $fn, $pfn, $stat_fn);
    my $command = join(" ", @args);
    # print "running command: $command\n";
    system($command) == 0 or die;
}

#######################################################################################
# create summary table
#######################################################################################

my $command;
if ($p{parse_source} eq "T") {
    $command = sprintf("%s -parse_src T -idir %s -src_contig_table %s -tgt_contig_table %s -read_length %d -src_ofn %s -tgt_ofn %s -src_summary_ofn %s -tgt_summary_ofn %s",
		       $p{profile_binary},
		       $parsed_dir,
		       $p{src_contig_table},
		       $p{tgt_contig_table},
		       $p{read_length},
		       $odir."/src_table",
		       $odir."/tgt_table",
		       $odir."/src_summary",
		       $odir."/tgt_summary");
} else {
    $command = sprintf("%s -parse_src F -idir %s -tgt_contig_table %s -read_length %d -tgt_ofn %s -tgt_summary_ofn %s",
		       $p{profile_binary},
		       $parsed_dir,
		       $p{tgt_contig_table},
		       $p{read_length},
		       $odir."/tgt_table",
		       $odir."/tgt_summary");
}
print "running summary command: $command\n";
system($command) == 0 or die;

#######################################################################################
# remove intermediate files
#######################################################################################

$command = sprintf("rm -rf %s %s", $parsed_dir, $raw_dir);
print "clearing intermediate files: $command\n";
system($command) == 0 or die;

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


