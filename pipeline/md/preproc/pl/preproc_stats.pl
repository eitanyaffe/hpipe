#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
    print STDERR "usage: $0 <input dir> <quality trimmed dir> <adaptor dir> <no human dir> <ofn>\n";
    exit 1;
}

my $input_dir = $ARGV[0];
my $trimmed_dir = $ARGV[1];
my $adaptor_dir = $ARGV[2];
my $human_dir = $ARGV[3];
my $ofn = $ARGV[4];

print "generating file: $ofn\n";
open(OUT, ">", $ofn);
print OUT "field\tcount\n";

line_count($input_dir."/R1", "input_r1");
line_count($input_dir."/R2", "input_r2");

line_count($trimmed_dir."/R1.fastq", "quality_trimmed_r1");
line_count($trimmed_dir."/R2.fastq", "quality_trimmed_r2");

line_count($adaptor_dir."/R1_*.fastq", "no_adaptor_r1");
line_count($adaptor_dir."/R2_*.fastq", "no_adaptor_r2");

line_count($human_dir."/R1_*.fastq_clean.fq", "no_human_r1");
line_count($human_dir."/R2_*.fastq_clean.fq", "no_human_r2");

close(OUT);

sub line_count
{
    my ($ifn, $title) = @_;
    print "counting lines for files matching pattern $ifn...\n";
    my @lines = `wc -l $ifn`;
    my $line = $lines[-1];
    my @f = split(/\s+/, $line);
    my $lc = $f[scalar(@lines) == 1 ? 0 : 1];
    print OUT $title, "\t", $lc, "\n";
    print $title, ": ", $lc, "\n";
}

