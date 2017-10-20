#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::stat;
use Time::localtime;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <usearch_binary> <query_ifn> <target_ifn> <target_db_file> <nthreads> <evalue> <ofn>\n";
	exit 1;
}

my ($usearch_binary,
    $query_ifn, $target_ifn, $target_db_file,
    $nthreads, $evalue, $ofn) = @ARGV;

###################################################################################################################
# usrearch
###################################################################################################################

# index file
if (-e $target_db_file) {
    print "usearch db file exists, no need to create\n";
} else {

    print "generating usearch db file: $target_db_file\n";
    my $command = "$usearch_binary -makeudb_ublast $target_ifn -output $target_db_file\n";
    print "running usearch ...\n";
    system($command)  == 0 or die;
}

my $command = "$usearch_binary -ublast $query_ifn -db $target_db_file -strand both -threads $nthreads -evalue $evalue -blast6out $ofn";
system($command)  == 0 or die $command;

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
