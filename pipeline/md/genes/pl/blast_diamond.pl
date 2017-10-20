#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::stat;
use Time::localtime;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <diamond_binary> <query_ifn> <target_ifn> <target_db_file> <nthreads> <evalue> <ofn>\n";
	exit 1;
}

my ($binary,
    $query_ifn, $target_ifn, $target_db_file,
    $nthreads, $evalue, $ofn) = @ARGV;

###################################################################################################################
# usrearch
###################################################################################################################

# index file
if (-e $target_db_file) {
    print "db file exists, no need to create\n";
} else {

    print "generating db file: $target_db_file\n";
    my $command = "$binary makedb -c 1 -b 20 --in $target_ifn -p $nthreads -d $target_db_file\n";
    print "diamond index: $command\n";
    # system($command)  == 0 or die;
}

my $command = "$binary blastx -c 1 -d $target_db_file --in $query_ifn -p $nthreads -e $evalue -a $ofn";
print "diamond blast: $command\n";
#system($command)  == 0 or die $command;

$command = "$binary view -c 1 -a $ofn -f sam -o $ofn";
print "diamond view: $command\n";
#system($command)  == 0 or die $command;

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
