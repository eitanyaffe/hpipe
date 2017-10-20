#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);

if ($#ARGV == -1) {
    print STDERR "usage: $0 <ifn table> <genebank table> <out dir> <failed log table>\n";
	exit 1;
}

my @pnames = ("ifn", "ifn_genebank", "odir", "ofn_failed");
my %p = parse_parameters(\@pnames, \@ARGV);

#######################################################################################
# read genome table
#######################################################################################

my %genomes;

print "reading table: $p{ifn}\n";
open(IN, $p{ifn}) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line, -1);
    my $id = $f[$h{accession}];
    $genomes{$id} = {};
}
close(IN);


#######################################################################################
# read genebank table
#######################################################################################

print "output failed log: $p{ofn_failed}\n";
open(OUT, ">", $p{ofn_failed}) || die;
print OUT "id\tasm\toriginal_asm\tpath\n";

print "reading genebank table: $p{ifn_genebank}\n";
open(IN, $p{ifn_genebank}) || die;

# header
$header = <IN>;
%h = parse_header($header);

while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line, -1);
    my $id = $f[$h{"# assembly_accession"}];
    next if (!defined($genomes{$id}));
    $genomes{$id}->{path} = $f[$h{ftp_path}];
    $genomes{$id}->{asm} = $f[$h{asm_name}];
}

my $N = scalar(keys %genomes);
print "downloading $N files to dir: $p{odir}\n";
system("rm -rf $p{odir}") == 0 or die;
system("mkdir -p $p{odir}") == 0 or die;

print "total ref genomes: ", scalar(keys(%genomes)), "\n";
foreach my $id (keys %genomes) {
    print ".";
    defined($genomes{$id}->{path}) or die;
    my $path = $genomes{$id}->{path};
    my $asm = $genomes{$id}->{asm};
    my $dir = $p{odir}."/".$id;
    my $log = $p{odir}."/.log_".$id;

    my $oasm = $asm;
    $asm =~ s/\#/_/g;
    $asm =~ s/[ ]+$//g;
    $asm =~ s/^[ ]+//g;
    $asm =~ s/ /_/g;
    $asm =~ s/\(//g;
    $asm =~ s/\)//g;
    $asm =~ s/\//_/g;
    $asm =~ s/,//g;

    my $fn = sprintf("%s/%s_%s_genomic.fna.gz", $path, $id, $asm);

    # print "id: $id, asm: $asm, original asm: $oasm\n";
    my $command = sprintf("wget -P %s %s > %s 2>&1", $dir, $fn, $log);
    # print "command: ", $command, "\n";

    my $rc = system($command);
    # second attempt, adding an extra '_'
    if ($rc != 0) {
	$fn = sprintf("%s/%s_%s__genomic.fna.gz", $path, $id, $asm);
	my $command = sprintf("wget -P %s %s > %s 2>&1", $dir, $fn, $log);
	# print "command: ", $command, "\n";
	$rc = system($command);
	exit($rc) if $rc == 2;
    }
    if ($rc != 0) {
	exit($rc) if $rc == 2;
	print OUT "$id\tasm\t$oasm\t$fn\n";
	next;
    }

    my $afile = $dir."/genome.fasta";
    $command = sprintf("gunzip -c %s/*_genomic.fna.gz > %s", $dir, $afile);
    # print "command: ", $command, "\n";
    system($command) == 0 or die;
}
print " done\n";

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

sub parse_parameters
{
    my ($p_pnames, $p_argv) = @_;
    my @pnames = @{$p_pnames};
    my @argv = @{$p_argv};

    my %p;
    @p{@pnames} = @argv;
    my %po;
    my $N = scalar(@argv);
    my @pi = (1..$N);
    @po{@pi} = @pnames;

    print "=============================================\n";
    foreach my $index (sort {$a <=> $b} keys %po) {
	my $key = $po{$index};
	defined($p{$key}) or die "parameter $key not defined (check if all parameters defined)";
	print $key, ": ", $p{$key}, "\n";
    }
    print "=============================================\n";

    return (%p);
}
