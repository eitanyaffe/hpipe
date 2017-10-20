#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);

if ($#ARGV == -1) {
    print STDERR "usage: $0 <ifn table> <genebank table> <out dir>\n";
	exit 1;
}

my @pnames = ("ifn", "ifn_genebank", "odir");
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

foreach my $id (keys %genomes) {
    print ".";
    defined($genomes{$id}->{path}) or die;
    my $path = $genomes{$id}->{path};
    my $asm = $genomes{$id}->{asm};
    my $dir = $p{odir}."/".$id;
    my $log = $p{odir}."/.log_".$id;

    $asm =~ s/\#/_/g;
    $asm =~ s/ /_/g;
    # print "id: $id, asm: $asm\n";
    my $command = sprintf("wget -P %s %s/%s_%s_genomic.fna.gz > %s 2>&1", $dir, $path, $id, $asm, $log);
    print "command: ", $command, "\n";
    system($command) == 0 or die;
    my $afile = $dir."/genome.fasta";
    $command = sprintf("gunzip -c %s/*_genomic.fna.gz > %s", $dir, $afile);
    # print "command: ", $command, "\n";
    system($command) == 0 or die;
}
print " done\n";

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
