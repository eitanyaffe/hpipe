#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <ofn> <field1 field2....>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ofn = $ARGV[1];
shift; shift;

my @fields = @ARGV;
print "binning fields: ", join(",",@fields), "\n";

#############################################################################################
# init
#############################################################################################

my %ht;
foreach my $field (@fields) {
    !defined($ht{$field}) or die;
    $ht{$field} = {};
    $ht{$field}->{bin} = 1;
    $ht{$field}->{field2bin} = {};
    $ht{$field}->{bin2field} = {};
}

#############################################################################################
# go over table
#############################################################################################

open(IN, $ifn) || die;
print STDERR "Reading file $ifn into hash...\n";
my $header = <IN>;
my %h = parse_header($header);

print STDERR "Writing output file: $ofn\n";
open(OUT, ">", $ofn) || die;
chomp($header);
print OUT $header;
foreach my $field (@fields) {
    print OUT "\t".$field."_bin";
}
print OUT "\n";

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);

    print OUT $line;
    foreach my $field (@fields) {
	my $value = $f[$h{$field}];
	my $bin = $ht{$field}->{bin};
	if (!defined($ht{$field}->{field2bin}->{$value})) {
	    $ht{$field}->{field2bin}->{$value} = $bin;
	    $ht{$field}->{bin2field}->{$bin} = $value;
	    $ht{$field}->{bin}++;
	}
	print OUT "\t", $ht{$field}->{field2bin}->{$value};
    }
    print OUT "\n";
}
close(IN);
close(OUT);


foreach my $field (@fields) {
    my $ofn_bins = $ofn.".".$field;
    print STDERR "Writing output bin file: $ofn_bins\n";
    open(OUT, ">", $ofn_bins) || die;
    print OUT "bin\t$field\n";
    foreach my $bin (sort {$a <=> $b} keys %{$ht{$field}->{bin2field}}) {
	my $value = $ht{$field}->{bin2field}->{$bin};
	print OUT "$bin\t$value\n";
    }
    close(OUT);
}

######################################################################################################
# Subroutines
######################################################################################################

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
