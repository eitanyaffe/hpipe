#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <input dir> <contig table> <similarity mask table> <offset> <output>\n";
	exit 1;
}

my ($idir, $icontig, $icoverage, $isim, $offset, $odir_contact, $odir_masked, $ofn) = @ARGV;

#######################################################################################
# read contig table
#######################################################################################

# contig length table
my %contigs;

print STDERR "reading contig table: $icontig\n";
open(IN, $icontig) || die $icontig;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $length = $f[$h{length}];
    $contigs{$contig} = {};
    $contigs{$contig}->{length} = $length;
}
close(IN);

#######################################################################################
# read coverage table
#######################################################################################

print STDERR "reading contig coverage table: $icoverage\n";
open(IN, $icoverage) || die $icoverage;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);

    my $contig = $f[$h{contig}];
    my $abundance = $f[$h{"relative.abundance"}];
    defined($contigs{$contig}) or die;
    $contigs{$contig}->{abundance} = $abundance;
}
close(IN);

#######################################################################################
# read similarity table
#######################################################################################

# similarity double hashtable
my %sim;

my @isims = <$isim/*table>;
scalar(@isims) > 0 or die "no similarity tables found in dir: $isim\n";

print STDERR "found ".scalar(@isims)." files in similarity dir: $isim\n";
foreach my $ifn (@isims) {
#    print STDERR "traversing file: $ifn\n";
    open(IN, $ifn) || die $ifn;
    my $header = <IN>;
    my %h = parse_header($header);
    while (my $line = <IN>) {
	chomp($line);
	my @f = split("\t", $line);

	# contig1
	my $contig1 = $f[$h{contig1}];
	my $start1 = $f[$h{start1}] - $offset;
	my $end1 = $f[$h{end1}] + $offset;

	# contig2
	my $contig2 = $f[$h{contig2}];
	my $start2 = $f[$h{start2}] - $offset;
	my $end2 = $f[$h{end2}] + $offset;

	$sim{$contig1} = {} if (!defined($sim{$contig1}));
	$sim{$contig2} = {} if (!defined($sim{$contig2}));

	if (!defined($sim{$contig1}->{$contig2})) {
	    $sim{$contig1}->{$contig2} = {};
	    $sim{$contig1}->{$contig2}->{rect_count} = 0;
	    $sim{$contig1}->{$contig2}->{rects} = {};
	}
	if (!defined($sim{$contig2}->{$contig1})) {
	    $sim{$contig2}->{$contig1} = {};
	    $sim{$contig2}->{$contig1}->{rect_count} = 0;
	    $sim{$contig2}->{$contig1}->{rects} = {};
	}

	my $index = $sim{$contig1}->{$contig2}->{rect_count};
	$sim{$contig1}->{$contig2}->{rects}->{$index} = [$start1, $end1, $start2, $end2];
	$sim{$contig1}->{$contig2}->{rect_count} = $index + 1;

	$index = $sim{$contig2}->{$contig1}->{rect_count};
	$sim{$contig2}->{$contig1}->{rects}->{$index} = [$start1, $end1, $start2, $end2];
	$sim{$contig2}->{$contig1}->{rect_count} = $index + 1;

    }
    close(IN);
}

#######################################################################################
# traverse reads
#######################################################################################

# contig length table
my %counts;

print STDERR "contact input dir: $idir\n";
my @ifns = <$idir/*>;
scalar(@ifns) > 0 or die "no files found in contact input dir";
print STDERR "traversing ".scalar(@ifns)." input files ";

my $read_count = 0;
foreach my $ifn (@ifns) {
    print STDERR ".";
    open(IN, $ifn) || die $ifn;
    my $header = <IN>;
    my %h = parse_header($header);

    my $ofn = $odir_contact."/".basename($ifn);
    open(OUT, ">", $ofn) || die $ofn;
    print OUT $header;

    # output masked contacts
    $ofn = $odir_masked."/".basename($ifn);
    open(OUTM, ">", $ofn) || die $ofn;
    print OUTM $header;

    while (my $line = <IN>) {
	$read_count++;
#	print STDERR "line: $read_count\n" if ($read_count % 100000 == 0);

	chomp $line;
	my @f = split("\t", $line);

	my $contig1 = $f[$h{contig1}];
	my $coord1 = $f[$h{coord1}];
	my $strand1 = $f[$h{strand1}];
	my $score1 = $f[$h{score1}];

	my $contig2 = $f[$h{contig2}];
	my $coord2 = $f[$h{coord2}];
	my $strand2 = $f[$h{strand2}];
	my $score2 = $f[$h{score2}];

	my $id = $f[$h{id}];

	next if (!defined($contigs{$contig1}) || !defined($contigs{$contig2}));

	my $length1 = $contigs{$contig1}->{length};
	my $length2 = $contigs{$contig2}->{length};
	my $abundance1 = $contigs{$contig1}->{abundance};
	my $abundance2 = $contigs{$contig2}->{abundance};
	my $area = $length1 * $length2;
	my $factor = $length1 * $length2 * $abundance1 * $abundance2;

	# check only one side since similirity is symmetric
	my $mask = 0;
	if (defined($sim{$contig1}) && defined($sim{$contig1}->{$contig2})) {
	    my $count = $sim{$contig1}->{$contig2}->{rect_count};
	    for (my $i = 0; $i < $count; $i++) {
		my ($s1, $e1, $s2, $e2) = @{$sim{$contig1}->{$contig2}->{rects}->{$i}};
		if ($s1 <= $coord1 && $coord1 <= $e1 && $s2 <= $coord2 && $coord2 <= $e2) {
		    $mask = 1;
		    last;
		}
	    }
	}

	if (!defined($counts{$contig1}->{$contig2})) {
	    $counts{$contig1}->{$contig2} = {};
	    $counts{$contig1}->{$contig2}->{contacts} = 0;
	    $counts{$contig1}->{$contig2}->{masked_contacts} = 0;
	    $counts{$contig1}->{$contig2}->{area} = $area;
	    $counts{$contig1}->{$contig2}->{factor} = $factor;
	}
	if (!defined($counts{$contig2}->{$contig1})) {
	    $counts{$contig2}->{$contig1} = {};
	    $counts{$contig2}->{$contig1}->{contacts} = 0;
	    $counts{$contig2}->{$contig1}->{masked_contacts} = 0;
	    $counts{$contig2}->{$contig1}->{area} = $area;
	    $counts{$contig2}->{$contig1}->{factor} = $factor;
	}

	if (!$mask) {
	    if ($contig1 ne $contig2) {
		$counts{$contig1}->{$contig2}->{contacts} += 1;
		$counts{$contig2}->{$contig1}->{contacts} += 1;
	    } else {
		$counts{$contig1}->{$contig2}->{contacts} += 1;
	    }
	    print OUT $line, "\n";

	} else {
	    if ($contig1 ne $contig2) {
		$counts{$contig1}->{$contig2}->{masked_contacts} += 1;
		$counts{$contig2}->{$contig1}->{masked_contacts} += 1;
	    } else {
		$counts{$contig1}->{$contig2}->{masked_contacts} += 1;
	    }
	    print OUTM $line, "\n";
	}
    }
    close(IN);
    close(OUT);
    close(OUTM);
}
print ".";

# output result matrix
print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "contig1\tcontig2\tcontacts\tmasked_contacts\tarea\tfactor\n";
foreach my $contig1 (keys %counts) {
    foreach my $contig2 (keys %{$counts{$contig1}}) {
	my $contacts = $counts{$contig1}->{$contig2}->{contacts};
	my $masked_contacts = $counts{$contig1}->{$contig2}->{masked_contacts};
	my $area = $counts{$contig1}->{$contig2}->{area};
	my $factor = $counts{$contig1}->{$contig2}->{factor};
	print OUT "$contig1\t$contig2\t$contacts\t$masked_contacts\t$area\t$factor\n";
	} }
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
