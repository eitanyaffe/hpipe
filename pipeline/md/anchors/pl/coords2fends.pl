#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Data::Dumper;
use Switch;

if ($#ARGV == -1) {
    print STDERR "usage: $0 <fends> <read length> <trim length> <discard facing threshold> <max segment length> <contig table> <output mat prefix> <correct coords of trimmed reads> <input file1> [input file2] ...\n";
    exit 1;
}

my $in_fends_fn = $ARGV[0];
my $read_len = $ARGV[1];
my $trim_len = $ARGV[2];
my $discard_facing_threshold = $ARGV[3];
my $max_len = $ARGV[4];
my $contig_table = $ARGV[5];
my $mat_fn_prefix = $ARGV[6];
my $correct_trimmed = $ARGV[7];
shift;shift;shift;shift;shift;shift;shift;shift;
my @ifns =  @ARGV;
print STDERR "Input files: ", join("\n", @ifns), "\n";

# table with fends
our %fends;

# hash table to quickly get from approximate coord to fend
our %coord2index;

# contig size table
our %contigs;

##########################################################################################
# read contig table
##########################################################################################

open(IN, $contig_table) || die $contig_table;
my $header = <IN>;
my %h = parse_header($header);
print STDERR "Reading contig table $contig_table into hash...\n";
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $contig = $f[$h{contig}];
    my $len = $f[$h{length}];
    $contigs{$contig} = $len;
}
close(IN);

##########################################################################################
# read fends file
##########################################################################################

print STDERR "Reading input file $in_fends_fn into hash...\n";
open(IN, $in_fends_fn) || die $in_fends_fn;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $fend = $f[$h{fend}];
    my $contig = $f[$h{contig}];
    my $coord = $f[$h{coord}];
    my $strand = $f[$h{strand}];
    my $frag = $f[$h{frag}];
    my $frag_len = $f[$h{frag_len}];

    !defined($fends{$fend}) or die "non-unique fend";
    $fends{$fend} = {};
    $fends{$fend}->{fend} = $fend;
    $fends{$fend}->{frag} = $frag;
    $fends{$fend}->{strand} = $strand;
    $fends{$fend}->{contig} = $contig;
    $fends{$fend}->{coord} = $coord;
    $fends{$fend}->{frag_len} = $frag_len;

    if (!defined($coord2index{$contig}))
    {
	$coord2index{$contig} = {};
	$coord2index{$contig}->{coords} = {};
    }
    $coord2index{$contig}->{coords}->{$coord} = $fend;
}
close(IN);

# compute sorted coords per contig
foreach my $contig (keys %coord2index)
{
    my @sorted = sort {$a <=> $b} keys %{$coord2index{$contig}->{coords}};
    $coord2index{$contig}->{sorted_coords} = \@sorted;
}

##########################################################################################
# parse pair files
##########################################################################################

my $read_count = 0;
my $one_side_cutter = 0;

my @types = ("s0", "s1", "s2");
my %fend_matrix;
my %read_stats;

# read types
my @fields = ("genomic_with_sites", "genomic_no_site", "genomic_side", "side", "self_ligation", "seglength_out_of_range", "cis_0_1k", "cis_1k_10k", "cis_10k_100k", "cis_100k_1m", "cis_1m_max", "trans");

foreach my $type (@types) {
    $fend_matrix{$type} = {};
}


# init stat containers first time hitting the contig

foreach my $contig (keys %contigs) {
    $read_stats{$contig} = {};
    $read_stats{$contig}->{total} = 0;
    foreach my $ttype (@types) {
	$read_stats{$contig}->{$ttype} = {};
	$read_stats{$contig}->{$ttype}->{count} = 0;
	foreach my $field (@fields) {
	    $read_stats{$contig}->{$ttype}->{$field} = 0;
	}
    }
}

foreach my $ifn (@ifns) {

    print STDERR "traversing file: $ifn\n";
    open(IN, $ifn) || die $ifn;

    $header = <IN>;
    %h = parse_header($header);

    while (my $line = <IN>) {
	$read_count++;
	print STDERR "line: $read_count\n" if ($read_count % 100000 == 0);
	chomp $line;
	my @f = split("\t", $line);

	my $contig1 = $f[$h{contig1}];
	my $coord1 = $f[$h{coord1}];
	my $strand1 = $f[$h{strand1}] == 1 ? "+" : "-";
	my $contig2 = $f[$h{contig2}];
	my $coord2 = $f[$h{coord2}];
	my $strand2 = $f[$h{strand2}] == 1 ? "+" : "-";

	if ($correct_trimmed eq "T") {
	    $coord1 -= $trim_len if ($strand1 eq "+");
	    $coord2 -= $trim_len if ($strand2 eq "+");
	    $coord1 += $read_len+$trim_len+1 if ($strand1 eq "-");
	    $coord2 += $read_len+$trim_len+1 if ($strand2 eq "-");
	}

	next if (!defined($contigs{$contig1}) || !defined($contigs{$contig2}));

	my ($fend1, $fwd_dist1, $bck_dist1) = coord_to_fend($contig1, $coord1, $strand1);
	my ($fend2, $fwd_dist2, $bck_dist2) = coord_to_fend($contig2, $coord2, $strand2);

	# identify type of restriction
	my $type;
	my $on_cutter_site1 = ($bck_dist1 == 2);
	my $on_cutter_site2 = ($bck_dist2 == 2);
	if ($on_cutter_site1 && $on_cutter_site2) {
	    $type = "s2";
	} elsif (!$on_cutter_site1 && !$on_cutter_site2) {
	    $type = "s0";
	} else {
	    $type = "s1";
	}

	# to estimate chance of one side to fall on cutter
	$one_side_cutter++ if ($on_cutter_site1);

	$read_stats{$contig1}->{total}++;
	$read_stats{$contig2}->{total}++ if ($contig1 ne $contig2);

	# total read count for type
	$read_stats{$contig1}->{$type}->{count}++;
	$read_stats{$contig2}->{$type}->{count}++ if ($contig1 ne $contig2);

	my $face_towards = (($coord1 < $coord2) && ($strand1 eq "+") && ($strand2 eq "-")) ||
	    (($coord2 < $coord1) && ($strand1 eq "-") && ($strand2 eq "+"));

	# first discard the close cis facing products - they are probably no restriction events
	if (($contig1 eq $contig2) && $face_towards && (abs($coord1-$coord2) < $discard_facing_threshold)) {
	    if (count_fends($contig1, $coord1, $coord2) == 0)
	    {
		$read_stats{$contig1}->{$type}->{genomic_no_site}++;
	    } else {
		$read_stats{$contig1}->{$type}->{genomic_with_sites}++;
	    }
	    next;
	}

	my $side_dist1 = ($strand1 eq "+") ? $contigs{$contig1} - $coord1 : $coord1;
	my $side_dist2 = ($strand2 eq "+") ? $contigs{$contig2} - $coord2 : $coord2;

	# genomic crossing contigs
	if (($contig1 ne $contig2) && ($side_dist1 + $side_dist2) < $discard_facing_threshold) {
	    $read_stats{$contig1}->{$type}->{genomic_side}++;
	    $read_stats{$contig2}->{$type}->{genomic_side}++;
	    next;
	}

	# side of contigs
	if (($fend1 == -1 || $fend2 == -1)) {
	    $read_stats{$contig1}->{$type}->{side}++;
	    $read_stats{$contig2}->{$type}->{side}++ if ($contig1 ne $contig2);
	    next;
	}

	# compute stats for different events
	my $frag1 = $fends{$fend1}->{frag};
	my $frag2 = $fends{$fend2}->{frag};

	# skip if ligation is within a single fragment
	if ($frag1 == $frag2) {
	    if ($face_towards) {
		$read_stats{$contig1}->{$type}->{genomic_with_sites}++;
	    } else {
		$read_stats{$contig1}->{$type}->{self_ligation}++;
	    }
	    next;
	}

	# finally if inferred segment distance is too long we discard read as well
        # this is mainly relevant for 6-cutters
	my $segment_len = $fwd_dist1 + $fwd_dist2;
	if ($max_len > 0 && $segment_len > $max_len)
	{
	    $read_stats{$contig1}->{$type}->{seglen_out_of_range}++;
	    $read_stats{$contig1}->{$type}->{seglen_out_of_range}++ if ($contig1 ne $contig2);
	    next;
	}

	# now, breakdown the normal ligation events according to distance
	my $dist = abs($fends{$fend1}->{coord} - $fends{$fend2}->{coord});
	if ($contig1 eq $contig2) {
	    if ($dist < 1000) {
		$read_stats{$contig1}->{$type}->{cis_0_1k}++;
	    } elsif ($dist < 10000) {
		$read_stats{$contig1}->{$type}->{cis_1k_10k}++;
	    } elsif ($dist < 100000) {
		$read_stats{$contig1}->{$type}->{cis_10k_100k}++;
	    } elsif ($dist < 1000000) {
		$read_stats{$contig1}->{$type}->{cis_100k_1m}++;
	    } else {
		$read_stats{$contig1}->{$type}->{cis_1m_max}++;
	    }
	} else {
	    $read_stats{$contig1}->{$type}->{trans}++;
	    $read_stats{$contig2}->{$type}->{trans}++;
	}

	# track number of covered fends
	my $fend_small = ($fend1 <= $fend2) ? $fend1 : $fend2;
	my $fend_large = ($fend1 <= $fend2) ? $fend2 : $fend1;
	$fend_matrix{$type}->{$fend_small} = {} if !defined($fend_matrix{$type}->{$fend_small});
	$fend_matrix{$type}->{$fend_small}->{$fend_large} = 0 if !defined($fend_matrix{$type}->{$fend_small}->{$fend_large});
	$fend_matrix{$type}->{$fend_small}->{$fend_large}++;
    }
    close(IN);
}

######################################################################################################
# write mat file
######################################################################################################

foreach my $type (@types)
{
    my $mat_fn = $mat_fn_prefix."_$type.mat";
    print STDERR "writing file: $mat_fn\n";
    open(OUT, ">", $mat_fn) || die;
    print OUT "fend1\tfend2\tcount\n";
    foreach my $fend1 (sort { $a <=> $b } keys %{$fend_matrix{$type}})
    {
	foreach my $fend2 (sort { $a <=> $b } keys %{$fend_matrix{$type}->{$fend1}})
	{
	    my $count = $fend_matrix{$type}->{$fend1}->{$fend2};
	    print OUT $fend1, "\t" ,$fend2, "\t", $count, "\n";
	}
    }
    close(OUT);
}

######################################################################################################
# Write stats
######################################################################################################

# global stats for all types
my $gstats_fn = $mat_fn_prefix."_all.mat.stats";
open(OUT, ">", $gstats_fn) || die;
print OUT "total_reads\tside1_on_cutter_site\tdiscarded_reads\n";
print OUT $read_count, "\t", $one_side_cutter, "\n";
close(OUT);

my $p = $one_side_cutter / $read_count;

foreach my $type (@types)
{
    my $read_stats_fn = $mat_fn_prefix."_$type.read.stats";
    open(OUT, ">", $read_stats_fn) || die;
    print OUT "contig\ttotal\texpected_for_type\t", join("\t", @fields), "\n";

    foreach my $contig (keys %contigs) {
	# compute expected sizes of datasets
	my $count = $read_stats{$contig}->{total};
	my %expected;
	$expected{"s2"} = $p * $p * $count;
	$expected{"s0"} = (1-$p) * (1-$p) * $count;
	$expected{"s1"} = 2 * (1-$p) * $p * $count;

	print OUT $contig, "\t";
	print OUT $read_stats{$contig}->{$type}->{count}, "\t", int($expected{$type});
	for my $field (@fields) {
	    print OUT "\t", $read_stats{$contig}->{$type}->{$field}};
	print OUT "\n";
    }
    close(OUT);
}

print STDERR "coord2fends done\n";

######################################################################################################
# Subroutines
######################################################################################################

# check if C is between and A,B
sub between
{
    my ($A, $B, $C) = @_;
    $B = $B - $A;
    $C = $C - $A;
    return (($B>$C && $C>0) || ($B<$C && $C<0));
}

sub add_fend
{
    my ($contig, $bin, $fend, $strand) = @_;
    $coord2index{$contig} = {} if !defined($coord2index{$contig});
    $coord2index{$contig}->{$bin} = {} if !defined($coord2index{$contig}->{$bin});

    # mark bin if multiple fends map to it
    if (!defined($coord2index{$contig}->{$bin}->{$strand}))
    {
	$coord2index{$contig}->{$bin}->{$strand} = $fend;
    }
    else
    {
	$coord2index{$contig}->{$bin}->{$strand} = -1;
    }
}

sub apprx_lines
{
    my ($fn) = @_;
    my $tmp = "/tmp/".getlogin()."_apprx_lines.tmp";
    system("head -n 100000 $fn > $tmp");
    my $size_head = -s $tmp;
    my $size_all = -s $fn;
    $size_head > 0 or die;
    return (int($size_all/$size_head*100000));
}

sub perc_str
{
    my ($n, $total) = @_;
    return ($n." (".(int(1000 * $n / $total)/10)."%)");
}

sub perc_str2
{
    my ($n, $total) = @_;
    return ((int(1000 * $n / $total)/10)."%");
}

# returns first element above/below value in sorted array
sub binary_search
{
    my $arr = shift;
    my $value = shift;
    my $above = shift;

    my $left = 0;
    my $right = $#$arr;

    while ($left <= $right) {
	my $mid = ($right + $left) >> 1;
	my $c = $arr->[$mid] <=> $value;
	return $mid if ($c == 0);
	if ($c > 0) {
	    $right = $mid - 1;
	} else {
	    $left  = $mid + 1;
	}
    }
    $left = -1 if ($left > $#$arr);
    $right = -1 if ($right < 0);
    return (($above eq "+") ? $left : $right);
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

# returns (fend, fwd distance, bck distance)
sub coord_to_fend
{
    my ($contig, $coord, $strand) = @_;
    defined($coord2index{$contig}) or return (-1, -1, -1);

    # opposite strand
    my $ostrand = ($strand eq "+") ? "-" : "+";

    # search fwd and bck in coord list
    my $index_p = binary_search($coord2index{$contig}->{sorted_coords}, $coord, $strand);
    my $index_m = binary_search($coord2index{$contig}->{sorted_coords}, $coord, $ostrand);

    # get fwd and bck coords
    my $fcoord_p = ($index_p != -1) ? $coord2index{$contig}->{sorted_coords}[$index_p] : -1;
    my $fcoord_m = ($index_m != -1) ? $coord2index{$contig}->{sorted_coords}[$index_m] : -1;

    # get fwd and bck distance
    my $dist_p = ($index_p != -1) ? abs($fcoord_p - $coord) : -1;
    my $dist_m = ($index_m != -1) ? abs($fcoord_m - $coord) : -1;

    # result fend
    my $fend = -1;

    # get fend only on fwd side
    if ($index_p != -1) {
	# sanity check
	($fcoord_p != -1) and defined ($coord2index{$contig}->{coords}->{$fcoord_p}) and defined($fcoord_p) or die;
	$fend = $coord2index{$contig}->{coords}->{$fcoord_p};
    }

    return ($fend, $dist_p, $dist_m);
}

sub count_fends
{
    my ($contig, $coord1, $coord2) = @_;
    defined($coord2index{$contig}) or return 0;

    return 0 if ($coord1 == $coord2);

    my $coord_l = ($coord1 < $coord2) ? $coord1 : $coord2;
    my $coord_r = ($coord1 < $coord2) ? $coord2 : $coord1;

    # search fwd and bck in coord list
    my $index_p = binary_search($coord2index{$contig}->{sorted_coords}, $coord_l, "+");
    my $index_m = binary_search($coord2index{$contig}->{sorted_coords}, $coord_r, "-");

    return 0 if (($index_p == -1) || ($index_m == -1));
    return ($index_m>$index_p) ? ($index_m-$index_p+1) : 0;
}
