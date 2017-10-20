#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Scalar::Util qw(looks_like_number);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ofn = $ARGV[1];

print STDERR "reading file: $ifn\n";
open(IN, $ifn) || die $ifn;

print STDERR "writing file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "id\tparent_ids\ttype\tdesc\troot\n";

my %term;
my $meta = 1;
my $in_term = 0;
while (my $line = <IN>) {
    chomp($line);
    $meta = 0 if ($line eq "" && $meta);
    next if ($meta);

    # print $line, "\n";
    if ($line eq "[Term]") {
	%term = ();
	$term{parent_ids} = {};
	$in_term = 1;
	next;
    }
    if ($line eq "" && defined($term{id})) {
	next if (!$in_term);
	$in_term = 0;
	defined($term{name}) and defined($term{namespace}) or die;
	if (!defined($term{is_obsolete})) {
	    my $root = $term{namespace} eq $term{name} ? "T" : "F";
	    print OUT $term{id}, "\t", join(";", keys %{$term{parent_ids}}), "\t", $term{namespace}, "\t", $term{name}, "\t", $root, "\n";
	}
	next;
    }

    my $ix = index($line, ":");
    next if ($ix == -1);
    my $field = substr($line, 0, $ix);
    my $value = substr($line, $ix+2);
    # print "> $field: $value\n";
    $term{id} = $value if ($field eq "id");
    $term{name} = $value if ($field eq "name");
    $term{namespace} = $value if ($field eq "namespace");
    $term{is_obsolete} = 1 if ($field eq "is_obsolete" && $value eq "true");

    if ($field eq "is_a") {
	my $ip = index($value, "!");
	$ip != -1 or die;
	my $pid = substr($value, 0, $ip-1);
	# print "> add parent: \'$pid\'\n";
	$term{parent_ids}->{$pid} = 1;
    }
}
close(IN);
close(OUT);
