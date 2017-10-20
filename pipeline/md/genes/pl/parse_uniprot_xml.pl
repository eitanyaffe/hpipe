#!/usr/bin/perl -w

use strict;
use XML::Parser;
use LWP::Simple;  # used to fetch the chatterbox ticker
use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <xml ifn> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ofn = $ARGV[1];

print "reading xml: $ifn\n";

my %entry;      # Hashref containing infos on a entry

my $fh = IO::File->new($ifn) or die;

my $parser = new XML::Parser ( Handlers => {   # Creates our parser object
    Start   => \&hdl_start,
    Char     => \&hdl_char,
    End     => \&hdl_end,
    Default => \&hdl_def
			       });

print "writing table: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT "uniprot\tGO_ids\n";
$parser->parse($fh);
close(OUT);

# The Handlers
sub hdl_start{
    my ($p, $elt, %atts) = @_;
    # print "start elt=$elt, ", join(";", keys %atts), "\n";
    if ($elt eq 'entry') {
	%entry = ();
	$entry{in_accession} = 0;
	$entry{in_GO} = 0;
	$entry{GO_ids} = {};
    }

    $entry{in_accession} = 1 if ($elt eq 'accession');
    if ($elt eq 'dbReference' && defined($atts{type}) && $atts{type} eq "GO") {
	$entry{in_go_reference} = 1;
	$entry{GO_ids}->{$atts{id}} = 1 if (defined($atts{id}));
    }
}

sub hdl_char {
    my ($p, $string) = @_;
    # print "char string=$string\n";
    $entry{prot_id} = $string if ($entry{in_accession});
}

sub hdl_end {
    my ($p, $elt) = @_;
    $entry{in_accession} = 0 if ($elt eq 'accession');
    $entry{in_go_reference} = 0 if ($elt eq 'dbReference');
    return if ($elt ne 'entry');
    print OUT $entry{prot_id}, "\t", join(";", keys %{$entry{GO_ids}}), "\n";
}

sub hdl_def { }  # We just throw everything else
