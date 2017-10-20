#!/usr/bin/perl -w

use strict;
use XML::Parser;
use LWP::Simple;  # used to fetch the chatterbox ticker
use strict;
# use warnings FATAL => qw(all);

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
print OUT "uniparc\tuniprots\n";
$parser->parse($fh, ProtocolEncoding => 'US-ASCII');
close(OUT);

# The Handlers
sub hdl_start{
    my ($p, $elt, %atts) = @_;
#    print "p=$p, elt=$elt, ", join(";", keys %atts), "\n";
    if ($elt eq 'entry') {
	%entry = ();
	$entry{id} = $atts{id};
	$entry{in_accession} = 0;
	$entry{uniprots} = {};
    }
    $entry{in_accession} = 1 if ($elt eq 'accession');
    if ($elt eq 'dbReference') {
	return 1 if (!defined($atts{type}) || ($atts{type} ne "UniProtKB/Swiss-Prot" && $atts{type} ne "UniProtKB/TrEMBL"));
	return 1 if (!defined($atts{active}) || $atts{active} ne "Y");
	return 1 if (!defined($atts{id}));
	$entry{uniprots}->{$atts{id}} = 1;
    }
}

sub hdl_char {
    my ($p, $string) = @_;
    # print "char string=$string\n";
    $entry{id} = $string if ($entry{in_accession});
}

sub hdl_end {
    my ($p, $elt) = @_;
    $entry{in_accession} = 0 if ($elt eq 'accession');
    return unless $elt eq 'entry';
    my @uniprots = keys %{$entry{uniprots}};
    print OUT $entry{id}, "\t", scalar(@uniprots) > 0 ? join(";", @uniprots): "NA", "\n";
}

sub hdl_def { }  # We just throw everything else
