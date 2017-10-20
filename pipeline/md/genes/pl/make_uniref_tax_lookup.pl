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
    End     => \&hdl_end,
    Default => \&hdl_def,
			       });

print "writing table: $ofn\n";
open(OUT, ">", $ofn) || die;
print OUT "uniref\ttax_ids\n";
$parser->parse($fh);
close(OUT);

# The Handlers
sub hdl_start{
    my ($p, $elt, %atts) = @_;
#    print "p=$p, elt=$elt, ", join(";", keys %atts), "\n";
    if ($elt eq 'entry') {
	%entry = ();
	$entry{ref_id} = $atts{id};
	$entry{tax_ids} = {};
    }
    if ($elt eq 'property' && defined($atts{type}) && $atts{type} eq "NCBI taxonomy") {
	defined($atts{value}) or die;
	$entry{tax_ids}->{$atts{value}} = 1;
    }
}

sub hdl_end{
    my ($p, $elt) = @_;
    return unless $elt eq 'entry' and scalar(keys %{$entry{tax_ids}}) > 0;
    print OUT $entry{ref_id}, "\t", join(";", keys %{$entry{tax_ids}}), "\n";
}

sub hdl_def { }  # We just throw everything else
