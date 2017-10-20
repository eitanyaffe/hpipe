my $prefix = $ARGV[0];
my $ifn = $ARGV[1];

my %genes;
open(IN, $ifn) || die;
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
    chomp($line);
    my @f = split("\t", $line);
    my $gene = $f[$h{gene}];
    next if ($f[$h{truncated}] eq "T");
    $genes{$gene} = 1;
}

my $skip = 0;
while (my $line = <STDIN>) {
    chomp($line);
    next if ($line eq "");
    if (substr($line, 0, 1) eq ">") {
	my $gene = $prefix.substr($line, index($line, "_"), index($line, "|")-index($line, "_"));
	$skip = !defined($genes{$gene});
	print ">$gene\n" if (!$skip);
    } else {
	print $line, "\n" if (!$skip);
    }

}

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
