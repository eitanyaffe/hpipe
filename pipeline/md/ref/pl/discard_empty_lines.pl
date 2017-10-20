
while (my $line = <STDIN>) {
    chomp($line);
    next if ($line eq "");
    print $line, "\n";
}
