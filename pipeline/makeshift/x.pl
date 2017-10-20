open OUT,">",$ARGV[0]; while (my $line = <STDIN>) { print $line; print OUT $line;} close(OUT);
