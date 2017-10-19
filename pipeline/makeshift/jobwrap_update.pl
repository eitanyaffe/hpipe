use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);
use Getopt::Long;

# title of job
my $title   = "";

# send email
my $email = "";

# send updates when regex met
my $regex = "";

GetOptions ("title=s" => \$title,
	    "regex=s" => \$regex,
	    "email=s" => \$email)
    or die("Error in command line arguments\n");

$regex ne "" or die;
$email ne "" or die;
$title ne "" or die;

my @lines;
while (my $line = <STDIN>)
{
    print $line;
    push(@lines, $line);
    if ($line =~ $regex) {
	my $email_command = sprintf("echo command line: %s | mail -s \"[JW DONE] Job %s ended %s\" -a %s %s", $command, $title, $desc, $ofn, $email);
	system($email_command) == 0 or die;
    }
}
