#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use List::Util qw(sum);
use Getopt::Long;
use Cwd 'abs_path';

my $path = abs_path($0);

if ($#ARGV == -1) {
    print STDERR "usage: $0 [options] -- <command>\n";
    print STDERR "Options:\n";
    print STDERR "   -title <string>: job title\n";
    print STDERR "   -ldir <path>: directory to save log file\n";
    print STDERR "   -update <regex>: send update whenever log meets regex pattern\n";
    print STDERR "   -email <email address>: send results to email\n";
    print STDERR "   -always_email <0|1>: send results to email (default is true)\n";
    print STDERR 'Example: ./jobwrap.pl -title ab -email eitany@stanford.edu -- ls *', "\n";
    exit 1;
}

# title of job
my $title   = "";

# log directory
my $ldir = "logs";

# send email
my $email = "";
my $always_email = 1;

# send updates when regex met
my @regexs;

GetOptions ("title=s" => \$title,
	    "ldir=s"  => \$ldir,
	    "regex=s" => \@regexs,
	    "always_email=s" => \$always_email,
	    "email=s" => \$email)
    or die("Error in command line arguments\n");

$title ne "" or die "job title must be specified (-title)";
$email ne "" or die "emaikl must be specified (-email)";

my $ofn = $ldir."/".$title;
my $command = join(" ", @ARGV);

print "## Job title: $title\n";
print "## Saving stdin/stdout to file: $ofn\n";
print "## Email: $email\n";
print "## Command: $command\n";
print "###############################################\n";

system("mkdir -p ".$ldir) == 0 or die;
my $lcommand = $command.' 2>&1 | perl -e \'open OUT,">","'.$ofn.'"; while (my $line = <STDIN>) { print $line; print OUT $line;} close(OUT);\'; exit ${PIPESTATUS[0]}';
# print "lcommand: $lcommand\n";
my $error = (system($lcommand) != 0);
my $jobrc = $? >> 8;

print "###############################################\n";
printf("## command return code: %d\n", $jobrc);

if ($email ne "" && ($always_email || $error)) {
    print "## emailing log file: $email\n";
    my $desc = $jobrc == 0 ? "successfully" : sprintf("with an error, return code: %d", $jobrc);
    my $email_command = sprintf("echo command line: %s | mail -s \"[JW DONE] Job %s ended %s\" -a %s %s", $command, $title, $desc, $ofn, $email);
    # print "email command: $email_command\n";
    system($email_command) == 0 or die;
}
