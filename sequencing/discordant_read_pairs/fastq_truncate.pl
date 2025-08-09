#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 fastq_input length > fastq_output\n\n";

my $file = shift or die ($usage);
my $length = shift or die ($usage);

my $FH = open_file($file);

while (my $line = <$FH>){
	print "$line";
	$line = <$FH>;
	chomp $line;
	print substr ($line, 0, $length), "\n";
	$line = <$FH>;
	print $line;
	$line = <$FH>;
	chomp $line;
	print substr ($line, 0, $length), "\n";
}

