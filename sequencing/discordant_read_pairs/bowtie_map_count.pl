#!/usr/bin/perl 

use strict;
use warnings;
use sloan;


my $usage = "\nUSAGE: perl $0 bowtie_err_file > output\n\n";

my $file = shift or die ($usage);
my @lines = file_to_array($file);

my $map_count= 0;

if ($lines[3] =~ /^\s+(\d+)\s+/){
	$map_count += $1;
}

if ($lines[4] =~ /^\s+(\d+)\s+/){
	$map_count += $1;
}

if ($lines[7] =~ /^\s+(\d+)\s+/){
	$map_count += $1;
}

print "$file\t$map_count\n";