#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE perl $0 samtools_depth_file window_size\n\n";

my $depth_file = shift or die ($usage);
my $window_size = shift or die ($usage);

my @lines = file_to_array($depth_file);

my $file_num = scalar (split (/\t/, $lines[0])) - 2;

my $base_count=0;
my $pos = 0;
my @window_sums;

print "Position";

for (my $i=1; $i <= $file_num; ++$i){
	print "\tFile_$i";
}
print "\n";

foreach (@lines){
	chomp $_;
	my @sl = split (/\t/, $_);
	$pos + 1 == $sl[1] or die ("Base positions are not continuous and starting from 1. Check that samtools depth was run with -a option\n");

	if ($base_count == $window_size){
		
		print $pos - $base_count/2;
		for (my $i = 1; $i <= $file_num; ++$i){
			print "\t", $window_sums[$i-1] / $base_count;
		}
		print "\n";
		$base_count = 0;
		@window_sums = map { 0 } @window_sums;
	}
	
	++$pos;
	++$base_count;
	for (my $i = 1; $i <= $file_num; ++$i){
		$window_sums[$i-1] += $sl[$i+1];
	}	
}
