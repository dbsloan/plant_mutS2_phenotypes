#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use sloan;

my $usage = "\nUSAGE: perl $0 bin_size max_bin breaks_filename(s) > output\n\n";

my $bin_size = shift or die ($usage);
my $max_bin = shift or die ($usage);

my @file_list;

while (my $file = shift){
	push (@file_list, $file);
}

@file_list or die ($usage);

my %bin_hash;

foreach (@file_list){
	my @lines = file_to_array($_);
	shift @lines;
	foreach my $break_line (@lines){
		my @sl = split (/\t/, $break_line);
		#print STDERR "bin_size\t$sl[1]\t$sl[2]\n";		
		my $pos1 = ceil ($sl[1] / $bin_size) * $bin_size;
		my $pos2 = ceil ($sl[2] / $bin_size) * $bin_size;
		++$bin_hash{$pos1};
		++$bin_hash{$pos2};
	}	
}

print "Position\tRead_Count\n";

for (my $i = $bin_size; $i <= $max_bin; $i += $bin_size){
	if (exists ($bin_hash{$i})){
		print "$i\t$bin_hash{$i}\n";
	}else{
		print "$i\t0\n";
	}
}
