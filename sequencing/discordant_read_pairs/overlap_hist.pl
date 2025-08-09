#!/usr/bin/perl 

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 breaks_file > output\n\n";


my $breaks_file = shift or die ($usage);

my @lines = file_to_array($breaks_file);

shift @lines;

my %overlap_hash;

foreach (@lines){
	chomp $_;
	my @sl = split (/\t/, $_);
	
	if ($sl[5] eq "broken"){
		++$overlap_hash{$sl[6]};
	}elsif ($sl[7] eq "broken"){
		++$overlap_hash{$sl[8]};	
	}
}

print "OverlapSize\tCount\n";
foreach (sort {$a <=> $b} keys %overlap_hash){
	print "$_\t$overlap_hash{$_}\n";
}
