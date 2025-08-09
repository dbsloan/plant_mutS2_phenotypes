#!/usr/bin/perl 

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 library_file sim_file max_overlap > output\n\n";


my $library_file = shift or die ($usage);
my $sim_file = shift or die ($usage);
my $max_overlap = shift or die ($usage);

my @sim_lines = file_to_array($sim_file);
shift @sim_lines;
my $sim_sum = 0;
my %sim_hash;

for (my $i=0; $i <= $max_overlap; ++$i){
	$sim_hash{$i} = 0;
}

foreach (@sim_lines){
	chomp $_;
	my @sl = split (/\t/, $_);
	$sim_hash{$sl[0]} += $sl[1];
	$sim_sum += $sl[1];
}


my @file_lines = file_to_array($library_file);

print "Treatment\tGenotype\tOverlap\tFrequency\tSimulation\n";

foreach (@file_lines){
	chomp $_;
	my $total_count;
	my %overlap_hash;
	for (my $i=0; $i <= $max_overlap; ++$i){
		$overlap_hash{$i} = 0;
	}
	
	my @sl = split (/\t/, $_);
	my @sl2 = split (/\,/, $sl[2]);
	
	foreach my $file (@sl2){
		my @hist_lines = file_to_array($file);
		shift @hist_lines;
		foreach my $line (@hist_lines){
			chomp $line;
			my @shl = split (/\t/ , $line);
			$overlap_hash{$shl[0]} += $shl[1];
			$total_count += $shl[1];
		}
	}
	for (my $i=0; $i <= $max_overlap; ++$i){
		print "$sl[0]\t$sl[1]\t$i\t", $overlap_hash{$i} / $total_count, "\t", $sim_hash{$i} / $sim_sum, "\n";
	}
}