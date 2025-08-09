#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 sam_file read_len\n\n";

my $file = shift or die ($usage); #sam file from bowtie2 mapping
my $read_len = shift or die ($usage); #length of sequence reads. Must all be the same after previous tuncation script

my $FH = open_file($file);

my %flags;

print "Name\tRead1Pos\tRead2Pos\tRead1Strand\tRead2Strand\n";

while (my $line = <$FH>){

	#skip header lines
	$line =~ /^\@/ and next;
	my @sl = split (/\t/, $line);

	my $flag = $sl[1];
	my $first_strand;
	my $second_strand;
	
	#Parse sam bit flag
	$flag >= 128 and $flag -= 128;
	$flag >= 64 and $flag -= 64;
	
	if ($flag >= 32){
		$second_strand = '-';
		$flag -= 32;
	}else{
		$second_strand = '+'
	}

	if ($flag >= 16){
		$first_strand = '-';
		$flag -= 16;
	}else{
		$first_strand = '+'
	}

	#exclude concordnant mapping reads
	$flag >= 2 and next;
	
	
	my $line2 = <$FH>;
	my @sl2 = split (/\t/, $line2);

	#incomplate mappings within indels
	unless ($sl[5] eq $read_len . "M" and $sl2[5] eq $read_len . "M"){
		next;
	}

	#exlcude mappings with substitutions
	$line =~ /XM\:i\:0/ or next;
	$line2 =~ /XM\:i\:0/ or next;
	
	print "$sl[0]\t$sl[3]\t$sl2[3]\t$first_strand\t$second_strand\n";	
	
}
