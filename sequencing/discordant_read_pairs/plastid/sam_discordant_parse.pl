#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 sam_file read_len ref_length ir_start ir_end boundary\n\n";

my $file = shift or die ($usage); #sam file from bowtie2 mapping
my $read_len = shift or die ($usage); #length of sequence reads. Must all be the same after previous tuncation script
my $ref_len = shift or die ($usage); #length of reference cp genome used in mapping (should have second copy of IR removed
my $ir_start = shift or die ($usage); #start position of IR in reference
my $ir_end = shift or die ($usage); #end position of IR in reference
my $boundary = shift or die ($usage); #length of flanking region used to scan for read pairs that are falsely identified as discordant because they represent real IR-SSC or IR-LSC connections not available in the reference used for mapping

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
	
	#filter IR-LSC connection. Including 5 bp of wiggle room at IR start because there appears to be some mapping tolerance there.
	if ($first_strand eq '-' and $second_strand eq '-'){
		if ($sl[3] <= $boundary - $read_len + 1 and $sl2[3] >= $ir_start -5 and $sl2[3] <= $ir_start + $boundary - $read_len){
			next;
		}elsif($sl[3] >= $ir_start - 5 and $sl[3] <= $ir_start + $boundary - $read_len and $sl2[3] <= $boundary - $read_len + 1){
			next;
		}
	}

	#filter IR-SSC connections
        if ($first_strand eq '+' and $second_strand eq '+'){
		if($sl[3] >= $ref_len - $boundary + 1 and $sl2[3] >= $ir_end - $boundary +1 and $sl2[3] <= $ir_end){
			next;
		}elsif($sl[3] >= $ir_end - $boundary +1 and $sl[3] <= $ir_end and $sl2[3] >= $ref_len - $boundary + 1){
			next;
		}
	}

	print "$sl[0]\t$sl[3]\t$sl2[3]\t$first_strand\t$second_strand\n";	
	
}
