#!/usr/bin/perl 

use strict;
use warnings;
use sloan;
use Bio::SearchIO;


my $usage = "\nUSAGE: perl $0 fasta_file min_pos max_pos frag1_len frag2_len PE_window reps > output\n\n";


my $fasta_file = shift or die ($usage);
my $min_pos = shift or die ($usage);
my $max_pos = shift or die ($usage);
my $frag1_len = shift or die ($usage);
my $frag2_len = shift or die ($usage); #should be smaller than frag1 to ensure frag1 is listed first in blast HSPs.
my $PE_window = shift or die ($usage);
my $reps = shift or die ($usage);

my %overlap_hash;

my %fasta = fasta2hash($fasta_file);

my @headers = sort keys %fasta;


my $FH_BLASTQ = open_output ("TEMP_BLAST_SIM_QUERY");

for (my $i = 1; $i <= $reps; ++$i){
	
	my $start1 = int(rand($max_pos - $min_pos + 1)) + $min_pos;
	my $start2 = int(rand($PE_window + 1)) + $start1;
	my $strand1;
	if (rand() > 0.5){
		$strand1 = "fwd";
	}else{
		$strand1 = "rev";
	}

	my $strand2;
	if (rand() > 0.5){
		$strand2 = "fwd";
	}else{
		$strand2 = "rev";
	}
	
	my $frag1;
	my $frag2;
	
	if ($strand1 eq "fwd"){
		$frag1 = substr ($fasta{$headers[0]}, $start1 - 1, $frag1_len);
	}elsif($strand1 eq "rev"){
		$frag1 = revcom (substr ($fasta{$headers[0]}, $start1 - $frag1_len - 1, $frag1_len));
	}else{
		die("\nERROR: stand1 is not fwd or rev\n\n");
	}

	if ($strand2 eq "fwd"){
		$frag2 = substr ($fasta{$headers[0]}, $start2 - 1, $frag2_len);
	}elsif($strand2 eq "rev"){
		$frag2 = revcom (substr ($fasta{$headers[0]}, $start2 - $frag2_len - 1, $frag2_len));
	}else{
		die("\nERROR: stand2 is not fwd or rev\n\n");
	}

	my $seq = $frag1 . $frag2;
	
	print $FH_BLASTQ ">Seq_$i\n$seq\n";

}

system("blastn -task blastn -evalue 0.000001 -db $fasta_file -query TEMP_BLAST_SIM_QUERY -out TEMP_BLAST_SIM_OUT");

my $SearchIO_obj = new Bio::SearchIO(-format => 'blast', -file => "TEMP_BLAST_SIM_OUT");
while (my $result_obj = $SearchIO_obj->next_result){

	my $query_name = $result_obj->query_name;
	if (my $hit_obj = $result_obj->next_hit){
	
		my $hsp_obj1 = $hit_obj->next_hsp;
		my $hsp_obj2 = $hit_obj->next_hsp;
	
		unless ($hsp_obj1 and $hsp_obj2){
			print STDERR "\nWARNING: $query_name does not have 2 HSPs\n";
			next;
		}
		
		$hsp_obj1->start('query') != 1 and print STDERR "\nWARNING: $query_name HSP1 does not start at 1\n" and next;
		$hsp_obj2->end('query') != $frag1_len + $frag2_len and print STDERR "\nWARNING: $query_name HSP2 does not reach end\n" and next;
		$hsp_obj1->end('query') < $hsp_obj2->start('query') - 1 and print STDERR "\nWARNING: $query_name HSP1 and HSP2 do not have full coverage\n" and next;

		++$overlap_hash{$hsp_obj1->end('query') - $hsp_obj2->start('query') + 1};
	
	}else{
		print STDERR "\nWARNING: not hit for $query_name\n";
	}

}

print "OverlapSize\tCount\n";
foreach (sort {$a <=> $b} keys %overlap_hash){
	print "$_\t$overlap_hash{$_}\n";
}


#unlink("TEMP_BLAST_SIM_QUERY");
#unlink("TEMP_BLAST_SIM_OUT");