#!/usr/bin/perl 

use strict;
use warnings;
use sloan;
use Bio::SearchIO;

my $usage = "\nUSAGE: perl $0 input_discord_file input_fq1 inputfq2 blast_db min_position max_position dist_cutoff > output\n\n";

my $file = shift or die ($usage);
my $file_out = shift or die ($usage);
my $fq1_file = shift or die ($usage);
my $fq2_file = shift or die ($usage);
my $blast_db = shift or die ($usage);
#use min and max mapping positions to cut of ends of reference genome (should already have one copy of IR removed for plastid genome.)
#this will avoid mismapping doing due to circularization or IR issues.
#values for Arabidopsis (no IR version) to avoid reads mapping to last 1000 bp: 1001 127214. 
my $min_pos = shift or die ($usage); 
my $max_pos = shift or die ($usage);
#reports a subset of mismatches where the two ends are within some distance cutoff. 
#Note for the individual break/overlap reporting, there is currently a hardcoded value of 2000 bp. The below value is for summary counts, but it will also filter inputs to the blast analysis, so it will supersede the 2000 value if set to smaller number.
my $dist_cutoff = shift or die ($usage);


my @lines = file_to_array($file);

my $header = shift @lines;
chomp $header;

my $FH_SUMM = open_output($file_out);

print $FH_SUMM "$header\tRead1Status\tRead1Overalap\tRead2Status\tRead2Overalap\n";

#load all the fastq sequences into a hash of arrays (key is seq ID; reads 1 and 2 are the two array elements )
my %seq_HoA;
my $FH1 = open_file($fq1_file);
my $FH2 = open_file($fq2_file);

while (<$FH1>){
	my @sl = split (/\s/, $_);
	my $name = substr ($sl[0], 1);
	my $seq = <$FH1>;
	chomp $seq;
	my $drop_line = <$FH1>;
	$drop_line = <$FH1>;	
	$seq_HoA{$name}[0] = $seq;
}

while (<$FH2>){
	my @sl = split (/\s/, $_);
	my $name = substr ($sl[0], 1);
	my $seq = <$FH2>;
	chomp $seq;
	my $drop_line = <$FH2>;
	$drop_line = <$FH2>;	
	$seq_HoA{$name}[1] = $seq;
}


#setup counter variables to keep track improperly mapped reads.
my $same_neg_R1_lower = 0;
my $same_neg_R1_higher = 0;
my $same_pos_R1_lower = 0;
my $same_pos_R1_higher = 0;;
my $outward = 0;
my $inward_too_far = 0;

#"cut" versions of counters only report mismatches with mapping points within the specified cutoff dist.
my $same_neg_R1_lower_cut = 0;
my $same_neg_R1_higher_cut = 0;
my $same_pos_R1_lower_cut = 0;
my $same_pos_R1_higher_cut = 0;;
my $outward_cut = 0;
my $inward_too_far_cut = 0;

my %hit_pos_HoH;


foreach (@lines){
	chomp $_;
	my @sl = split (/\t/, $_);

	#exclude reads where either of the pair maps outside the specified min-max range.
	if ($sl[1] > $max_pos or $sl[1] < $min_pos or $sl[2] > $max_pos or $sl[2] < $min_pos){
		next;
	}	

	#only do blast search to look for reads with mapping position within the specified range.
	my $do_blast = 0;
	my $dist = abs ($sl[1] - $sl[2]);

	if ($sl[3] eq '-' and $sl[4] eq '-'){
		if ($sl[1] <= $sl[2]){
			++$same_neg_R1_lower;
			$dist_cutoff >= $dist and ++$same_neg_R1_lower_cut and $do_blast = 1;
		}else{
			++$same_neg_R1_higher;	
			$dist_cutoff >= $dist and ++$same_neg_R1_higher_cut and $do_blast = 1;
		}
	}elsif ($sl[3] eq '+' and $sl[4] eq '+'){
		if ($sl[1] <= $sl[2]){
			++$same_pos_R1_lower;
			$dist_cutoff >= $dist and ++$same_pos_R1_lower_cut and $do_blast = 1;
		}else{
			++$same_pos_R1_higher;	
			$dist_cutoff >= $dist and ++$same_pos_R1_higher_cut and $do_blast = 1;
		}
	}else{
		if ($sl[3] eq '+'){
			if ($sl[1] >= $sl[2]){
				++$outward;
				$dist_cutoff >= $dist and ++$outward_cut and $do_blast = 1;
			}else{
				++$inward_too_far;
				$dist_cutoff >= $dist and ++$inward_too_far_cut;
			}
		}else{
			if ($sl[1] > $sl[2]){
				++$inward_too_far;
				$dist_cutoff >= $dist and ++$inward_too_far_cut;
			}else{
				++$outward;
				$dist_cutoff >= $dist and ++$outward_cut and $do_blast = 1;
			}			
		}
	}

	#run blast search for discordant reads, report status and whether this any sequence overlap ("microhomology") for broken read mapping points.	
	if ($do_blast){
		unless (exists($hit_pos_HoH{$sl[1]}->{$sl[2]})){
			$hit_pos_HoH{$sl[1]}->{$sl[2]} = 1;
			my $FHO = open_output("TEMP_BLAST_QUERY.txt");
			print $FHO ">read1\n$seq_HoA{$sl[0]}[0]\n>read2\n$seq_HoA{$sl[0]}[1]\n";	
			system("blastn -task blastn -evalue 0.000001 -db $blast_db -query TEMP_BLAST_QUERY.txt -out TEMP_BLAST_OUT.txt");
			
			
			my $SearchIO_obj = new Bio::SearchIO(-format => 'blast', -file   => "TEMP_BLAST_OUT.txt");
			my $result_obj1 = $SearchIO_obj->next_result;
			my $result_obj2 = $SearchIO_obj->next_result;
			
			my ($read1_status, $read1_overlap) = check_for_breaks($result_obj1);
			my ($read2_status, $read2_overlap) = check_for_breaks($result_obj2);
			
			#print output to breaks/overlap file
			print $FH_SUMM "$_\t$read1_status\t$read1_overlap\t$read2_status\t$read2_overlap\n";
			
			close $FHO;
			unlink ("TEMP_BLAST_QUERY.txt");
			unlink ("TEMP_BLAST_OUT.txt");
		}
	}
}

#print summary counts
print "File\tsame_neg_R1_lower\tsame_neg_R1_higher\tsame_pos_R1_lower\tsame_pos_R1_higher\toutward\tinward_too_far\tsame_neg_R1_lower_cut\tsame_neg_R1_higher_cut\tsame_pos_R1_lower_cut\tsame_pos_R1_higher_cut\toutward_cut\tinward_too_far_cut\n";
print "$file\t$same_neg_R1_lower\t$same_neg_R1_higher\t$same_pos_R1_lower\t$same_pos_R1_higher\t$outward\t$inward_too_far\t$same_neg_R1_lower_cut\t$same_neg_R1_higher_cut\t$same_pos_R1_lower_cut\t$same_pos_R1_higher_cut\t$outward_cut\t$inward_too_far_cut\n";


#subroutine to parse blast output and classify read based on mapping completeness and/or fragment location and report whether break points have any overlap ("microhomology")
sub check_for_breaks {

	use strict;
	use warnings;

	my $max_dist = 2000; ## NOTICE HARDCODING HERE

    my ($result_obj) = @_;
	my $broken = 0; #full-length mapping in two different pieces
	my $full = 0; #read maps end-to-end in one locations
	my $complex = 0; #complex mapping involving multiple pieces that is not a single simple break
	my $partial = 0; #Incomplete mapping. Only one hit but not full length
	my $no_hit = 0; #No blast hit for read.
	my $multiple_hsps = 0; #counter to track whether there is >1 hsps. Then used in classification purposes.
	my $overlap = "NA"; #variable that will store overlap length if there is a "broken" read.


	my $query_length = $result_obj->query_length;
	if (my $hit_obj = $result_obj->next_hit){
		my $first_hsp = 1;
		my $start_break = 0;
		my $end_break = 0;
		my $genome_coord = 0;
		while (my $hsp_obj = $hit_obj->next_hsp){
			my $hsp_length = $hsp_obj->length('query');
			$hsp_length == $query_length and $full = 1 and last;
						
			if ($first_hsp){
				if ($hsp_obj->start('query') > 1 and $hsp_obj->end('query') < $query_length){
					$complex = 1;
					last;
				}elsif ($hsp_obj->start('query') > 1){
					$start_break = $hsp_obj->start('query');
					$genome_coord = $hsp_obj->start('hit')
				}elsif($hsp_obj->end('query') < $query_length){
					$end_break = $hsp_obj->end('query');
					$genome_coord = $hsp_obj->end('hit')
				}else{
					die ("\nSomething went wrong in parsing position possibilities for first hsp\n\n");
				}	
			}else{
				$multiple_hsps = 1;
				if($start_break){
					if ($hsp_obj->start('query') == 1 and $hsp_obj->end('query') >= $start_break - 1){
						if (abs($genome_coord - $hsp_obj->end('hit')) < $max_dist){
							$overlap = $hsp_obj->end('query') - $start_break + 1;
							$broken = 1;
							last;
						}
					}
				}elsif($end_break){
					if ($hsp_obj->start('query') <= $end_break + 1 and $hsp_obj->end('query') == $query_length){
						if (abs($genome_coord - $hsp_obj->start('hit')) < $max_dist){
							$overlap = $end_break - $hsp_obj->start('query') + 1;
							$broken = 1;
							last;
						}
					}			
				}else{
					die ("\nSomething went wrong in assigning break positions\n\n")
				}
			}
			$first_hsp = 0;
		}
	}else{
		$no_hit = 1;
	}
	
	unless ($broken or $full or $complex or $no_hit){
		if ($multiple_hsps){
			$complex = 1;
		}else{
			$partial = 1;
		}
	}
	$broken + $complex + $full + $partial + $no_hit > 1 and die ("\nMultiple states returned for $_\n\n$broken + $complex + $full + $partial + $no_hit\n\n");

	my $status = "NA";
	
	if ($broken){
		$status = "broken"
	}elsif ($full){
		$status = "full"
	}elsif ($complex){
		$status = "complex"
	}elsif ($partial){
		$status = "partial"
	}elsif ($no_hit){
		$status = "no_hit"
	}
	
	return ($status, $overlap);
}

exit;
