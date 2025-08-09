#!/bin/bash

#trim reads from Illumina fastq files for quality and adapters

for file in *1.fq.gz; do cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 20 --minimum-length 75 --cores 48 -e 0.15 -o ${file%fq.gz}trim.fq -p ${file%1.fq.gz}2.trim.fq $file ${file%1.fq.gz}2.fq.gz > ${file%_1.fq.gz}.cutadapt_log.txt; done


#truncate reads to the first 45 bp for mapping

for file in *trim.fq; do perl ../fastq_truncate.pl $file 45 > ${file%fq}45.fq; done

#map to plastid reference genome with only one copy of the inverted repeat

for file in *1.trim.45.fq; do bowtie2 --no-unal -p 48 -X 1000 -x cp_genome.noIR.fas -1 $file -2 ${file%1.trim.45.fq}2.trim.45.fq -S ${file%_1.trim.45.fq}.sam  >> ${file%_1.trim.45.fq}.bowtie2.log.txt 2>> ${file%_1.trim.45.fq}.bowtie2.err.txt; done


#parse discordant read pairs from sam files

for file in *sam; do perl sam_discordant_parse.pl $file 45 128214 84171 110434 600 > ${file%sam}discord.txt; done

for file in *discord.txt; do perl ../discord_summary.pl $file ${file%txt}breaks.txt ${file%.discord.txt}_1.trim.fq ${file%.discord.txt}_2.trim.fq cp_genome.noIR.fas 1001 127214 1000 >> discord.summary.txt; done

#analyze sequence overlaps at break points for "microhomologies"

for file in *breaks.txt; do perl ../overlap_hist.pl $file > ${file%breaks.txt}overlap_hist.txt; done

#simulate overlaps for a null/random expectation

perl ../overlap_sim.pl cp_genome.noIR.fas 1501 126714 160 150 1000 100000 > overlap_sim.out.txt 2> overlap_sim.err.txt

perl ../overlap_summary.pl hist_libraries.txt overlap_sim.out.txt 30 > overlap_freq.txt

#Generated Overlap_plot.pdf with Overlap_plot.R


#count read mapping data from bowtie stderr

for file in *err.txt; do perl ../bowtie_map_count.pl $file >> bowtie_counts.txt; done


#manually combined and modified discord.summary.txt and bowtie_counts.txt to discord.summary.mod.csv. Analyzed and plotted with Discord_plot_and_stats.R


#analyzing hotspots

perl ../discord_reads_by_pos.pl 50 128250 C4_CKDN230028844-1A_HGKLFDSX7_L1.discord.breaks.txt C8_CKDN230028848-1A_HGKLFDSX7_L1.discord.breaks.txt C12_CKDN230028852-1A_HGKLFDSX7_L1.discord.breaks.txt C20_CKDN230028856-1A_HGKLFDSX7_L1.discord.breaks.txt C24_CKDN230028860-1A_HGKLFDSX7_L1.discord.breaks.txt > hotspots.mutS2.txt

perl ../discord_reads_by_pos.pl 50 128250 C3_CKDN230028843-1A_HGKLFDSX7_L1.discord.breaks.txt C7_CKDN230028847-1A_HGKLFDSX7_L1.discord.breaks.txt C11_CKDN230028851-1A_HGKLFDSX7_L1.discord.breaks.txt C19_CKDN230028855-1A_HGKLFDSX7_L1.discord.breaks.txt C23_CKDN230028859-1A_HGKLFDSX7_L1.discord.breaks.txt > hotspots.WT.txt 

#Generated plot with Hotspot_plot.R (Hotspot_plot.pdf)

