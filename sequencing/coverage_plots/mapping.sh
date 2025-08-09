#!/bin/bash

#trim reads from Illumina fastq files for quality and adapters

for file in *1.fq.gz; do cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 20 --minimum-length 75 --cores 48 -e 0.15 -o ${file%fq.gz}trim.fq -p ${file%1.fq.gz}2.trim.fq $file ${file%1.fq.gz}2.fq.gz > ${file%_1.fq.gz}.cutadapt_log.txt; done


#map reads to Arabidopsis genomes. TAIR 10 assembly with orgaenlle genomes and the large numt removed from Chr2. Only one copy of IR included in plastid genome reference

for file in *1.trim.fq; do bowtie2 --no-unal -p 48 -X 1000 -x full_genome_map/gus_ref.no_chr2.fas -1 $file -2 ${file%1.trim.fq}2.trim.fq -S full_genome_map/${file%_1.trim.fq}.sam  >> full_genome_map/${file%_1.trim.45.fq}.bowtie2.log.txt 2>> full_genome_map/${file%_1.trim.45.fq}.bowtie2.err.txt; done

#summarize depth at each nucleotide position and in sliding windows across the organelle genome

for file in *sam; do samtools sort -o ${file%sam}bam $file; samtools index ${file%sam}bam; done

samtools depth -d 200000 -r mito -a -o depth.mito.txt *bam
samtools depth -d 200000 -r plastid_noIR -a -o depth.plastid.txt *bam

perl sliding_window_depth.cipro.pl depth.mito.txt 100 > depth.mito.sw.txt
perl sliding_window_depth.cipro.pl depth.plastid.txt 100 > depth.plastid.sw.txt

#manually modified headers to contain info on genotype, treatment and rep.
#depth.mito.sw.mod.txt
#depth.plastid.sw.mod.txt
