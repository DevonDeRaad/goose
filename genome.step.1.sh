#!/bin/sh
#
#SBATCH --job-name=ref.align.gooses             # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --cpus-per-task=15               # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/work/goose	# Set working d$
#SBATCH --mem-per-cpu=5gb            # memory requested
#SBATCH --time=5000

#set the name of each read file 
rawfiles="CCBT1ANXX_s1_1_iTru7_11_02-iTru5_01_A_SL318134
CCBT1ANXX_s1_2_iTru7_11_02-iTru5_01_A_SL318134
ERR2193528_1
ERR2193528_2
ERR2193529_1
ERR2193529_2
ERR2193530_1
ERR2193530_2
ERR2193534_1
ERR2193534_2
ERR2193535_1
ERR2193535_2
ERR2193536_1
ERR2193536_2"


module load java
#directory that contains a folder with sample data in fasta format, plus the reference genome 'zebra.finch.ref.fa'
src=/home/d669d153/work/goose

#pre-process all read files to remove adapters and phix spike in
#for sample in $rawfiles
#do 
#/home/d669d153/work/bbmap/bbduk.sh in=${sample}.fastq.gz out=${sample}.clean.fastq.gz ref=/home/d669d153/work/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
#done

#for sample in $rawfiles
#do 
#/home/d669d153/work/bbmap/bbduk.sh in=${sample}.clean.fastq.gz out=${sample}.cleaned.fastq.gz ref=/home/d669d153/work/bbmap/resources/phix_adapters.fa.gz k=31 hdist=1 stats=${sample}.stats.txt
#done
#move on to alignment

#index reference genome using bwa so that we can align our rad reads to it
#/panfs/pfs.local/work/bi/bin/bwa/bwa index swan.goose.ref.fna


for sample in $rawfiles
do 
/panfs/pfs.local/work/bi/bin/bwa/bwa mem -M -t 15 -R '@RG\tID:${sample}\tSM:${sample}\tPL:Illumina' swan.goose.ref.fna ${sample}.cleaned.fastq.gz > ${sample}.sam
done

#combine mystery goose
/home/d669d153/work/gatk-4.1.4.0/gatk MergeSamFiles \
      -I CCBT1ANXX_s1_1_iTru7_11_02-iTru5_01_A_SL318134.sam \
      -I CCBT1ANXX_s1_2_iTru7_11_02-iTru5_01_A_SL318134.sam \
      -O mysterygoose.sam

#combine goose1
/home/d669d153/work/gatk-4.1.4.0/gatk MergeSamFiles \
      -I ERR2193528_1.sam \
      -I ERR2193528_2.sam \
      -I ERR2193529_1.sam \
      -I ERR2193529_2.sam \
      -I ERR2193530_1.sam \
      -I ERR2193530_2.sam \
      -O goose1.sam

#combine goose2
/home/d669d153/work/gatk-4.1.4.0/gatk MergeSamFiles \
      -I ERR2193534_1.sam \
      -I ERR2193534_2.sam \
      -I ERR2193535_1.sam \
      -I ERR2193535_2.sam \
      -I ERR2193536_1.sam \
      -I ERR2193536_2.sam \
      -O goose2.sam

files="goose1
goose2
mysterygoose"
      
#deduplicate and sort sam files using gatk
for sample in $files
do 
/home/d669d153/work/gatk-4.1.4.0/gatk MarkDuplicatesSpark \
      --remove-sequencing-duplicates true \
      -I ${sample}.sam \
      -O ${sample}.dedup.bam \
      --create-output-bam-index true
done

#if you need to free space
#rm *.sam
#this should give indexed, sorted, and de-duplicated, analysis-ready .bam files