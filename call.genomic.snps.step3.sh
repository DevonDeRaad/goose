#!/bin/sh
#
#SBATCH --job-name=ref.align.gooses             # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --cpus-per-task=20               # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/temp/30day/bi/d669d153/goose	# Set working d$
#SBATCH --mem-per-cpu=3gb            # memory requested
#SBATCH --time=50000

#set the name of each read file 
rawfiles="mysterygoose
ERR3849979
ERR3849980
ERR3849981
ERR3849982
ERR3849983
ERR3849984
ERR3849985
ERR3849986
ERR3849987
ERR3849988
ERR3849989
ERR3849990
ERR3849991
ERR3849992
ERR3849993
ERR3849994
ERR3849995
ERR3849996
"

module load java

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

#for sample in $rawfiles
#do 
#/panfs/pfs.local/work/bi/bin/bwa/bwa mem -R "@RG\tID:${sample}\tSM:${sample}\tPL:Illumina" -t 20 -M swan.goose.ref.fna ${sample}_1.fastq.gz ${sample}_2.fastq.gz > ${sample}.sam
#done

#deduplicate and sort sam files using gatk
#index ref with picard
#/home/d669d153/work/gatk-4.1.4.0/gatk CreateSequenceDictionary -R swan.goose.ref.fna -O swan.goose.ref.dict

#for sample in $rawfiles
#do 
#/home/d669d153/work/gatk-4.1.4.0/gatk MarkDuplicatesSpark \
#      --remove-sequencing-duplicates true \
#      -R swan.goose.ref.fna \
#      -I ${sample}.sam \
#      -O ${sample}.dedup.bam \
#      --verbosity ERROR \
#      --create-output-bam-index true \
#      --conf 'spark.executor.cores=20'
#done

#index bams with samtools
#for sample in $rawfiles
#do
#/panfs/pfs.local/work/bi/bin/samtools-1.3.1/bin/samtools index ${sample}.dedup.bam
#done

#if you need to free space
#rm *.sam
#this should give indexed, sorted, and de-duplicated, analysis-ready .bam files

#index ref with samtools
#/panfs/pfs.local/work/bi/bin/samtools-1.3.1/bin/samtools faidx swan.goose.ref.fna
#index ref with picard
#/home/d669d153/work/gatk-4.1.4.0/gatk CreateSequenceDictionary -R swan.goose.ref.fna -O swan.goose.ref.dict


#AT THIS STEP, USE AN ARRAY TO RUN HAPLOTYPECALLER AS A SEPARATE JOB FOR EACH SAMPLE

#this step will take forever
#loop to generate a g.vcf for each bam file
#for sample in $rawfiles
#do 
#/home/d669d153/work/gatk-4.1.4.0/gatk HaplotypeCaller \
#   -R swan.goose.ref.fna \
#   -I ${sample}.dedup.bam \
#   --emit-ref-confidence GVCF \
#   -O ${sample}.g.vcf
#done

#NOW COME BACK AFTER ALL YOUR SAMPLES HAVE BEEN HAPLOTYPED, FINISH BY GENOTYPING, COMBINING INTO A SINGLE VCF, AND FILTERING

#Get the names of the vcf files to be used in the next step
ls -d -1 *.g.vcf > gvcf.list

# make .fai file 
samtools faidx swan.goose.ref.fna
# make .bed file
awk '{print $1 "\t0\t" $2}' swan.goose.ref.fna.fai > swan.goose.ref.fna.bed

#Genotyping with GVCF in all the variant files produced by HaplotypeCaller
#combine gvcfs
/home/d669d153/work/gatk-4.1.4.0/gatk GenomicsDBImport \
-V gvcf.list \
--genomicsdb-workspace-path my_database \
--reader-threads 20 \
--intervals swan.goose.ref.fna.bed

#genotype combined gvcfs and export single vcf
/home/d669d153/work/gatk-4.1.4.0/gatk GenotypeGVCFs \
    -R swan.goose.ref.fna \
    -V gendb://my_database \
    -O genotyped_X_samples.vcf

 #Extract the SNPs from the call set
/home/d669d153/work/gatk-4.1.4.0/gatk SelectVariants \
-R swan.goose.ref.fna \
-V genotyped_X_samples.vcf \
-select-type SNP \
-O genotyped_X_samples_snps.vcf

#Extract the indels from the call set
/home/d669d153/work/gatk-4.1.4.0/gatk SelectVariants \
-R swan.goose.ref.fna \
-V genotyped_X_samples.vcf \
-select-type INDEL \
-O genotyped_X_samples_indels.vcf

#filter SNP calls around indels and apply quality filters following Faircloth https://gist.github.com/brantfaircloth/4315737 and http://gatkforums.broadinstitute.org/discussion/3286/quality-score-recalibration-for-non-model-organisms   

/home/d669d153/work/gatk-4.1.4.0/gatk VariantFiltration \
-R swan.goose.ref.fna \
-V genotyped_X_samples_snps.vcf \
--mask genotyped_X_samples_indels.vcf \
--mask-extension 5 \
--mask-name InDel \
--cluster-window-size 10 \
--filter-expression "DP < 8" \
--filter-name "LowDepth" \
--filter-expression "QUAL < 30.0" \
--filter-name "LowQual" \
--filter-expression "QD < 5.0" \
--filter-name "LowVQCBD" \
--filter-expression "FS > 60.0" \
--filter-name "FisherStrand" \
--output genotyped_X_samples_filtered_1st.vcf

# get only pass snps in recal_0.filtered.vcf
cat genotyped_X_samples_filtered_1st.vcf | grep 'PASS\|^#' > filtered.goose.vcf

#final step: use vcftools to filter the vcf file to retain only SNPs present in all individuals
vcftools --vcf filtered.goose.vcf --max-missing-count 0 --recode --out complete.goose

vcftools --vcf complete.goose.recode.vcf --maf .05 --recode --out complete.goose.maf


