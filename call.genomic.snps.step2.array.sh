#!/bin/sh
#
#SBATCH --job-name=ref.align.gooses             # Job Name
#SBATCH --nodes=1             # nodes
#SBATCH --cpus-per-task=1               # CPU allocation per Task
#SBATCH --array=79-96
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/temp/30day/bi/d669d153/goose	# Set working d$
#SBATCH --mem-per-cpu=3gb            # memory requested
#SBATCH --time=10000

module load java

#this step will take forever
#loop to generate a g.vcf for each bam file
/home/d669d153/work/gatk-4.1.4.0/gatk HaplotypeCaller \
   -R swan.goose.ref.fna \
   -I ERR38499${SLURM_ARRAY_TASK_ID}.dedup.bam \
   --emit-ref-confidence GVCF \
   -O ERR38499${SLURM_ARRAY_TASK_ID}.g.vcf
done

