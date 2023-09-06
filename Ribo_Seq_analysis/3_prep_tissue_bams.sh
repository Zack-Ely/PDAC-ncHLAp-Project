#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-type=END
#SBATCH --mail-user=zack_ely@mit.edu

module load perl/5.24.1

module load star/2.5.3a

module load samtools/1.10

#note this example is for one tissue. in this case, Brain_Bams is a list of all the ribo-seq BAM files processed for brain samples
samtools merge -b Brain_Bams.txt merged_Brain_tissues_NO_H.bam

samtools index merged_Brain_tissues_NO_H.bam
