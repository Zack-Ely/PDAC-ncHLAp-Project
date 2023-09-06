#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-type=END
#SBATCH --array=1-57

module load perl/5.24.1

module load star/2.5.3a

module load samtools/1.10

FASTQ_ID=$(cat /net/Ribo_seq_normal_tissue/chothani_atlas_paper/fastq_IDs.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

samtools view -S -b "$FASTQ_ID"Aligned.out.sam > "$FASTQ_ID"Aligned.out.bam

samtools sort "$FASTQ_ID"Aligned.out.bam -o "$FASTQ_ID"Aligned.out.sorted.bam

samtools index "$FASTQ_ID"Aligned.out.sorted.bam
