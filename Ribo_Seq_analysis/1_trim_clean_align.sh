#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-type=END
#SBATCH --array=1-57

#path to STAR
PATH=/net/create_proteome_search_spaces/nuORF_mutation/rnaSeq_analysis/new_star/STAR-2.6.1d/bin/Linux_x86_64_static:$PATH

module load perl/5.24.1

module load r/3.6.0

module load rsem/1.3.1

FASTQ_ID=$(cat /chothani_atlas_paper/fastq_IDs.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

trimpath="/net/Ribo_seq_normal_tissue/trim_tool/Trimmomatic-0.36"

#trim adapters, etc
java -jar "$trimpath"/trimmomatic-0.36.jar SE -phred33 "$FASTQ_ID".fastq "$FASTQ_ID".trimmed.fastq ILLUMINACLIP:"$trimpath"/adapters/All_TruSeqForTrimmomatic.fa:2:30:10 MAXINFO:20:0.5 MINLEN:20


module load bowtie2/2.3.5.1
#separate contaminant RNA from useful RNA
bowtie2 -L 20 --un="$FASTQ_ID".cleaned.fastq -x contaminant_RNA/ContamRNA "$FASTQ_ID".trimmed.fastq > "$FASTQ_ID".contamRNA.fastq

#Align the filtered reads to genome w nuORF transcriptome annotation
module load star/2.5.3a

STAR --runMode alignReads \
	--runThreadN 8 \
	--genomeDir /Cabili_V26_genome_directory \
	--readFilesIn "$FASTQ_ID".cleaned.fastq \
	--outFileNamePrefix "$FASTQ_ID" \
	--alignSJDBoverhangMin 1 \
	--alignSJoverhangMin 51 \
	--outFilterMismatchNmax 2 \
	--alignEndsType EndToEnd \
	--alignIntronMin 20 \
	-–alignIntronMax 100000 \
	-–outFilterType BySJout \
	--outFilterMismatchNoverLmax 0.04 \
	--twopassMode Basic \
	--outSAMattributes MD NH
