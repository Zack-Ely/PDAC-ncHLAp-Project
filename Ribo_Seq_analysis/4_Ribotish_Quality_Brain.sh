#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=70G
#SBATCH --mail-type=END

module load perl/5.24.1

module load star/2.5.3a

module load r/4.0.4

module load python3/3.6.4

ribotish quality -b merged_Brain_tissues_NO_H.bam -g resolved_Cabili_only_stranded_V26_combined_gencode_mitranscriptome_unannotated.gtf 


