#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=END
#SBATCH --mail-user=zack_ely@mit.edu
#SBATCH --array=10-14

module add perl/5.24.1
export PATH="/net/gelbart/data/jacks/zacks_area/tools/tabix-0.2.6:$PATH";

module load bcftools/1.10.2

module load vcftools/0.1.13

module load singularity/3.5.0

PAT_ID=$(cat /net/gelbart/data/jacks/zacks_area/intron_pipeline/create_proteome_search_spaces/physical_sample_IDs.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

cp ../"$PAT_ID"/revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf.gz ./

gzip -d revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf.gz

grep -v "#" revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf | cut -f 1-10 > "$PAT_ID"_tail_g

grep "#" revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf | cut -f 1-10 > "$PAT_ID"_head_g

cat "$PAT_ID"_head_g "$PAT_ID"_tail_g > fix_revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf

bgzip -c fix_revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf > fix_revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf.gz

tabix -p vcf fix_revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf.gz

singularity exec pvactools-3.1.1.sif pvacseq generate_protein_fasta fix_revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf.gz 30 "$PAT_ID"_ONLY_germline_UNphased.fasta -d full

rm "$PAT_ID"_tail_g "$PAT_ID"_head_g
