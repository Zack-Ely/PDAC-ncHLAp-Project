#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=END
#SBATCH --array=1-5


PAT_ID=$(cat /path/jacks/zacks_area/GeneFusion/Final_5_Organoid_list.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

export SINGULARITY_CACHEDIR=/path/jacks/zacks_area/GeneFusion/new_cache

module load singularity/3.5.0

cd /path/jacks/zacks_area/GeneFusion

mkdir "$PAT_ID"_fastqs

mkdir ./outputs/"$PAT_ID"_out

cp -r /path/jacks/zacks_area/intron_pipeline/create_proteome_search_spaces/nuORF_mutation/rnaSeq_analysis/new_star/"$PAT_ID"_*f*gz ./"$PAT_ID"_fastqs/


singularity exec --bind /path/jacks/zacks_area/GeneFusion/easyfuse_ref:/ref docker://tronbioinformatics/easyfuse:1.3.6 python /code/easyfuse/processing.py -i /path/jacks/zacks_area/GeneFusion/"$PAT_ID"_fastqs/*fast*gz -o /path/jacks/zacks_area/GeneFusion/outputs/"$PAT_ID"_out/


rm ./"$PAT_ID"_fastqs/*gz

