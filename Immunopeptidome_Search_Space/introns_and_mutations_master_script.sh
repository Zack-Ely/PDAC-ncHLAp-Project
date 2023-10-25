#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-type=END
#SBATCH --mail-user=zack_ely@mit.edu
#SBATCH --array=10-14

export PATH="/path/intron_pipeline/myseq:$PATH";
module load python/2.7.13
module load r/4.0.4
module load bowtie2/2.3.5.1 
module load samtools/1.10
module load emboss/6.6.0
module load picard/2.18.26
module load kallisto/0.45.0
module load ucsc-tools/20170321

module load gatk/4.1.2.0

module add perl/5.24.1
module load vcftools/0.1.13
export PATH="/path/tools/tabix-0.2.6:$PATH";
module load bcftools/1.10.2

#this script is for generating a search space of patient-specific germline and somatic (phased) variants) and retained introns
#note this script assumes certain processes described in the ms have already been run, for example intron retention predction with
#KMA/IRFinder and somatic mutation calling
#please reach out to authors if you have questions or need supplementary scripts

PAT_ID=$(cat /path/intron_pipeline/create_proteome_search_spaces/physical_sample_IDs.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

mkdir /path/intron_pipeline/create_proteome_search_spaces/"$PAT_ID"


######call germline variants using Strelka
cd /path/OrganoidWGS/OrganoidWGS/
n_var="$(ls | grep "$PAT_ID"_Ger | head -n 1)"
var_norm="$(echo "$n_var".bam)"
cd /path/intron_pipeline/create_proteome_search_spaces/"$PAT_ID"
mkdir results_strelka

#prep strelka run directory
python /path/tools/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
--bam /path/OrganoidWGS/OrganoidWGS/"$n_var"/*/"$var_norm" \
--referenceFasta /path/another_grch37/CORRECT_pancseq_reference/Homo_sapiens_assembly19.fasta \
--runDir /path/intron_pipeline/create_proteome_search_spaces/"$PAT_ID"/results_strelka/
python /path/intron_pipeline/create_proteome_search_spaces/"$PAT_ID"/results_strelka/runWorkflow.py -m local -j 8

#get pass only
gzip -d /path/intron_pipeline/create_proteome_search_spaces/"$PAT_ID"/results_strelka/results/variants/variants.vcf.gz
#should probably use VCF recode here to get only pass variants

#grep -v "#" /path/intron_pipeline/create_proteome_search_spaces/"$PAT_ID"/results_strelka/results/variants/variants.vcf | grep PASS > "$PAT_ID"_germline_pass.vcf
grep -E "#|PASS" /path/intron_pipeline/create_proteome_search_spaces/"$PAT_ID"/results_strelka/results/variants/variants.vcf > "$PAT_ID"_germline_pass.vcf

#rewrite sample name to match somatic vcf file
sed -i 's/_Germline//g' "$PAT_ID"_germline_pass.vcf


######retrieve and prep the somatic variant vcf file
cp /path/fourth_organoid_analysis/merge_and_annotate_vcf_files/merged_union_indel_snvs_vcfs/"$PAT_ID"_merged.vcf ./

#remove non-canonical contig/chromosome references and the 'normal sample' column from the VCF (to make merging with germline possible)
grep -v "contig=" "$PAT_ID"_merged.vcf > v2_"$PAT_ID"_merged.vcf
grep "##" v2_"$PAT_ID"_merged.vcf > top
grep "#" v2_"$PAT_ID"_merged.vcf | tail -n 1 | cut -f 1-9,11 > top_2
grep -v "#" v2_"$PAT_ID"_merged.vcf | cut -f 1-9,11 > bottom

cat top top_2 bottom > v3_"$PAT_ID"_merged.vcf

rm top top_2 bottom v2_"$PAT_ID"_merged.vcf "$PAT_ID"_merged.vcf

sed -i "s/TUMOR/$PAT_ID/g" v3_"$PAT_ID"_merged.vcf


######Merge the germline and somatic variants; note will retrieve somatic indels at a later step
gatk MergeVcfs \
-R /path/another_grch37/CORRECT_pancseq_reference/Homo_sapiens_assembly19.fasta \
-I "$PAT_ID"_germline_pass.vcf \
-I v3_"$PAT_ID"_merged.vcf \
-D /path/intron_pipeline/analysis/pancreas_T_N_organoid_IR_comparison/analysis_determine_pt_specific_RIs/pilot_for_P413/phase_work/Homo_sapiens_assembly19.dict \
-O "$PAT_ID"_combined_somatic_plus_germline.vcf

#sort the combined variants
java -jar /path/picard_jar/picard.jar SortVcf \
I="$PAT_ID"_combined_somatic_plus_germline.vcf \
O=sorted_"$PAT_ID"_combined_somatic_plus_germline.vcf \
SEQUENCE_DICTIONARY=/path/intron_pipeline/analysis/pancreas_T_N_organoid_IR_comparison/analysis_determine_pt_specific_RIs/pilot_for_P413/phase_work/Homo_sapiens_assembly19.dict

######Phase the Germline-only and combined sets of variants
module unload gatk/4.1.2.0
module load gatk/3.7


#change the column headers to match the name of the bam files being used for phasing (otherwise cmd will fail)
bam_name_somatic="$(ls /path/OrganoidWGS/OrganoidWGS/ | grep $PAT_ID | grep -v Germ)"

sed -i "s/$PAT_ID/$bam_name_somatic/g" sorted_"$PAT_ID"_combined_somatic_plus_germline.vcf
sed -i "s/$PAT_ID/"$PAT_ID"_Germline/g" "$PAT_ID"_germline_pass.vcf


java -jar /home/software/gatk/gatk-3.7/GenomeAnalysisTK.jar -T ReadBackedPhasing \
-R /path/another_grch37/CORRECT_pancseq_reference/Homo_sapiens_assembly19.fasta \
-I /path/OrganoidWGS/OrganoidWGS/"$PAT_ID"_T*/*/*.bam  \
--variant sorted_"$PAT_ID"_combined_somatic_plus_germline.vcf \
-L sorted_"$PAT_ID"_combined_somatic_plus_germline.vcf \
-o phased_"$PAT_ID"_combined_somatic_plus_germline.vcf


java -jar /home/software/gatk/gatk-3.7/GenomeAnalysisTK.jar -T ReadBackedPhasing \
-R /path/another_grch37/CORRECT_pancseq_reference/Homo_sapiens_assembly19.fasta \
-I /path/OrganoidWGS/OrganoidWGS/"$PAT_ID"_Germline/*/"$PAT_ID"_Germline.bam  \
--variant "$PAT_ID"_germline_pass.vcf \
-L "$PAT_ID"_germline_pass.vcf \
-o phased_"$PAT_ID"_ONLY_germline_variants.vcf


######Annotate the combined set and germline-only set of variants using Ensembl VEP
/path/tools/ensembl-vep/vep --input_file phased_"$PAT_ID"_combined_somatic_plus_germline.vcf --output_file FINAL_FS_Plug_phased_"$PAT_ID"_somatic_and_germline_variants_annotated.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta /path/another_grch37/GRCh37.p13.genome.fa --offline --cache --plugin Frameshift --plugin Wildtype --dir_cache /path/ensembl_vep_cache_37 --dir_plugins /path/5_apr_att/VEP_plugins --transcript_version --pick --pick_order rank,canonical,appris,tsl,biotype,ccds,length,mane --force_overwrite

/path/tools/ensembl-vep/vep --input_file phased_"$PAT_ID"_ONLY_germline_variants.vcf --output_file FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta /path/another_grch37/GRCh37.p13.genome.fa --offline --cache --plugin Frameshift --plugin Wildtype --dir_cache /path/ensembl_vep_cache_37 --dir_plugins /path/5_apr_att/VEP_plugins --transcript_version --pick --pick_order rank,canonical,appris,tsl,biotype,ccds,length,mane --force_overwrite

##Retrieve and annotate missing somatic indels as well
cp /path/fourth_organoid_analysis/merge_and_annotate_vcf_files/union_indel_vcfs/"$PAT_ID"_union.vcf ./
/path/tools/ensembl-vep/vep --input_file "$PAT_ID"_union.vcf --output_file FINAL_FS_Plug_"$PAT_ID"_union_annotated.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta /path/another_grch37/GRCh37.p13.genome.fa --offline --cache --plugin Frameshift --plugin Wildtype --dir_cache /path/ensembl_vep_cache_37 --dir_plugins /path/5_apr_att/VEP_plugins --transcript_version --pick --pick_order rank,canonical,appris,tsl,biotype,ccds,length,mane --force_overwrite



######REVISE variants by excluding germline indels and recalling certain somatic variants that were mis-called by Mutect

cat FINAL_FS_Plug_phased_"$PAT_ID"_somatic_and_germline_variants_annotated.vcf > temp.vcf
grep "#" temp.vcf > top
#get all PASSed variants except germline indels and certain genes to be addressed next

#grep -v -E "#|CIGAR" temp.vcf | grep "PASS" | grep -v -E "\|KRAS\||\|TP53\||\|GNAS\||\|RNF43\||\|PLEC\||\|FLG\||\|AHNAK\||\|APOB\||\|CSMD1\||\|PLXNA1\||\|MCM6\||\|MKI67\||\|SIPA1\|" > m1

#previously, using 'grep PASS' caused many artifactual variants to be retained because they had the amino acid sequence, 'PASS'
#refer to Freed-Pastor et al 2021 cancer cell for explanation
grep -v -E "#|CIGAR" temp.vcf | awk -F '\t' '$7 == "PASS"' | grep -v -E "\|KRAS\||\|TP53\||\|GNAS\||\|RNF43\||\|PLEC\||\|FLG\||\|AHNAK\||\|APOB\||\|CSMD1\||\|PLXNA1\||\|MCM6\||\|MKI67\||\|SIPA1\|" > m1
grep "|KRAS|" temp.vcf | grep -E "PASS|clustered_event|homologous_mapping_event|panel_of_norm" > m2
grep "|TP53|" temp.vcf | grep -E "PASS|panel_of_norm|bSeq" > m3
grep "|GNAS|" temp.vcf | grep -E "PASS|clustered_events|panel_of_norm" > m4
grep "|RNF43|" temp.vcf | grep -E "PASS|clustered_events" > m5
grep "|PLEC|" temp.vcf | grep -E "PASS|clustered_events|panel_of_norm" > m6
grep "|FLG|" temp.vcf | grep -E "PASS|panel_of_norm" > m7
grep "|AHNAK|" temp.vcf | grep -E "PASS|panel_of_norm" > m8
grep "|APOB|" temp.vcf | grep -E "PASS|panel_of_norm" > m9
grep "|CSMD1|" temp.vcf | grep -E "PASS|panel_of_norm" > m10
grep "|PLXNA1|" temp.vcf | grep -E "PASS|clustered_events" > m11
grep "|MCM6|" temp.vcf | grep -E "PASS|clustered_events" > m12
grep "|MKI67|" temp.vcf | grep -E "PASS|clustered_events" > m13
grep "|SIPA1|" temp.vcf | grep -E "PASS|clustered_events|homologous_mapping" > m14
rm temp.vcf
cat top m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 > temp.vcf
vcf-sort -c temp.vcf > revised_FINAL_FS_Plug_phased_"$PAT_ID"_somatic_and_germline_variants_annotated.vcf
rm temp.vcf top m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14


cat FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf > tempz.vcf
grep "#" tempz.vcf > topz
#fixed the below by using the Awk ... PASS function
grep -v -E "#|CIGAR" tempz.vcf | awk -F '\t' '$7 == "PASS"' > mz1
rm tempz.vcf
cat topz mz1 > tempz.vcf
vcf-sort -c tempz.vcf > revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf
rm mz1 topz tempz.vcf



######add in somatic indels and zip and index both vcf files
#output somatic indels for inclusion
grep -v "#" FINAL_FS_Plug_"$PAT_ID"_union_annotated.vcf > to_add

cat revised_FINAL_FS_Plug_phased_"$PAT_ID"_somatic_and_germline_variants_annotated.vcf to_add > temp_added.vcf

#keep obsolete vcf in case helpful for debugging later
mv revised_FINAL_FS_Plug_phased_"$PAT_ID"_somatic_and_germline_variants_annotated.vcf ./obsolete.vcf

vcf-sort -c temp_added.vcf > revised_FINAL_FS_Plug_phased_"$PAT_ID"_somatic_and_germline_variants_annotated.vcf

bgzip -c revised_FINAL_FS_Plug_phased_"$PAT_ID"_somatic_and_germline_variants_annotated.vcf > revised_FINAL_FS_Plug_phased_"$PAT_ID"_somatic_and_germline_variants_annotated.vcf.gz
tabix -p vcf revised_FINAL_FS_Plug_phased_"$PAT_ID"_somatic_and_germline_variants_annotated.vcf.gz

rm temp_added.vcf to_add


#now zip and index gerrmline
bgzip -c revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf > revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf.gz
tabix -p vcf revised_FINAL_FS_Plug_phased_"$PAT_ID"_ONLY_germline_variants_annotated.vcf.gz



######run haplotype-maker script to perform intron-recoding
#retrieve the PASSed SNP variants' coordinates (Germline and somatic SNPs!)


####Excluding germline and somatic indels from vcf file
grep -v -E "CIGAR=|SomaticEVS" revised_FINAL_FS_Plug_phased_"$PAT_ID"_somatic_and_germline_variants_annotated.vcf | grep -E "#|PASS" > temp_var.vcf

#obsolete: grep -E "#|SNVHPO" /path/intron_pipeline/create_proteome_search_spaces/"$PAT_ID"/results_strelka/results/variants/variants.vcf | grep -E "#|PASS" > temp_var.vcf

grep -v "#" temp_var.vcf | cut -f 1-3 > PASS_coordinates.bed
cut -f 1-2 PASS_coordinates.bed > part_1.txt

cut -f 2 PASS_coordinates.bed > part_2.txt

paste part_1.txt part_2.txt | bedtools sort -i stdin > tmpf
rm PASS_coordinates.bed part_1.txt part_2.txt
mv tmpf ./PASS_coordinates.bed



#Reformat retained intron coordinates as BED files
#extract coordinates of P413 introns
cp /path/intron_pipeline/organoid_outs/"$PAT_ID"_*/headermap.txt ./
cat headermap.txt | cut -d "|" -f 1 | cut -f 2 | cut -d ">" -f 2 > "$PAT_ID"_introns.txt

#reformat to BED format
sed 's/:/\t/g' "$PAT_ID"_introns.txt | sed 's/-/\t/g' > "$PAT_ID"_introns.bed

#sort bed
bedtools sort -i "$PAT_ID"_introns.bed > tmpf
rm "$PAT_ID"_introns.bed
mv tmpf ./"$PAT_ID"_introns.bed
sed -i 's/chr//g' "$PAT_ID"_introns.bed

#subset to all germline and somatic SNPs that overlap all retained introns detected in the sample
bedtools intersect -wa -a PASS_coordinates.bed -b "$PAT_ID"_introns.bed -sorted | cut -f 1-2 > key_ALL_intron_variant_coords.txt



grep -f key_ALL_intron_variant_coords.txt revised_FINAL_FS_Plug_phased_"$PAT_ID"_somatic_and_germline_variants_annotated.vcf > retained_variants
grep "#" revised_FINAL_FS_Plug_phased_"$PAT_ID"_somatic_and_germline_variants_annotated.vcf > retained_header
cat retained_header retained_variants > retained_intron_variants.vcf
rm retained_header retained_variants

bgzip -c retained_intron_variants.vcf > retained_intron_variants.vcf.gz
tabix -p vcf retained_intron_variants.vcf.gz



######Retrieve retained intron sequences and recode with germline SNPs
cp /path/intron_pipeline/organoid_outs/"$PAT_ID"_*/kma_results.txt ./
python /path/intron_pipeline/create_proteome_search_spaces/v7_kmaToPeptideSeqs.py ./kma_results.txt 9 ./

sed -i 's/chr//g' CodingIntronSeqsFASTA.txt

bcftools consensus -H A -f CodingIntronSeqsFASTA.txt retained_intron_variants.vcf.gz > recoded_introns_all_ALT.fasta

bcftools consensus -H 1 -f CodingIntronSeqsFASTA.txt retained_intron_variants.vcf.gz > recoded_introns_H1.fasta

bcftools consensus -H 2 -f CodingIntronSeqsFASTA.txt retained_intron_variants.vcf.gz > recoded_introns_H2.fasta


cut -f 1 Intronheadermap.txt > replace_part_1

sed 's/|\-/|\-|Reference_Frame_1/g' replace_part_1 > ref_frame_1_part_1

sed -i 's/|+/|+|Reference_Frame_1/g' ref_frame_1_part_1

#should be able to just take this and alter the header to make 'Frame_2' etc after translating in python


sed 's/|\-/|\-|Alt_Variants_Frame_1/g' replace_part_1 > alt_frame_1_part_1

sed -i 's/|+/|+|Alt_Variants_Frame_1/g' alt_frame_1_part_1

sed 's/|\-/|\-|Haplotype_1_Frame_1/g' replace_part_1 > H1_frame_1_part_1

sed -i 's/|+/|+|Haplotype_1_Frame_1/g' H1_frame_1_part_1

sed 's/|\-/|\-|Haplotype_2_Frame_1/g' replace_part_1 > H2_frame_1_part_1

sed -i 's/|+/|+|Haplotype_2_Frame_1/g' H2_frame_1_part_1


grep ">" CodingIntronSeqsFASTA.txt > orig_headers

x=0

while read line; do
	let x=x+1
	repl="$(head -n $x ref_frame_1_part_1 | tail -n 1)"
	sed -i "s/$line/$repl/g" CodingIntronSeqsFASTA.txt ;
done < orig_headers


x=0

while read line; do
        let x=x+1
        repl="$(head -n $x alt_frame_1_part_1 | tail -n 1)"
        sed -i "s/$line/$repl/g" recoded_introns_all_ALT.fasta ;
done < orig_headers


x=0

while read line; do
        let x=x+1
        repl="$(head -n $x H1_frame_1_part_1 | tail -n 1)"
        sed -i "s/$line/$repl/g" recoded_introns_H1.fasta ;
done < orig_headers

x=0

while read line; do
        let x=x+1
        repl="$(head -n $x H2_frame_1_part_1 | tail -n 1)"
        sed -i "s/$line/$repl/g" recoded_introns_H2.fasta ;
done < orig_headers


rm orig_headers alt_frame_1_part_1 H1_frame_1_part_1 H2_frame_1_part_1 replace_part_1



######run re-frame/translation script to generate all possible retained intron-encoded protein sequences

python /path/intron_pipeline/create_proteome_search_spaces/MAIN_frameshift_and_translate_recoded_seqs.py CodingIntronSeqsFASTA.txt Haplotype_Reference "$PAT_ID"

python /path/intron_pipeline/create_proteome_search_spaces/MAIN_frameshift_and_translate_recoded_seqs.py recoded_introns_all_ALT.fasta Haplotype_All_Variants "$PAT_ID"

python /path/intron_pipeline/create_proteome_search_spaces/MAIN_frameshift_and_translate_recoded_seqs.py recoded_introns_H1.fasta Haplotype_H1 "$PAT_ID"

python /path/intron_pipeline/create_proteome_search_spaces/MAIN_frameshift_and_translate_recoded_seqs.py recoded_introns_H2.fasta Haplotype_H2 "$PAT_ID"


#concatenate all retained intron fasta files

cat recoded*fa CodingIntronSeqsFASTA*fa > "$PAT_ID"_retained_intron_search_space.fa


#####Remove all entries with sequence <9 amino acids
