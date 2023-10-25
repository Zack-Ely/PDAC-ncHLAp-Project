#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-type=END
#SBATCH --array=10-14

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
export PATH="/path/zacks_area/tools/tabix-0.2.6:$PATH";
module load bcftools/1.10.2
module load fastxtoolkit/0.0.13


#this script IDs the phased seqeunces involving somatic variants that are truly distinct from the unphased somatic variant sequences.
#the output is a fasta file containing the resulting phased protein sequences
#this script does not retrieve distinct phased sequences only comprised of germline variants.

PAT_ID=$(cat /path/zacks_area/intron_pipeline/create_proteome_search_spaces/physical_sample_IDs.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

cd /path/zacks_area/intron_pipeline/create_proteome_search_spaces/mut_format/"$PAT_ID"

mv ../"$PAT_ID"_somatic_and_germline_phased.fasta ./

awk '/^>/ {P=index($0,">WT")==0} {if(P) print} ' "$PAT_ID"_somatic_and_germline_phased.fasta > no_wt_expt.fasta

fasta_formatter -i no_wt_expt.fasta -w 0 -o wrapped_no_wt_expt.fasta

grep ">" RENAMED_WRAPPED_NO_WT_"$PAT_ID"_somatic_UNphased.fasta | cut -d "_" -f 3 | sort | uniq > real_genes

grep -f real_genes somatic_headers > real_somatic_headers

grep -A 1 -f real_somatic_headers wrapped_no_wt_expt.fasta | grep -v "\-\-" > confined_phase_search

grep -v ">" confined_phase_search > only_confined_seqs


#if phased seqeunces have perfect match in unphase FASTA, then the phased sequence is not distinct. 
#thus, these somatic variants do not need a phased entry in the search space
grep -o -f only_confined_seqs RENAMED_WRAPPED_NO_WT_"$PAT_ID"_somatic_UNphased.fasta > seqs_phasing_doesnt_apply

#only unique sequences below will be the ones affected by phasing
cat seqs_phasing_doesnt_apply only_confined_seqs | sort | uniq -u > seqs_phasing_matters

#subset the phased fasta to only sequences where phasing matters
grep -B 1 -f seqs_phasing_matters wrapped_no_wt_expt.fasta | grep -v "\-\-" > prelim_phased_search_space_"$PAT_ID".fasta

###now must reformat heads

grep ">" prelim_phased_search_space_"$PAT_ID".fasta > headerz
cp prelim_phased_search_space_"$PAT_ID".fasta ./"$PAT_ID"_phased_all.fasta

while read line ; do
	x_count="$(grep -n "$line" headerz | cut -d ":" -f 1)"
	var="$(grep "$line" "$PAT_ID"_phased_all.fasta | cut -d "." -f 3- | tr '.' '_')"
	var2=$(echo ">somatic_phased_$var|entry_somatic_mutation_phased_$x_count|$PAT_ID" | sed 's,/,_,g')
	sed -i '' -e "s,$line,$var2,g" "$PAT_ID"_phased_all.fasta ;
done < headerz

rm headerz
mv "$PAT_ID"_phased_all.fasta ./RENAMED_WRAPPED_NO_WT_"$PAT_ID"_ONLY_somatic_phased.fasta

rm seqs_phasing_matters seqs_phasing_doesnt_apply only_confined_seqs confined_phase_search real_genes real_somatic_headers

rm no_wt_expt.fasta wrapped_no_wt_expt.fasta

rm prelim_phased_search_space_"$PAT_ID".fasta

#>somatic_unphased_CDK11B_ENST00000407249_3_missense_603K_T|somatic_unphased_mutation_1|PANFR0029
