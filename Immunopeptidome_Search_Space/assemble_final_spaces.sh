#create final search space file with all types
#of potential antigen sources

while read line ; do
	cat RENAMED_WRAPPED_NO_WT_"$line"_ONLY_germline_UNphased.fasta RENAMED_WRAPPED_NO_WT_"$line"_somatic_UNphased.fasta short_intron_"$line"_formatted.fasta RENAMED_WRAPPED_NO_WT_"$line"_ONLY_somatic_phased.fasta Final_GeneFusion_"$line"*.fasta "$line"_nuorf_mutation.fa "$line"_frame_1_protein.fa > ./2023_Complete_"$line"_search_space.fasta ;
done < ../physical_sample_IDs.txt
