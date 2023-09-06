head -n 333 HumanRNAContaminants.fa >> All_RNA_Contaminants.fa

head -n 1194 HumanRNAContaminants.fa | tail -n 861 >> All_RNA_Contaminants.fa

tail -n 49 HumanRNAContaminants.fa | head -n 46 >> All_RNA_Contaminants.fa

cat hg19-tRNAs.fa >> All_RNA_Contaminants.fa
