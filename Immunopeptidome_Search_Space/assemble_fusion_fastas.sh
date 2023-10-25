#module load seqtk/1.3 

for i in {1..5} ; do
	cd /path/outputs
	zvar="$(head -n $i ../Final_5_Organoid_list.txt | tail -n 1)"
	xvar=$(echo "$zvar" | cut -d "_" -f 1)
	cd "$zvar"_out/FusionSummary
	cut -d ";" -f 3,20 "$zvar"_fusRank_1.pred.csv | grep -v "\;0" | sort -u -t \; -k 2 | grep -v "neo_pep" | sed 's/"//g' > part_1
	while read line ; do
		echo ">$line" | sed 's/;/\n/g' >> GeneFusion_"$zvar".fasta ;
	done < part_1
	rm part_1
	grep ">" GeneFusion_"$zvar".fasta > headerz
	d=$((0))
	while read line ; do
		d=$((d+1))
		var="$(grep "$line" GeneFusion_"$zvar".fasta | sed 's/:/_/g' | sed 's,>,,g')"
		var2=$(echo ">gene_fusion_$var|gene_fusion_$d|$xvar")
		sed -i '' -e "s,$line,$var2,g" GeneFusion_"$zvar".fasta ;
	done < headerz
	seqtk seq -L 8 GeneFusion_"$zvar".fasta > Final_GeneFusion_"$zvar".fasta
	rm headerz GeneFusion_"$zvar".fasta ;
done
