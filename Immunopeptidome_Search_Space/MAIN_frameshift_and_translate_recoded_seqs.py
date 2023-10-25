# Import necessary packages

#!/usr/bin/python
import sys
import numpy as np
import subprocess
from Bio.Seq import Seq
import bisect
from Bio import SeqIO

# Function: manualTranslate
# Inputs: FASTA sequence to translate to protein
# Returns: List of protein sequences corresponding to specific FASTA sequence
# Summary: Uses manual codon-->AA dictionary to walk through and translate FASTA sequence. 
def manualTranslate(fastasequence):
        # Initialize codon table and list of peptides
        codontable = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

        # Translate full length fasta sequence all the way through
        fulllengthprotein = ''
        for i in xrange(0, len(fastasequence), 3):
                codon = fastasequence[i:i+3]
                # Account for bizarre edge cases that should really never happen
                if len(codon) != 3 or codon not in codontable:
                        break
                if codontable[codon] == '*':  # Stop translating when we hit a stop codon
                        break
                AA = codontable[codon]
                fulllengthprotein += AA
	
	return fulllengthprotein


fasta_file = sys.argv[1]

fasta_prefix = fasta_file.split('.')[0]

id_prefix = sys.argv[2]

sample_id = sys.argv[3]

outfile = open(fasta_prefix + '_frame_1_protein.fa','a')

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

#output original reading frame
for fasta in fasta_sequences:
	if fasta.id.split('|')[1] == ('-'):
		new_seq = str(fasta.seq.reverse_complement())
		peptide_seq = manualTranslate(new_seq)
		new_id = '>' + id_prefix + '_Frame_1' + '|' + 'Retained_intron_' + fasta.id.split('|')[0] + '|' + fasta.id.split('|')[1] + '|' + sample_id
		#new_id = '>' + fasta.id
		outfile.write(new_id+'\n')
		outfile.write(peptide_seq+'\n')
	else:
		new_seq = str(fasta.seq)
		peptide_seq = manualTranslate(new_seq)
		new_id = '>' + id_prefix + '_Frame_1' + '|' + 'Retained_intron_' + fasta.id.split('|')[0] + '|' + fasta.id.split('|')[1] + '|' + sample_id
		#new_id = '>' + fasta.id
		outfile.write(new_id+'\n')
		outfile.write(peptide_seq+'\n')


#output second reading frame
outfile_f2 = open(fasta_prefix + '_frame_2_protein.fa','a')

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

for fasta in fasta_sequences:
        if fasta.id.split('|')[1] == ('-'):
		new_seq = str(fasta.seq.reverse_complement())
                new_frame_seq = new_seq[1:]
		peptide_seq = manualTranslate(new_frame_seq)
                new_id = '>' + id_prefix + '_Frame_2' + '|' + 'Retained_intron_' + fasta.id.split('|')[0] + '|' + fasta.id.split('|')[1] + '|' + sample_id
		#new_id = '>' + fasta.id.split('|')[0] + '|' + fasta.id.split('|')[1] + '|' + id_prefix + '_Frame_2'
                outfile_f2.write(new_id+'\n')
                outfile_f2.write(peptide_seq+'\n')
	else:
		frame_seq = str(fasta.seq[1:])
             	new_seq = frame_seq
                peptide_seq = manualTranslate(new_seq)
                new_id = '>' + id_prefix + '_Frame_2' +	'|' + 'Retained_intron_' + fasta.id.split('|')[0] + '|' + fasta.id.split('|')[1] + '|' + sample_id
		#new_id = '>' + fasta.id.split('|')[0] + '|' + fasta.id.split('|')[1] + '|' + id_prefix + '_Frame_2'
                outfile_f2.write(new_id+'\n')
                outfile_f2.write(peptide_seq+'\n')

#output third reading frame
outfile_f3 = open(fasta_prefix + '_frame_3_protein.fa','a')

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

for fasta in fasta_sequences:
        if fasta.id.split('|')[1] == ('-'):
                new_seq = str(fasta.seq.reverse_complement())
                new_frame_seq = new_seq[2:]
                peptide_seq = manualTranslate(new_frame_seq)
                new_id = '>' + id_prefix + '_Frame_3' +	'|' + 'Retained_intron_' + fasta.id.split('|')[0] + '|' + fasta.id.split('|')[1] + '|' + sample_id
		#new_id = '>' + fasta.id.split('|')[0] + '|' + fasta.id.split('|')[1] + '|' + id_prefix + '_Frame_3'
                outfile_f3.write(new_id+'\n')
                outfile_f3.write(peptide_seq+'\n')
        else:
                frame_seq = str(fasta.seq[2:])
                new_seq = frame_seq
                peptide_seq = manualTranslate(new_seq)
		new_id = '>' + id_prefix + '_Frame_3' +	'|' + 'Retained_intron_' + fasta.id.split('|')[0] + '|' + fasta.id.split('|')[1] + '|' + sample_id
                #new_id = '>' + fasta.id.split('|')[0] + '|' + fasta.id.split('|')[1] + '|' + id_prefix + '_Frame_3'
                outfile_f3.write(new_id+'\n')
                outfile_f3.write(peptide_seq+'\n')
