library('biomaRt')
library(data.table)
library(plyr)
library(ggplot2)
library("RColorBrewer")
library(dplyr)
library("gplots")
library(estimate)
library(DESeq2)

###The below code starts off with counts matrices that were generated with the RNA-Seq analysis
###procedure (RSEM+STAR, etc) described in the methods of the manuscript.

#load bulk gene counts
counts_matrix = fread("./ENS_only_aggregated_bulk_gene_counts.txt")
classic_genes=fread("./classic_full_name_protein_coding_ensembl_ids.txt",header=FALSE)

#subset matrix by only  the classic genes
rownames(counts_matrix) <- counts_matrix$gene_id
v2_counts_matrix <- counts_matrix[rownames(counts_matrix) %in% classic_genes$V1, ]

#convert gene IDs to canonical gene symbols where possible
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- v2_counts_matrix$gene_id
genes2 <- gsub('\\..+$', '', genes)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes2,mart=mart)

v2_counts_matrix$gene_id <- genes2

#format matrix
df <- merge(x = G_list, y = v2_counts_matrix, by.x="ensembl_gene_id", by.y="gene_id")
df3 <- df[,3:51]
df4 <- as.matrix(df3)
rownames(df4) <- df$hgnc_symbol

#log-transform
gene.mean.exp <- rowMeans(log2(df4 + 1))

#next 2 commands jsut a plotting reference
hist(gene.mean.exp,50)
abline(v=2)

#filter down to expressed genes
exp.genes <- names(which(gene.mean.exp > 2))
df5 <- df4[rownames(df4) %in% exp.genes, ]

column_data <- data.frame(colnames(df5))
colnames(column_data) <- c("Bulk")


###now do all of the above for organoid data
org_counts_matrix = fread("./ENS_only_aggregated_organoid_gene_counts.txt")
#extract the 27 samples available for publication
True_PDAC_samples <- read.table("./True_PDAC_samples.txt")
True_PDAC_samples <- append(True_PDAC_samples$V1,"gene_id")
org_counts_matrix <- subset(org_counts_matrix, select = True_PDAC_samples)

#subset matrix by only  the classic genes
rownames(org_counts_matrix) <- org_counts_matrix$gene_id

#prepare matrix
v2_org_counts_matrix <- org_counts_matrix[rownames(org_counts_matrix) %in% classic_genes$V1, ]
v2_org_counts_matrix$gene_id <- genes2
org_df <- merge(x = G_list, y = v2_org_counts_matrix, by.x="ensembl_gene_id", by.y="gene_id")
org_df3 <- org_df[,3:29]
org_df4 <- as.matrix(org_df3)
rownames(org_df4) <- org_df$hgnc_symbol

org.gene.mean.exp <- rowMeans(log2(org_df4 + 1))

#plot
hist(org.gene.mean.exp,28)
abline(v=2)

#filter down to expressed genes
org.exp.genes <- names(which(org.gene.mean.exp > 2))
org_df5 <- org_df4[rownames(org_df4) %in% org.exp.genes, ]
org_column_data <- data.frame(colnames(org_df5))
colnames(org_column_data) <- c("Organoid")



###Now, combine organoid and bulk counts and normalize jointly
#in this case, do not exclude non-expressed organoid genes - these may include stromal genes that are informative 
#for the stromal comparison analysis
#subset org matrix by genes > 0 in bulk 

org_joint <- org_df4[rownames(org_df4) %in% rownames(df5), ]

#combine datasets
joint_matrix <- cbind(df5,org_joint)

#prep matrix
p1 <- rep(c("Bulk"),49)
p2 <- rep(c("Organoid"),27)
joint_column_data <- data.frame(c(p1,p2))
colnames(joint_column_data) <- c("Cohort")
rownames(joint_column_data) <- colnames(joint_matrix)

#load into Deseq for further data processing
joint_dds <- DESeqDataSetFromMatrix(countData = round(joint_matrix), colData = joint_column_data, design = ~ Cohort)

joint_dds <- estimateSizeFactors(joint_dds)

log_normalized_joint_counts <- log2(counts(joint_dds,normalize=T)+1)

write.table(log_normalized_joint_counts, file="./log_norm_counts_joint.txt", quote = F, sep = "\t",col.names = NA)



#Now perform analysis with ESTIMATE
#note: counts matrices (from above steps) loaded below were slightly edited (nano) in unix prior to re-loading in R
filterCommonGenes(input.f="./formatted_log_norm_counts_joint.txt",output.f="./lognorm_joint_PANFRgenes.gct",id="GeneSymbol")

#calculate and export all ESTIMATE (stromal, immune, etc) signature score and export for further analysis
#the resulting data are used to construct Fig 1C and SF1A
estimateScore("./lognorm_joint_PANFRgenes.gct", "./lognorm_joint_PANFR_estimate_score.gct", platform="illumina")


###Determine stromal and immune genes missing from the organoid matrix
#The two files below are the base ESTIMATE modules listed in the supplementary data of that publication
stromal_genes=fread("./stromal_gene_module.txt",header=FALSE)
immune_genes=fread("./immune_gene_module.txt",header=FALSE)


#Generate genes lists based on expression patterns and the ESTIMATE modules. These results then feed into
#SF1B-C.
write.table(setdiff(stromal_genes$V1,org.exp.genes), file="./stromal_genes_missing_in_organoids.txt")
write.table(setdiff(immune_genes$V1,org.exp.genes), file="./immune_genes_missing_in_organoids.txt")

write.table(intersect(stromal_genes$V1,org.exp.genes), file="./stromal_genes_expressed_in_organoids.txt")
write.table(intersect(immune_genes$V1,org.exp.genes), file="./immune_genes_expressed_in_organoids.txt")

##the same steps are then repeated for bulk tumors