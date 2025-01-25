library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

#The code below illustrates our process of deriving and plotting PSGMs on the scRNA-Seq dataset of 
#primary tumors depicted in Supp Fig 6. This process was applied identically to the other
#two human scRNA-seq datasets from the manuscript

#load organoid immunopeptidome data frame
data <- read_csv("./pdo_ip_df.csv")

#fxn to extract gene names from the 'entry_name' column
extract_gene_name <- function(entry_name) {
  gene_name <- sub(".*=([^ ]*).*", "\\1", entry_name)
  return(gene_name)
}

#filter out multi-mapping peptides and other patterns
data2 <- data %>%
  filter(species == "Human", !grepl("\\|", accession_numbers), length > 7, length < 12)

data2 <- data2 %>%
  filter(!grepl("_2D_", filename))

#relic
data3 <- data2 %>%
  filter(sapply(strsplit(accession_numbers, "\\|"), function(x) sum(grepl("^ENSP", x)) == 1))

#patterns to filter out
patterns_to_remove <- c("germ", "somat", "ENST", "Retained", "TCONS_", "ZNF", "T3", "smORFT0")

#relic
data4 <- data3 %>%
  filter(!grepl(paste(patterns_to_remove, collapse = "|"), accession_numbers))

#create a list of gene names for each sample (note, some redundant filtering code below)
gene_names_list <- data4 %>%
  filter(!grepl("\\|", accession_numbers)) %>%
  filter(!grepl(paste(patterns_to_remove, collapse = "|"), accession_numbers)) %>%
  filter(sapply(strsplit(accession_numbers, "\\|"), function(x) sum(grepl("^ENSP", x)) == 1)) %>%
  #remove rows where 'entry_name' does not contain 'GN='
  filter(grepl("GN=", entry_name)) %>%
  group_by(org_line) %>%
  summarize(gene_names = list(unique(extract_gene_name(entry_name))), .groups = 'drop')


gene_names_result <- setNames(gene_names_list$gene_names, gene_names_list$org_line)



##now assemble bulk tumor PSGMs
d625 <- read.csv("./C3L-00625_PDAC_HLA-I_BI_20220519_PSMexport.1.ssv", sep = ";")
d51 <- read.csv("./C3L-01051_PDAC_HLA-I_BI_20220519_PSMexport.1.ssv", sep = ";")
d31 <- read.csv("./C3L-01031_PDAC_HLA-I_BI_20220519_PSMexport.1.ssv", sep = ";")


d625 <- d625 %>%
  mutate(sample = "p0625")

d51 <- d51 %>%
  mutate(sample = "p1051")

d31 <- d31 %>%
  mutate(sample = "p1031")

merged_data <- bind_rows(d625, d51, d31)

bdata <- merged_data

#filter bulk tumor data frame
bdata <- bdata %>%
  filter(species == "Human")

#relic
bdata <- bdata %>%
  filter(!grepl("\\SAAP=", entry_name))

bdata2 <- bdata %>%
  filter(species == "Human", !grepl("\\|", accession_numbers), sequenceLength > 7, sequenceLength < 12)

#relic
bdata3 <- bdata2 %>%
  filter(sapply(strsplit(accession_numbers, "\\|"), function(x) sum(grepl("^ENSP", x)) == 1))

#patterns to filter (relic)
patterns_to_remove <- c("germ", "somat", "ENST", "Retained", "TCONS_", "ZNF", "T3", "smORFT0")

#relic
bdata4 <- bdata3 %>%
  filter(!grepl(paste(patterns_to_remove, collapse = "|"), accession_numbers))


JA_bulk_gene_names_list <- bdata4 %>%
  #remove rows where 'entry_name' doesnt contain 'GN='
  filter(grepl("GN=", entry_name)) %>%
  group_by(sample) %>%
  summarize(gene_names = list(unique(extract_gene_name(entry_name))), .groups = 'drop')

JA_bulk_gene_names_result <- setNames(JA_bulk_gene_names_list$gene_names, JA_bulk_gene_names_list$sample)



##################LOAD SEURAT OBJECT CONTAINING SINGLE CELL DATA OF PRIMARY PDAC TUMORS
#NOTE: this object was developed using all of the pre-processing steps described in the
#manuscript
pdacFP <- readRDS("./frp_combined_bc.RDS")

#make primary umap with annotated clusters
#this is after manual checks with feature plots using marker genes referenced in methods (see Peng et al)
new.cluster.ids <- c("Fibroblast","Malignant","Fibroblast","Macrophage","T cell","T cell","Endothelial","Ductal","B cell","Stellate","Fibroblast","Endocrine","Malignant","Mast cell","CDH19+ SCN7A+ cell")

names(new.cluster.ids) <- levels(pdacFP)
pdacFP_2 <- RenameIdents(pdacFP, new.cluster.ids)

#create Supplementary Figure 6A (left panel)
postscript("./FreedPastor_scData/1h_FP_annotated_umap.eps")
DimPlot(pdacFP_2, label = TRUE, label.size=5.5) + NoLegend()
dev.off()


####load PSG modules into the Seurat object - organoids, then bulk tumors
#as mentioned in the manuscript, some genes will be lost if they are unexpressed or 
#not present in the single cell data
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(gene_names_result[["PANFR0049"]]), ctrl = 8, name = "P049_Peptide_Source_Genes")
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(gene_names_result[["PANFR0069"]]), ctrl = 8, name = "P069_Peptide_Source_Genes")
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(gene_names_result[["PANFR0071"]]), ctrl = 8, name = "P071_Peptide_Source_Genes")
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(gene_names_result[["PANFR0087"]]), ctrl = 8, name = "P087_Peptide_Source_Genes")
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(gene_names_result[["PANFR0123"]]), ctrl = 8, name = "P123_Peptide_Source_Genes")
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(gene_names_result[["PANFR0151"]]), ctrl = 8, name = "P151_Peptide_Source_Genes")
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(gene_names_result[["PANFR0177"]]), ctrl = 8, name = "P177_Peptide_Source_Genes")
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(gene_names_result[["PANFR0290"]]), ctrl = 8, name = "P290_Peptide_Source_Genes")
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(gene_names_result[["PANFR0359"]]), ctrl = 8, name = "P359_Peptide_Source_Genes")
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(gene_names_result[["PANFR0402"]]), ctrl = 8, name = "P402_Peptide_Source_Genes")
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(gene_names_result[["PANFR0413"]]), ctrl = 8, name = "P413_Peptide_Source_Genes")


pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(JA_bulk_gene_names_result[["p0625"]]), ctrl = 8, name = "JA_OP2_P00625_Peptide_Source_Genes")
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(JA_bulk_gene_names_result[["p1031"]]), ctrl = 8, name = "JA_OP2_P01031_Peptide_Source_Genes")
pdacFP_2 <- AddModuleScore(object = pdacFP_2, features = list(JA_bulk_gene_names_result[["p1051"]]), ctrl = 8, name = "JA_OP2_P01051_Peptide_Source_Genes")


#output DotPlot - bulk + organoids
DotPlot(pdacFP_2, features=c("JA_OP2_P01051_Peptide_Source_Genes1","JA_OP2_P01031_Peptide_Source_Genes1","JA_OP2_P00625_Peptide_Source_Genes1","P049_Peptide_Source_Genes1","P069_Peptide_Source_Genes1","P071_Peptide_Source_Genes1","P087_Peptide_Source_Genes1","P123_Peptide_Source_Genes1","P151_Peptide_Source_Genes1","P177_Peptide_Source_Genes1","P290_Peptide_Source_Genes1","P359_Peptide_Source_Genes1","P402_Peptide_Source_Genes1","P413_Peptide_Source_Genes1"))+scale_colour_gradient2(low = "cadetblue1", mid="cadetblue3", high = "orangered4")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

#output 2 feature plots contrasting bulk and organoid modules in representative samples
#this covers the rest of 
#SFig 6A
postscript("./2025_OPTION2_JA_1h_FP_P0625_bulk_FeaturePlot.eps")
FeaturePlot(pdacFP_2, features = c("JA_OP2_P00625_Peptide_Source_Genes1"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1)+scale_colour_gradient(low = "lightcyan", high = "lightsalmon3")
dev.off()

postscript("./ULTIMATE_2025_FINAL_1h_FP_P413_organoid_FeaturePlot.eps")
FeaturePlot(pdacFP_2, features = c("P413_Peptide_Source_Genes1"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1)+scale_colour_gradient(low = "lightcyan", high = "lightsalmon3")
dev.off()







