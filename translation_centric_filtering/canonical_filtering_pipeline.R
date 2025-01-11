library('biomaRt')
library(BiocManager)
library(BiocFileCache)
library(data.table)
library(tidyverse)
library(stringr)
library(rebus)
library(readxl)
library(BiocManager)


tidy_canonical_ip_data <- function(ip_data){
  tmp <- ip_data %>% 
    filter(ncmap == FALSE) %>% 
    mutate(sequence = toupper(sequence)) %>% 
    dplyr::select(sequence, hgnc_symbol) %>% 
    unique()
  
  tmp$num_org_lines <- "x"
  tmp$org_lines <- "x"
  tmp$alleles <- "x"
  
  pb = txtProgressBar(min = 0, max = nrow(tmp), initial = 0, style = 3) 
  for(i in seq(nrow(tmp))){
    seq <- tmp$sequence[i]
    tmp$num_org_lines[i] <- ip_data %>% 
      filter(seq_upper == seq) %>% 
      dplyr::select(org_line) %>% 
      unique() %>% 
      nrow()
    tmp$org_lines[i] <- ip_data %>% 
      filter(seq_upper == seq) %>% 
      dplyr::select(org_line) %>% 
      unique() %>% 
      as.list() %>% 
      .$org_line %>% 
      str_flatten(., ",")
    tmp$alleles[i] <- ip_data %>%
      filter(seq_upper == seq) %>% 
      dplyr::select(bestAllele) %>% 
      unique() %>% 
      as.list() %>% 
      .$bestAllele %>% 
      str_split(., ",") %>% 
      unlist() %>% 
      unique() %>% 
      str_flatten(., ",")
    setTxtProgressBar(pb,i)
    close(pb)
    
  }
  return(tmp)
}

read_ribo_output <- function(input_path){
  output <- read_csv(input_path, col_names = FALSE) %>% 
    mutate(peptide = str_extract(X1, one_or_more("[:ALPHA:]")),
           file = str_remove_all(str_extract(X1, "\\:" %R% one_or_more("[:GRAPH:]")), "\\:"),
           file = str_remove_all(file, pattern = "filtered_ribo_"),
           file = str_remove_all(file, pattern = ".txt")) %>% 
    dplyr::select(peptide, file) %>% 
    group_by(peptide) %>% 
    mutate(ribo_tissues = as.character(file),
           ribo_tissues = ifelse(length(ribo_tissues)>=2, 
                                 paste(ribo_tissues, collapse=","), ribo_tissues)) %>%
    dplyr::select(peptide, ribo_tissues) %>% 
    unique() %>% 
    mutate(num_ribo_tissues = str_count(ribo_tissues, pattern = boundary("word")))
  return(output)
}

read_final_peptide_data <- function(path, num_samples){
  library(tidyverse, quietly = TRUE)
  library(rebus, quietly = TRUE)
  library(stringr, quietly = TRUE)
  warning("\nATTACHING NECESSARY PACKAGES\n", call. = FALSE)
  dataset <- read_delim(path)
  
  dataset <- dataset %>%
    mutate(culture = ifelse(str_detect(filename, "2D"), "2D", "ORG"),
           ncmap = ifelse(species == "nuORFdb_human", TRUE, FALSE),
           ncmap = ifelse(str_detect(accession_number, pattern = START %R% "Retained_intron_"), TRUE, ncmap),
           ncmap_class = ifelse(ncmap == TRUE, str_remove(entry_name, pattern = START %R% one_or_more("[:GRAPH:]") %R% "\\s"), NA),
           ncmap_class = ifelse(str_detect(accession_number, pattern =  START %R% "Retained_intron_"), "Retained_intron", ncmap_class),
           ncmap = ifelse(str_detect(accession_numbers, pattern = "\\|") & ncmap_class != "Retained_intron", FALSE, ncmap),
           ncmap = ifelse(str_detect(species, pattern = "\\|"), FALSE, ncmap))
  return(dataset)
}

##### Read in IP data #####
final_peptide_data <- read_final_peptide_data("Partners HealthCare Dropbox/Zachary Kulstad/PDAC_Immunopeptidome_Project/2023_Manuscript_Folder/dataframes/final_df_withalleles.csv")

##### Filter out ncHLAp ##### 
org_data_with_pep_id <- final_peptide_data %>% 
  filter(ncmap == FALSE) %>% 
  mutate(ensembl_peptide_id = gsub('\\..+$', '', accession_number))

##### Load in mart dataset #####
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

##### Get gene names for Org IP #####
genes_org <- unique(org_data_with_pep_id$ensembl_peptide_id)
g_list_org <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id", "ensembl_gene_id", "hgnc_symbol"),values=genes_org,mart=mart)
g_list_org <- g_list_org %>% mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl_peptide_id, hgnc_symbol))
org_data_with_gene_id <- merge(x = g_list_org, y = org_data_with_pep_id, by.x="ensembl_peptide_id", by.y="ensembl_peptide_id")

##### Get gene names for HH IP #####
hh_data_with_pep_id <- healthy_human %>% 
  mutate(ensembl_peptide_id = gsub('\\..+$', '', accession_number)) %>% 
  filter(!str_detect(directory, pattern = "Testis_HLAI"))
genes_hh <- unique(hh_data_with_pep_id$ensembl_peptide_id)
g_list_hh <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id", "ensembl_gene_id", "hgnc_symbol"),values=genes_hh,mart=mart)
g_list_hh <- g_list_hh %>% mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl_peptide_id, hgnc_symbol))
hh_data_with_gene_id <- merge(x = g_list_hh, y = hh_data_with_pep_id, by.x="ensembl_peptide_id", by.y="ensembl_peptide_id")

##### Get gene names for Thymus IP #####
thy_data_with_pep_id <- thymus %>% 
  mutate(ensembl_peptide_id = gsub('\\..+$', '', accession_number))
genes_thy <- unique(thy_data_with_pep_id$ensembl_peptide_id)
g_list_thy <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id", "ensembl_gene_id", "hgnc_symbol"),values=genes_thy,mart=mart)
g_list_thy <- g_list_thy %>% mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl_peptide_id, hgnc_symbol))
thy_data_with_gene_id <- merge(x = g_list_thy, y = thy_data_with_pep_id, by.x="ensembl_peptide_id", by.y="ensembl_peptide_id")

##### Tidy Org IP dataframe #####
tidy_org_data_with_gene_id <- tidy_canonical_ip_data(org_data_with_gene_id)

##### Referencing back to HH and Thymus IP, appending to Org dataframe #####
tidy_org_data_with_gene_id <- tidy_org_data_with_gene_id %>% 
  mutate(hh_found = ifelse(hgnc_symbol %in% unique(hh_data_with_gene_id$hgnc_symbol), TRUE, FALSE),
         thymus_found = ifelse(hgnc_symbol %in% unique(thy_data_with_gene_id$hgnc_symbol), TRUE, FALSE)) 

##### Load in RiboTISH dataframe #####
all_tish <- read_csv("Dropbox (MIT)/PDAC_Immunopeptidome_Project/2023_Manuscript_Folder/Misc/files_for_ZK/all_tish_orfs.csv")

all_tish <- all_tish %>% 
  select(AASeq, tissue) %>% 
  unique()

##### Make empty dataframe #####
empty_matrix <- data.frame(sequence = tidy_org_data_with_gene_id$sequence,
                           brain = rep(0, nrow(tidy_org_data_with_gene_id)),
                           fat = rep(0, nrow(tidy_org_data_with_gene_id)),
                           fibroblast = rep(0, nrow(tidy_org_data_with_gene_id)),
                           ha_ec = rep(0, nrow(tidy_org_data_with_gene_id)),
                           hcaec = rep(0, nrow(tidy_org_data_with_gene_id)),
                           hepatocyte = rep(0, nrow(tidy_org_data_with_gene_id)),
                           kidney = rep(0, nrow(tidy_org_data_with_gene_id)),
                           vsmc = rep(0, nrow(tidy_org_data_with_gene_id)))

##### Iterate through Org IP dataframe, append to empty dataframe #####
for (i in seq(nrow(tidy_org_data_with_gene_id))) {
  print(i)
  progress_bar = txtProgressBar(min=1, max=nrow(tidy_org_data_with_gene_id), style = 3, char="=")
  count_output <- all_tish %>% 
    filter(str_detect(AASeq, tidy_org_data_with_gene_id$sequence[i])) %>% 
    group_by(tissue) %>% 
    count()
  if (!identical(integer(0), count_output$n)) {
    for (j in seq(nrow(count_output))){
      empty_matrix[i, count_output[j,1]$tissue] <- count_output[j,2]$n
    } 
  }
  setTxtProgressBar(progress_bar, value = i)
}

##### Rename empty dataframe #####
ribo_orf_counts <- read_csv("Downloads/ribo_counts.csv")

##### Determine number of tissues each peptide is detected in #####
numbers <- ribo_orf_counts %>% 
  pivot_longer(cols = -c(1)) %>% 
  filter(value != 0) %>% 
  group_by(sequence) %>% 
  count() %>% 
  select(sequence, num_ribo_tissues = n)


##### Create final DF with all data sources #####
canonical_tcfp <- left_join(tidy_org_data_with_gene_id, ribo_orf_counts, by = "sequence") %>% 
  mutate(total = brain + fat + fibroblast + ha_ec + hcaec + hepatocyte + kidney + vsmc) %>% 
  left_join(., numbers, by = "sequence") %>% 
  select(sequence, 
         num_org_lines, 
         num_ribo_tissues,
         thymus_found,
         total_ribo_orfs = total, 
         everything()) %>%
  select(-c(brain,fat,fibroblast,ha_ec,hcaec,hepatocyte,vsmc, kidney)) %>% 
  mutate(num_ribo_tissues = ifelse(is.na(num_ribo_tissues), 0, num_ribo_tissues)) %>% 
  mutate(cancer_restricted = ifelse(hh_found == FALSE & thymus_found == FALSE & total_ribo_orfs == 0, "CR", "nonCR"))
