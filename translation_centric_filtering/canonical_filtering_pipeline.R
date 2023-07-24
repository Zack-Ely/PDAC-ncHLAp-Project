library('biomaRt')
library(data.table)
library(tidyverse)
library(stringr)
library(rebus)
library(readxl)

tidy_canonical_ip_data <- function(ip_data){
  tmp <- ip_data %>% 
    filter(ncmap == FALSE) %>% 
    mutate(sequence = str_to_upper(sequence)) %>% 
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

final_peptide_data <- read_final_peptide_data("PATH")

org_data_with_pep_id <- final_peptide_data %>% 
  filter(ncmap == FALSE) %>% 
  mutate(ensembl_peptide_id = gsub('\\..+$', '', accession_number))

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

genes_org <- unique(org_data_with_pep_id$ensembl_peptide_id)
g_list_org <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id", "ensembl_gene_id", "hgnc_symbol"),values=genes_org,mart=mart)
g_list_org <- g_list_org %>% mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl_peptide_id, hgnc_symbol))
org_data_with_gene_id <- merge(x = g_list_org, y = org_data_with_pep_id, by.x="ensembl_peptide_id", by.y="ensembl_peptide_id")

hh_data_with_pep_id <- healthy_human %>% 
  mutate(ensembl_peptide_id = gsub('\\..+$', '', accession_number)) %>% 
  filter(!str_detect(directory, pattern = "Testis_HLAI"))
genes_hh <- unique(hh_data_with_pep_id$ensembl_peptide_id)
g_list_hh <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id", "ensembl_gene_id", "hgnc_symbol"),values=genes_hh,mart=mart)
g_list_hh <- g_list_hh %>% mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl_peptide_id, hgnc_symbol))
hh_data_with_gene_id <- merge(x = g_list_hh, y = hh_data_with_pep_id, by.x="ensembl_peptide_id", by.y="ensembl_peptide_id")

thy_data_with_pep_id <- thymus %>% 
  mutate(ensembl_peptide_id = gsub('\\..+$', '', accession_number))
genes_thy <- unique(thy_data_with_pep_id$ensembl_peptide_id)
g_list_thy <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id", "ensembl_gene_id", "hgnc_symbol"),values=genes_thy,mart=mart)
g_list_thy <- g_list_thy %>% mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl_peptide_id, hgnc_symbol))
thy_data_with_gene_id <- merge(x = g_list_thy, y = thy_data_with_pep_id, by.x="ensembl_peptide_id", by.y="ensembl_peptide_id")

tidy_org_data_with_gene_id <- tidy_canonical_ip_data(org_data_with_gene_id)

tidy_org_data_with_gene_id <- tidy_org_data_with_gene_id %>% 
  mutate(hh_found = ifelse(hgnc_symbol %in% unique(hh_data_with_gene_id$hgnc_symbol), TRUE, FALSE),
         thymus_found = ifelse(hgnc_symbol %in% unique(thy_data_with_gene_id$hgnc_symbol), TRUE, FALSE)) 

###For analysis with HPC
tidy_org_data_with_gene_id %>%
  filter(hh_found == FALSE & thymus_found == FALSE) %>%
  dplyr::select(sequence) %>%
  mutate(command = paste(
    "grep", 
    " -l ",
    sequence,
    ' filtered_ribo_* | awk \'{print ', 
    paste('"', sequence, ':"', sep = ""),
    " $0}\'",
    " >> output_count_true_cr_only",
    sep = ""
  )) %>%
  dplyr::select(command) %>%
  write_csv("PATH", col_names = FALSE)

ribo <- read_ribo_output("PATH")

ribo %>%
  arrange(desc(num_ribo_tissues))
  group_by(num_ribo_tissues) %>% 
  count()

final_df <- left_join(tidy_org_data_with_gene_id, ribo, by =c("sequence"= "peptide"))

