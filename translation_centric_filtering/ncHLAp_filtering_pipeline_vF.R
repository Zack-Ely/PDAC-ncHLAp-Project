library(tidyverse)
library(rebus)

##### Step 0: Combine TISH outputs #####
files <- list.files(path = ".", pattern = "*filtered.txt")

combined_tish_output <- data.frame()

for (file_path in files){
  tissue_name <- str_extract(file_path, "Brain|Fat|Fibroblast|HA_EC|HCAEC|Hepatocyte|Kidney|VSMC")
  tmp <- read_csv(file = file_path, col_names = FALSE, show_col_types = FALSE) %>%
    mutate(tissue = tissue_name) %>% 
    select(AASeq = X1, tissue)
  
  combined_tish_output <- rbind(combined_tish_output, tmp)
}

##### Step 1: Define paths #####
path_to_carr_output <- "./final_df_withalleles_LOCKED.csv"
path_to_hh_ip <- "./combined_hh_ip.csv"
path_to_thymus_ip <- "./thymus.txt"
output_file_directory <- "."
if (str_sub(output_file_directory, -1, -1) != "/"){
  output_file_directory <- paste(output_file_directory, "/", sep = "")
}

##### Step 2: Function to read in PDAC IP data #####
read_final_peptide_data <- function(path){
  library(tidyverse, quietly = TRUE)
  library(rebus, quietly = TRUE)
  library(stringr, quietly = TRUE)
  dataset <- suppressWarnings(read_delim(path, show_col_types = FALSE))
  
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

##### Step 3: Function to measure overlap given two sets of coordinates #####
perc_overlap = function(x.start, x.end, y.start, y.end) {
  if (x.start == y.start & x.end == y.end) {
    return(100)
  }
  if (x.start > x.end) {
    tmp <- x.end
    x.end <- x.start
    x.start <- tmp
  }
  if (y.start > y.end) {
    tmp <- y.end
    y.end <- y.start
    y.start <- tmp
  }
  x.len = abs(x.end - x.start)
  # largest start
  max.start = max(c(x.start, y.start))
  min.end = min(c(x.end, y.end))
  overlap = min.end - max.start
  overlap = ifelse(overlap <= 0, 0, overlap)
  perc_overlap = overlap / x.len * 100
  return(perc_overlap)
}

##### Step 4: Function to make tidy version of PDAC IP DF #####
tidy_ncmap_ip_data <- function(ip_data){
  tmp <- ip_data %>% 
    filter(ncmap == TRUE) %>% 
    mutate(sequence = str_to_upper(sequence)) %>% 
    select(sequence, org_line, culture, accession_number, ncmap, ncmap_class, 
           previous_aa, next_aa, length, bestAllele, binder = Binder) 
  
  tmp$num_org_lines <- "x"
  tmp$org_lines <- "x"
  tmp$alleles <- "x"
  
  for(i in seq(nrow(tmp))){
    seq <- tmp$sequence[i]
    tmp$num_org_lines[i] <- tmp %>% 
      filter(sequence == seq) %>% 
      select(org_line) %>% 
      unique() %>% 
      nrow()
    tmp$org_lines[i] <- tmp %>% 
      filter(sequence == seq) %>% 
      select(org_line) %>% 
      unique() %>% 
      as.list() %>% 
      .$org_line %>% 
      str_flatten(., ",")
    tmp$alleles[i] <- tmp %>%
      filter(sequence == seq) %>% 
      select(bestAllele) %>% 
      unique() %>% 
      as.list() %>% 
      .$bestAllele %>% 
      str_split(., ",") %>% 
      unlist() %>% 
      unique() %>% 
      str_flatten(., ",")
    
  }
  
  tmp <- tmp %>% 
    select(-c(bestAllele, binder, org_line)) %>% 
    unique()
  
  return(tmp)
}

##### Step 5: Read in PDAC IP #####
final_peptide_data <- tidy_ncmap_ip_data(read_final_peptide_data(path_to_carr_output))

##### Step 6: Read in HH IP #####
healthy_human <- read_csv(path_to_hh_ip)

##### Step 7: Read in Thymus IP #####
thymus <- read_delim(path_to_thymus_ip)

##### Step 8: Creating PDAC IP DF with coordinates #####

final_peptide_data <- final_peptide_data %>%
  filter(
    ncmap == TRUE,
    culture == "ORG"
  ) %>%
  mutate(
    start_coordinate = str_extract(
      accession_number,
      pattern = "\\:" %R% one_or_more("[:digit:]") %R% "\\-"
    ),
    end_coordinate = str_extract(
      accession_number,
      pattern = "\\-" %R% one_or_more("[:digit:]") %R% "\\:"
    )
  ) %>%
  mutate(start_coordinate = as.numeric(str_sub(
    start_coordinate, start = 2, end = -2
  )),
  end_coordinate = as.numeric(str_sub(
    end_coordinate, start = 2, end = -2
  ))) %>%
  mutate(chromosome = str_extract(
    str_extract(
      accession_number,
      pattern = "[:PUNCT:]" %R% one_or_more("[:ALNUM:]") %R% "\\:"
    ),
    one_or_more("[:ALNUM:]")
  ))

##### Step 9: Creating HH IP DF with coordinates #####
hh_coord <- healthy_human %>%
  filter(species == "nuORFdb_human",!str_detect(accession_numbers, pattern = "\\|")) %>%
  filter(!str_detect(directory, pattern = "Testis_HLAI")) %>%
  mutate(
    start_coordinate = str_extract(
      accession_numbers,
      pattern = "\\:" %R% one_or_more("[:digit:]") %R% "\\-"
    ),
    end_coordinate = str_extract(
      accession_numbers,
      pattern = "\\-" %R% one_or_more("[:digit:]") %R% "\\:"
    )
  ) %>%
  mutate(start_coordinate = as.numeric(str_sub(
    start_coordinate, start = 2, end = -2
  )),
  end_coordinate = as.numeric(str_sub(
    end_coordinate, start = 2, end = -2
  ))) %>%
  mutate(chromosome = str_extract(
    str_extract(
      accession_numbers,
      pattern = "[:PUNCT:]" %R% one_or_more("[:ALNUM:]") %R% "\\:"
    ),
    one_or_more("[:ALNUM:]")
  ))

##### Step 10: Creating Thymus IP DF with coordinates #####
thymus_pep_coord <- thymus %>%
  filter(species == "nuORFdb_human") %>% 
  select(directory, sequence, accession_number) %>% 
  mutate(start_coordinate = str_extract(accession_number, 
                                        pattern = "\\:" %R% one_or_more("[:digit:]") %R% "\\-"),
         end_coordinate = str_extract(accession_number, 
                                      pattern = "\\-" %R% one_or_more("[:digit:]") %R% "\\:")) %>% 
  mutate(start_coordinate = as.numeric(str_sub(start_coordinate, start = 2, end = -2)),
         end_coordinate = as.numeric(str_sub(end_coordinate, start = 2, end = -2))) %>% 
  mutate(chromosome = str_extract(str_extract(accession_number, pattern = "[:PUNCT:]" %R% one_or_more("[:ALNUM:]") %R% "\\:"), one_or_more("[:ALNUM:]")))


##### Step 11: Measuring ORF overlap with PDAC IP and HH IP #####
list_of_overlap <- c()
list_of_peps <- c()
overlap_list <- c()
list_of_indeces <- c()
for (i in seq(nrow(final_peptide_data))) {
  progress_bar = txtProgressBar(min=1, max=nrow(final_peptide_data), style = 3, char="=")
  overlap_list <- c()
  for (j in seq(nrow(hh_coord))) {
    if (final_peptide_data$chromosome[i] == hh_coord$chromosome[j]) {
      overlap_list <-
        c(
          overlap_list,
          perc_overlap(
            final_peptide_data$start_coordinate[i],
            final_peptide_data$end_coordinate[i],
            hh_coord$start_coordinate[j],
            hh_coord$end_coordinate[j]
          )
        )
    } else {
      overlap_list <- c(overlap_list, 0)
    }
  }
  list_of_overlap <-
    c(list_of_overlap, max(overlap_list, na.rm = TRUE))
  list_of_peps <-
    c(list_of_peps, final_peptide_data$sequence[i])
  list_of_indeces <- 
    c(list_of_indeces, list(which(overlap_list > 0)))
  setTxtProgressBar(progress_bar, value = i)
}


##### Step 12: Determining if peptide is found (exact match) in HH IP #####

hh_no_testis <- healthy_human %>% filter(!str_detect(directory, pattern = "Testis_HLAI"))

df_with_overlap <-
  data.frame(sequence = list_of_peps, overlap = list_of_overlap)

df_with_overlap$found <-
  df_with_overlap$sequence %in% str_to_upper(healthy_human$sequence)

##### Step 13: Determining number of HH IP tissues ORF is detected in #####
hh_coord$patient <- str_extract(hh_coord$directory, pattern = START %R% one_or_more("[:ALNUM:]"))

num_hh_patients <- c()
for (i in seq(nrow(df_with_overlap))){
  num_hh_patients[i] <- length(unique(hh_coord[list_of_indeces[[i]],]$patient))
}

df_with_overlap$num_hh_patients <- num_hh_patients

df_with_overlap <- unique(df_with_overlap)

final_peptide_data <- final_peptide_data %>% 
  left_join(., df_with_overlap, by = "sequence") 

##### Step 14: Measuring ORF overlap with PDAC IP and Thymus IP #####
list_of_overlap_thymus <- c()
list_of_peps_thymus <- c()
overlap_list_thymus <- c()
list_of_indeces_thymus <- c()
for (i in seq(nrow(final_peptide_data))) {
  progress_bar = txtProgressBar(min=1, max=nrow(final_peptide_data), style = 3, char="=")
  overlap_list_thymus <- c()
  for (j in seq(nrow(thymus_pep_coord))) {
    if (final_peptide_data$chromosome[i] == thymus_pep_coord$chromosome[j]) {
      overlap_list_thymus <-
        c(
          overlap_list_thymus,
          perc_overlap(
            final_peptide_data$start_coordinate[i],
            final_peptide_data$end_coordinate[i],
            thymus_pep_coord$start_coordinate[j],
            thymus_pep_coord$end_coordinate[j]
          )
        )
    } else {
      overlap_list_thymus <- c(overlap_list_thymus, 0)
    }
  }
  list_of_overlap_thymus <-
    c(list_of_overlap_thymus, max(overlap_list_thymus))
  list_of_peps_thymus <-
    c(list_of_peps_thymus, final_peptide_data$sequence[i])
  list_of_indeces_thymus <- 
    c(list_of_indeces_thymus, list(which(overlap_list_thymus > 0)))
  setTxtProgressBar(progress_bar, value = i)
}

##### Step 15: Determining number of thymus samples ORF is detected in #####
thymus_overlap <- data.frame(sequence = list_of_peps_thymus,
                             thymic_overlap = list_of_overlap_thymus)

thymus_pep_coord$patient <- str_remove(thymus_pep_coord$directory, pattern = START %R% "THY" %R% "\\_")

num_thy_patients <- c()
for (i in seq(nrow(thymus_overlap))){
  num_thy_patients[i] <- length(unique(thymus_pep_coord[list_of_indeces_thymus[[i]],]$patient))
}

thymus_overlap$num_thy_patients <- num_thy_patients

thymus_overlap <- unique(thymus_overlap)

##### Step 16: Creating final DF with data from HH and IP overlaps #####
final_peptide_data <- left_join(final_peptide_data, thymus_overlap, by = "sequence")

##### Step 17: Read in ribo-seq, query for ORF translation in HH RiboSeq for all peptides #####
num_nchlap <- nrow(final_peptide_data)

empty_matrix <- data.frame(sequence = final_peptide_data$sequence,
                           Brain = rep(0, num_nchlap),
                           Fat = rep(0, num_nchlap),
                           Fibroblast = rep(0, num_nchlap),
                           HA_EC = rep(0, num_nchlap),
                           HCAEC = rep(0, num_nchlap),
                           Hepatocyte = rep(0, num_nchlap),
                           Kidney = rep(0, num_nchlap),
                           VSMC = rep(0, num_nchlap))

for (i in seq(nrow(final_peptide_data))) {
  progress_bar = txtProgressBar(min=1, max=nrow(final_peptide_data), style = 3, char="=")
  count_output <- combined_tish_output %>% 
    filter(str_detect(AASeq, final_peptide_data$sequence[i])) %>% 
    group_by(tissue) %>% 
    count()
  if (!identical(integer(0), count_output$n)) {
    for (j in seq(nrow(count_output))){
      empty_matrix[i, count_output[j,1]$tissue] <- count_output[j,2]$n
    } 
  }
  setTxtProgressBar(progress_bar, value = i)
}

ribo_orf_counts <- empty_matrix

rm(empty_matrix)

##### Step 18: Determine number of tissues each peptide is detected in #####
numbers <- ribo_orf_counts %>% 
  pivot_longer(cols = -c(1)) %>% 
  filter(value != 0) %>% 
  group_by(sequence) %>% 
  count() %>% 
  select(sequence, num_ribo_tissues = n)


##### Step 19: Create final DF with all data sources #####
left_join(final_peptide_data, ribo_orf_counts, by = "sequence") %>% 
  mutate(total = Brain + Fat + Fibroblast + HA_EC + HCAEC + Hepatocyte + Kidney + VSMC) %>% 
  left_join(., numbers, by = "sequence") %>% 
  select(sequence, 
         num_org_lines, 
         hh_overlap = overlap, 
         found_hh = found,
         num_hh_patients, 
         num_ribo_tissues,
         total_ribo_orfs = total, 
         thymic_overlap,
         num_thy_patients,
         everything()) %>%
  select(-colnames(ribo_orf_counts[-1])) %>% 
  mutate(num_ribo_tissues = ifelse(is.na(num_ribo_tissues), 0, num_ribo_tissues)) %>% 
  mutate(cancer_restricted = ifelse(hh_overlap == 0 & found_hh == FALSE & thymic_overlap == 0 & total_ribo_orfs == 0, "CR", "nonCR")) %>% 
  write_csv(paste(output_file_directory, "v1_withkidney_nolongest_tcfp_output.csv", sep = ""))



