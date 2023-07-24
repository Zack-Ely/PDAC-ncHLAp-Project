library(tidyverse)
library(stringr)
library(rebus)
library(readxl)

### This reads in IP into a single DF for the FINAL IP DF
### Creates columns for cell culture, HLA class, whether it is a nuORF, and what nuORF class, etc
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

### Reads in healthy human data from dropbox.
read_healthy_human <- function(file){
  healthy_human <- NULL
  for (i in seq(33)) {
    x <- dir(file)
    if (str_detect(x[i], pattern = START %R% "DN")) {
      healthy_human <-
        rbind(healthy_human, read_delim(paste(file, x[i], sep = "")))
    }
  }
  return (healthy_human)
}

### Function to measure overlap
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

### Function to make a tidy version of the IP data
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

### Org IP data
final_peptide_data <- tidy_ncmap_ip_data(read_final_peptide_data("PATH"))

### Healthy Human
healthy_human <- read_healthy_human("PATH")

### Thymus data
thymus <- read_delim("PATH")

### Creating data frame with coordinates, chromosome for each nuorf peptide
### Only uses peptides which are uniquely mapped to one transcript

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

### Creating data frame with coordinates, chromosome for each nuorf peptide
### Only uses peptides which are uniquely mapped to one transcript
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

thymus_pep_coord <- thymus %>%
  filter(species == "nuORFdb_human") %>% 
  dplyr::select(directory, sequence, accession_number) %>% 
  mutate(start_coordinate = str_extract(accession_number, 
                                        pattern = "\\:" %R% one_or_more("[:digit:]") %R% "\\-"),
         end_coordinate = str_extract(accession_number, 
                                      pattern = "\\-" %R% one_or_more("[:digit:]") %R% "\\:")) %>% 
  mutate(start_coordinate = as.numeric(str_sub(start_coordinate, start = 2, end = -2)),
         end_coordinate = as.numeric(str_sub(end_coordinate, start = 2, end = -2))) %>% 
  mutate(chromosome = str_extract(str_extract(accession_number, pattern = "[:PUNCT:]" %R% one_or_more("[:ALNUM:]") %R% "\\:"), one_or_more("[:ALNUM:]")))



### Measures overlap for a peptide against ALL nuorf peptide-encoding transcripts from
### normal tissue immunopeptidomics, and then takes the maximum overlap, and appends this to a list
list_of_overlap <- c()
list_of_peps <- c()
overlap_list <- c()
list_of_indeces <- c()
for (i in seq(nrow(final_peptide_data))) {
  print(i)
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
}

### Output dataframe from overlap loop, also determines here if each peptide
### is found (identically) in normal tissue IP
df_with_overlap <-
  data.frame(sequence = list_of_peps, overlap = list_of_overlap)
df_with_overlap$found <-
  df_with_overlap$sequence %in% str_to_upper(healthy_human$sequence)

### Identify the number of patients that each ORF is detected in
hh_coord$patient <- str_extract(hh_coord$directory, pattern = START %R% one_or_more("[:ALNUM:]"))

num_hh_patients <- c()
for (i in seq(nrow(df_with_overlap))){
  num_hh_patients[i] <- length(unique(hh_coord[list_of_indeces[[i]],]$patient))
}

df_with_overlap$num_hh_patients <- num_hh_patients

df_with_overlap <- unique(df_with_overlap)

final_peptide_data <- final_peptide_data %>% 
  left_join(., df_with_overlap, by = "sequence") 

### Below is the overlap analysis for thymus. Same structure as above
list_of_overlap_thymus <- c()
list_of_peps_thymus <- c()
overlap_list_thymus <- c()
list_of_indeces_thymus <- c()
for (i in seq(nrow(final_peptide_data))) {
  print(i)
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
}

thymus_overlap <- data.frame(sequence = list_of_peps_thymus,
                             thymic_overlap = list_of_overlap_thymus)

thymus_pep_coord$patient <- str_remove(thymus_pep_coord$directory, pattern = START %R% "THY" %R% "\\_")

num_thy_patients <- c()
for (i in seq(nrow(thymus_overlap))){
  num_thy_patients[i] <- length(unique(thymus_pep_coord[list_of_indeces_thymus[[i]],]$patient))
}

thymus_overlap$num_thy_patients <- num_thy_patients

thymus_overlap <- unique(thymus_overlap)

### Creates final data frame with data from healthy tissue and thymus

final_peptide_data <- left_join(final_peptide_data, thymus_overlap, by = "sequence")


### Take CSV, remove header, add info to make batch script to HPC
final_peptide_data %>%
  dplyr::select(sequence) %>%
  mutate(command = paste(
    "grep", 
    " -c ",
    sequence,
    ' filtered_ribo_* | awk \'{print ', 
    paste('"', sequence, ':"', sep = ""),
    " $0}\'",
    " >> output_count",
    sep = ""
  )) %>%
  select(command) %>% 
  write_csv("PATH", col_names = FALSE)

### Parse the HPC output
count <-
  read_delim(
    "PAHT",
    col_names = FALSE
  ) %>% 
  mutate(X2 = str_remove(X2, pattern = "filtered_ribo_"),
         X2 = str_remove(X2, pattern = ".txt")) %>% 
  pivot_wider(names_from = X2, values_from = X3) %>% 
  select(sequence = X1, everything())

numbers <- read_delim(
  "PATH",
  col_names = FALSE) %>% 
  filter(X3 > 0) %>% 
  mutate(X2 = str_remove(X2, pattern = "filtered_ribo_"),
         X2 = str_remove(X2, pattern = ".txt")) %>% 
  select(-X3) %>% 
  group_by(X1) %>% 
  count() %>% 
  select(sequence = X1, num_ribo_tissues = n)
  

### This is the amalgamation of all data sources. Creates DF with peptide, overlap with normal IP,
### whether the peptide is found in normal IP, whether the peptide is found in riboseq, and
### the overlap with thymus ORFs.

final_df <- left_join(final_peptide_data, count, by = "sequence") %>% 
  mutate(total = brain + fat + fibroblast + ha_ec + hcaec + hepatocyte + vsmc) %>% 
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
  mutate(num_ribo_tissues = ifelse(is.na(num_ribo_tissues), 0, num_ribo_tissues))

  