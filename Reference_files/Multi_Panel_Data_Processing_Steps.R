st_path <- "PATH TO STs"

##### Figure 2C #####

read_xlsx(st_path, sheet = 4, skip = 1) %>% 
  filter(species %in% c("nuORFdb_human", "Human")) %>% 
  select(species, org_line, seq_upper) %>% 
  unique() %>% 
  group_by(species, org_line) %>% 
  count() %>% 
  ungroup() %>% 
  group_by(org_line) %>% 
  mutate(total = sum(n),
         perc = n/total)

##### Figure 2D #####

read_xlsx(st_path, sheet = 7, skip = 1) %>%
  select(seq_upper, biotype, org_line) %>% 
  unique() %>% 
  group_by(biotype) %>% 
  count() %>% 
  arrange(desc(n))

##### Figure 2E #####

read_xlsx(st_path, sheet = 7, skip = 1) %>%
  filter(str_detect(assign.MSi_allele, "A0201")) %>% 
  select(org_line, seq_upper) %>% 
  unique() %>% 
  group_by(seq_upper) %>% 
  count() %>% 
  group_by(n) %>% 
  count()

read_xlsx(st_path, sheet = 7, skip = 1) %>%
  filter(str_detect(assign.MSi_allele, "B4402")) %>% 
  select(org_line, seq_upper) %>% 
  unique() %>% 
  group_by(seq_upper) %>% 
  count() %>% 
  group_by(n) %>% 
  count()

read_xlsx(st_path, sheet = 7, skip = 1) %>%
  filter(str_detect(assign.MSi_allele, "C0501")) %>% 
  select(org_line, seq_upper) %>% 
  unique() %>% 
  group_by(seq_upper) %>% 
  count() %>% 
  group_by(n) %>% 
  count()

##### Figure 2F #####

read_xlsx(st_path, sheet = 4, skip = 1) %>%
  filter(species %in% c("Human", "nuORFdb_human")) %>% 
  select(org_line, seq_upper, species) %>%
  unique() %>% 
  mutate(length = nchar(seq_upper)) %>% 
  filter(between(length, 8, 11)) %>% 
  group_by(org_line, species, length) %>% 
  count() %>% 
  ungroup() %>% 
  group_by(org_line, species) %>% 
  mutate(total = sum(n),
         perc_1 = n/total) %>% 
  group_by(length, species) %>% 
  mutate(perc = mean(perc_1),
            stdev = sd(perc_1)) %>% 
  select(species, length, perc, stdev) %>% 
  unique()

read_xlsx(st_path, sheet = 4, skip = 1) %>%
  filter(species %in% c("Human", "nuORFdb_human")) %>% 
  select(org_line, seq_upper, species) %>%
  unique() %>% 
  mutate(length = nchar(seq_upper)) %>% 
  filter(between(length, 8, 11)) %>% 
  group_by(org_line, species, length) %>% 
  count() %>% 
  ungroup() %>% 
  group_by(org_line, species) %>% 
  mutate(total = sum(n),
         perc = n/total)

##### Figure 2J #####

read_xlsx(st_path, sheet = 4, skip = 1) %>% 
  filter(!str_detect(filename, "2D")) %>% 
  mutate(seq_upper = str_to_upper(sequence)) %>% 
  filter(!str_detect(assign.MSi_allele, pattern = "\\,")) %>% 
  filter(species %in% c("Human", "nuORFdb_human")) %>% 
  select(seq_upper, assign.MSi_allele, species, org_line) %>% 
  filter(assign.MSi_allele != "unknown") %>% 
  filter(!is.na(assign.MSi_allele)) %>% 
  filter(assign.MSi_allele != "N/A") %>% 
  filter(!is.na(species)) %>% 
  group_by(assign.MSi_allele, species) %>% 
  count() %>% 
  group_by(species) %>% 
  mutate(total = sum(n),
         perc = n/total)

##### Figure 3B #####
  
read_xlsx(st_path, sheet = 7, skip = 1) %>%
    select(seq_upper, hh_ip_detection, ribo_detection, thymus_ip_detection) %>% 
    unique() %>% 
    count()

read_xlsx(st_path, sheet = 7, skip = 1) %>%
    select(seq_upper, hh_ip_detection, ribo_detection, thymus_ip_detection) %>% 
    unique() %>% 
    filter(hh_ip_detection == FALSE) %>% 
    count()
  
read_xlsx(st_path, sheet = 7, skip = 1) %>%
    select(seq_upper, hh_ip_detection, ribo_detection, thymus_ip_detection) %>% 
    unique() %>% 
    filter(hh_ip_detection == FALSE & ribo_detection == FALSE) %>% 
    count()
  
read_xlsx(st_path, sheet = 7, skip = 1) %>%
    select(seq_upper, hh_ip_detection, ribo_detection, thymus_ip_detection) %>% 
    unique() %>% 
    filter(hh_ip_detection == FALSE & ribo_detection == FALSE & thymus_ip_detection == FALSE) %>% 
    count()

##### Figure 3D #####

read_xlsx(st_path, sheet = 7, skip = 1) %>%
  filter(cancer_restricted == TRUE) %>% 
  select(seq_upper, biotype) %>% 
  unique() %>% 
  group_by(biotype) %>% 
  count() %>% 
  arrange(desc(n))

##### Figure 3E #####

a <- read_xlsx(st_path, sheet = 7, skip = 1) %>% 
  filter(cancer_restricted == TRUE) %>% 
  select(seq_upper, biotype) %>% 
  unique() %>% 
  group_by(biotype) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  mutate(origin = "CR")

b <- read_xlsx(st_path, sheet = 7, skip = 1) %>% 
  select(seq_upper, biotype) %>% 
  unique() %>% 
  group_by(biotype) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  mutate(origin = "total")

rbind(a,b) %>% 
  pivot_wider(names_from = "origin", values_from = "n") %>% 
  mutate(perc = CR/total) %>% 
  arrange(desc(perc))

##### Figure 3F #####

cr_peps <- read_xlsx(st_path, sheet = 7, skip = 1) %>% 
  filter(cancer_restricted == TRUE) %>%  
  pull(seq_upper) %>% unique()

read_xlsx(st_path, sheet = 7, skip = 1) %>%
  select(seq_upper, assign.MSi_allele, org_line) %>% 
  unique() %>% 
  filter(seq_upper %in% cr_peps) %>% 
  filter(str_detect(assign.MSi_allele, "A0201")) %>% 
  group_by(seq_upper) %>% 
  count() %>% 
  group_by(n) %>% 
  count() %>% 
  select(recurrence = n,
         num_peps = nn)

##### Figure 3G #####

cr_peps <- read_xlsx(st_path, sheet = 7, skip = 1) %>% 
  filter(cancer_restricted == "TRUE") %>% 
  pull(seq_upper) %>% unique()

read_xlsx(st_path, sheet = 7, skip = 1) %>% 
  select(seq_upper, org_line) %>% 
  unique() %>% 
  filter(seq_upper %in% cr_peps) %>% 
  group_by(org_line) %>% 
  count() %>% 
  arrange(desc(n))

##### Figure S5E #####

read_xlsx(st_path, sheet = 7, skip = 1) %>%
  select(seq_upper, num_org_lines, biotype) %>% 
  filter(num_org_lines > 1) %>%
  select(seq_upper, num_org_lines, biotype) %>% 
  unique() %>% 
  group_by(biotype) %>% 
  count() %>% 
  arrange(desc(n))

##### Figure S6A #####

read_xlsx(st_path, sheet = 8, skip = 1) %>%
  mutate(seq_upper = str_to_upper(sequence),
         length = nchar(seq_upper)) %>%
  filter(!str_detect(accession_numbers, "\\|")) %>% 
  filter(between(length, left = 8, right = 11)) %>% 
  select(sample, seq_upper) %>% 
  unique() %>% 
  group_by(sample) %>% 
  count()

##### Figure S6B #####

cptac_nchlap <- read_xlsx(st_path, sheet = 8, skip = 1) %>%
  filter(species == "nuORFdb_human") %>% 
  mutate(seq_upper = str_to_upper(sequence),
         length = nchar(seq_upper)) %>% 
  transmute(seq_upper = str_to_upper(sequence), 
            sample,
            length = nchar(seq_upper),
            species, 
            accession_numbers) %>% 
  filter(between(length, left = 8, right = 11)) %>% 
  unique() %>% 
  group_by(sample) %>% 
  count() %>% 
  mutate(class = "nchlap")

cptac_canonical <- read_xlsx(st_path, sheet = 5, skip = 1) %>%
  mutate(seq_upper = str_to_upper(sequence),
         length = nchar(seq_upper)) %>% 
  filter(species == "Human") %>% 
  transmute(seq_upper = str_to_upper(sequence), 
            sample = Sample,
            length = nchar(seq_upper),
            species, accession_numbers) %>% 
  filter(between(length, left = 8, right = 11)) %>% 
  unique() %>% 
  group_by(sample) %>% 
  count() %>% 
  mutate(class = "canonical")

rbind(cptac_nchlap, cptac_canonical) %>% 
  pivot_wider(names_from = "class", values_from = "n") %>% 
  mutate(perc_nchlap = nchlap/(nchlap+canonical),
         perc_canonical = canonical/(nchlap+canonical)) %>% 
  pivot_longer(cols = c(4,5), names_to = "ident", values_to = "perc") %>% 
  group_by(sample) %>% 
  mutate(max = max(perc)) %>% 
  mutate(ident = factor(ident, levels = c("perc_nchlap", "perc_canonical")))


##### Figure S6C #####

nc_length <- read_xlsx(st_path, sheet = 8, skip = 1) %>% 
  transmute(seq_upper = str_to_upper(sequence),
            length = nchar(seq_upper),
            species, 
            accession_numbers) %>%
  filter(!str_detect(accession_numbers, "\\|")) %>% 
  unique() %>% 
  filter(between(length, left = 8, right = 11)) %>% 
  group_by(length) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(total = sum(n),
         perc = n/total,
         ident = "nchlap")
  
can_length <- read_xlsx(st_path, sheet = 5, skip = 1) %>%
  transmute(seq_upper = str_to_upper(sequence),
            length = nchar(seq_upper),
            species, 
            accession_numbers) %>%
  filter(!str_detect(accession_numbers, "\\|")) %>% 
  unique() %>% 
  filter(between(length, left = 8, right = 11)) %>% 
  group_by(length) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(total = sum(n),
         perc = n/total,
         ident = "canonical")

rbind(nc_length, can_length)

##### Figure S6D #####

read_xlsx(st_path, sheet = 8, skip = 1) %>% 
  filter(!str_detect(accession_numbers, "\\|")) %>% 
  select(seq_upper, biotype) %>% 
  unique() %>% 
  group_by(biotype) %>% 
  count() %>% 
  arrange(desc(n))


  