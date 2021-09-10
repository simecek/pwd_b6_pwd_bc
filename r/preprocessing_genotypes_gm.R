library(tidyverse)

cols = cols(
  .default = col_character(),
  Chromosome = col_character(),
  Position = col_character()
)
gm_data <- read_csv("raw_data/MegaMUGA_13DecBC1Sona.csv", col_types = cols)
  dim(gm_data)
names(gm_data)[1] <- "Marker"


# Adding map ----------------------------------------

col_types = cols(
  marker = col_character(),
  chr = col_character(),
  bp_mm10 = col_double(),
  cM_cox = col_double(),
  cM_g2f1 = col_double(),
  strand = col_character(),
  snp = col_character(),
  unique = col_logical(),
  multi = col_logical(),
  unmapped = col_logical(),
  n_blast_hits = col_double(),
  n_blast_chr = col_double(),
  probe = col_character()
)

snp_map <- read_csv("raw_data/gm_uwisc_v1.csv", col_types = col_types) %>%
  select(marker, chr, bp_mm10, cM_g2f1, strand, unique, unmapped)
gm_data <- filter(gm_data, Marker %in% snp_map$marker)

all(gm_data$Marker %in% snp_map$marker)
all(!duplicated(gm_data$Marker))


mapped_data = gm_data %>%
  filter(`PWD.07` == `PWD.12`) %>% # consistent
  filter(`PWD.07` != `B6`) %>%  # informative
  filter(`PWD.07` != 'N', `B6` != 'N') %>% # non-mising
  filter(`PWD.07` != 'H', `B6` != 'H', ((`PWDxB6` == 'H') | (Chromosome == 'X')&(`PWDxB6` == `PWD.07`))) %>% 
  left_join(snp_map, by=c("Marker" = "marker")) %>%
  filter(`chr` %in% c(as.character(1:19), "X") )

# Recoding ----------------------------------------------------------------

genotype_data = select(mapped_data, -unique, -unmapped, -strand, -Chromosome, -Position, -SNP)
male_cols = names(genotype_data)[grep("^[0-9]+$", names(genotype_data))]
autosomes = as.character(1:19)

for (c in male_cols) {
  new_geno = rep("N", nrow(genotype_data))
  
  p_idx = (mapped_data[[c]] == mapped_data[["PWD.07"]]) & (mapped_data$chr %in% autosomes)
  new_geno[p_idx] <- "P"
  h_idx = (mapped_data[[c]] == "H") & (mapped_data$chr %in% autosomes)
  new_geno[h_idx] <- "H"
  
  p_idx = (mapped_data[[c]] == mapped_data[["PWD.07"]]) & (mapped_data$chr == "X")
  new_geno[p_idx] <- "P"
  b_idx = (mapped_data[[c]] == mapped_data[["B6"]]) & (mapped_data$chr == "X")
  new_geno[b_idx] <- "B"
  
  genotype_data[c] = new_geno
}

head(genotype_data)
sort(apply(genotype_data[,male_cols] == "N", 2, mean))  # 2114 is clearly off and should be removed
genotype_data <- select(genotype_data, -`2114`)
male_cols <- setdiff(male_cols, "2114")
names(genotype_data)

# Final preprocessing -----------------------------------------------------

output_data <- select(genotype_data, Marker, chr, bp_mm10, cM_g2f1, all_of(male_cols))

# reorder
output_data[['chr2']] <- as.numeric(output_data$chr)  # warning expected
output_data$chr2[output_data$chr == "X"] = 20

output_data = output_data %>% arrange(chr2, bp_mm10)

# save

output_data = select(output_data, -chr2)

sample_cols = names(output_data)[grep("^[0-9]+$", names(output_data))]
pct_n = apply(output_data[,sample_cols] == "N", 1, mean)
hist(pct_n)
table(pct_n > 0.1)
output_data <- output_data[pct_n <= 0.1, ]
table(output_data$chr)

write_csv(output_data, "data/genotypes_gm.csv")
