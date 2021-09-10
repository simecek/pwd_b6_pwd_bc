library(tidyverse)

# read some GigaMUGA data to 
cols = cols(
  .default = col_character(),
  Chromosome = col_character(),
  Position = col_character()
)
gm_map_data <- read_csv("raw_data/MegaMUGA_13DecBC1Sona.csv", col_types = cols) %>% 
  select(Name, Chromosome, Position)
names(gm_map_data)[1] <- "Marker"

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
gm_snp_map <- read_csv("raw_data/gm_uwisc_v1.csv", col_types = col_types) %>%
  select(marker, chr, bp_mm10, cM_g2f1, strand, unique, unmapped)

cm2bp_data = gm_map_data %>%
  inner_join(gm_snp_map, by=c("Marker" = "marker")) %>%
  filter(`unique`) %>%
  select(Marker, chr, cM_g2f1, bp_mm10) %>%
  arrange(chr, cM_g2f1)

cm2bp_onepos <- function(chrom, pos_cm) {
  chr_data <- filter(cm2bp_data, chr == as.character(chrom))
  first_less <- max(which(chr_data$cM_g2f1 <= pos_cm))
  first_more <- min(which(chr_data$cM_g2f1 >= pos_cm))
  approx(x = chr_data$cM_g2f1[c(first_less, first_more)],
         y = chr_data$bp_mm10[c(first_less, first_more)],
         xout = pos_cm)$y
}

cm2bp <- function(chrom, pos_cm) {
  stopifnot(length(chrom) == length(pos_cm))
  
  output <- rep(NA, length(chrom))
  
  for (c in unique(as.character(chrom))) {
    chr_data <- filter(cm2bp_data, chr == c & !is.na(bp_mm10))
    
    for (i in which(chrom == c)) {
    
      first_less <- max(which(chr_data$cM_g2f1 <= pos_cm[i]))
      first_more <- min(which(chr_data$cM_g2f1 >= pos_cm[i]))
      
      if (chr_data$cM_g2f1[first_less] == chr_data$cM_g2f1[first_more]) {
        output[i] <- chr_data$bp_mm10[first_less]/2 + chr_data$bp_mm10[first_more]/2
      } else {
        output[i] <- approx(x = chr_data$cM_g2f1[c(first_less, first_more)],
                          y = chr_data$bp_mm10[c(first_less, first_more)],
                          xout = pos_cm[i])$y
      }  
    }
  }
  output
}
