library(tidyverse)
library(qtl)


# Joining MM a GM genotypes -----------------------------------------------

gm_data <- read_csv("data/genotypes_gm.csv", col_types=cols(chr = col_character())) %>%
  filter(chr %in% c(as.character(1:19), "X"))
dim(gm_data)

small_data <- select(gm_data, -bp_mm10)
rm(gm_data)

# Checking phenotypes -----------------------------------------------------

col_types = cols(
  ID = col_double(),
  `D17/334` = col_character(),
  BW = col_number(),
  TW = col_double(),
  SC = col_double(),
  MegaMuga = col_double()
)
pheno <- read_csv("raw_data/PWDx(B6xPWD)MUGA_phenotypes.csv", col_types = col_types)

# those are missing
missing_pheno_ids <- names(small_data)[-(1:3)][!(as.numeric(names(small_data)[-(1:3)]) %in% pheno$ID)]
missing_pheno_ids

# samples with missing phenotypes should be excluded
dim(small_data)
small_data <- small_data[,!(names(small_data) %in% missing_pheno_ids)]
dim(small_data)
names(small_data)[1] <- "id"
names(small_data)[2:3] <- ""
write_csv(small_data, "data/genotypes_qtl.csv")

pheno_sel <- pheno[match(names(small_data)[-(1:3)], pheno$ID),]
names(pheno_sel)[1] <- "id"
dim(pheno_sel)
pheno_sel$infertility_cat <- as.numeric((pheno_sel$TW < 80) & (pheno_sel$SC < 5))

write_csv(pheno_sel, "data/phenotypes.csv")

# Transforming phenotypes for QTL mapping ---------------------------------

# TW
rpheno <- t(pheno_sel %>%
              select(TW, id)) 
write.table(as.data.frame(rpheno), "data/pheno_TW.csv", col.names = FALSE, row.names=TRUE, sep=",")

#SC
rpheno <- pheno_sel  %>%
  select(SC, id) 
rpheno$SC <- log(1+rpheno$SC)
write.table(as.data.frame(t(rpheno)), "data/pheno_SC.csv", col.names = FALSE, row.names=TRUE, sep=",")

#Infertility (cat.)
rpheno <- pheno_sel  %>%
  select(infertility_cat, id) 
write.table(as.data.frame(t(rpheno)), "data/pheno_infertility_cat.csv", col.names = FALSE, row.names=TRUE, sep=",")


