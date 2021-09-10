library(tidyverse)
library(qtl)


# TW ----------------------------------------------------------------------

bcp <- read.cross(format = "csvsr", dir="data/", "genotypes_qtl.csv", "pheno_TW.csv", genotypes=c("P","H","B","D","C"), na="N")
bcp <- calc.genoprob(bcp, step=1, error.prob=0.01)
bcp.scanone <- scanone(bcp, method="hk")
bcp.perm <- scanone(bcp, method="hk", n.perm=1000, n.cluster = 7)

save(bcp, bcp.scanone, bcp.perm, file = "data/scanone_TW.rdata")

summary(bcp.perm)
plot(bcp.scanone, ylim = c(0, max(bcp.scanone$lod, summary(bcp.perm)[[1]])), main="TW")
abline(h=summary(bcp.perm)[[1]], col="red", lty=2)
abline(h=summary(bcp.perm)[[2]], col="blue", lty=2)

# SC  ----------------------------------------------------------------------

bcp <- read.cross(format = "csvsr", dir="data/", "genotypes_qtl.csv", "pheno_SC.csv", genotypes=c("P","H","B","D","C"), na="N")
bcp <- calc.genoprob(bcp, step=1, error.prob=0.01)
bcp.scanone <- scanone(bcp, method="hk")
bcp.perm <- scanone(bcp, method="hk", n.perm=1000, n.cluster = 7)

save(bcp, bcp.scanone, bcp.perm, file = "data/scanone_SC.rdata")

summary(bcp.perm)
plot(bcp.scanone, ylim = c(0, max(bcp.scanone$lod, summary(bcp.perm)[[1]])), main="logSC")
abline(h=summary(bcp.perm)[[1]], col="red", lty=2)
abline(h=summary(bcp.perm)[[2]], col="blue", lty=2)


# INFERTILITY (cat.) ------------------------------------------------------

bcp <- read.cross(format = "csvsr", dir="data/", "genotypes_qtl.csv", "pheno_infertility_cat.csv", genotypes=c("P","H","B","D","C"), na="N")
bcp <- calc.genoprob(bcp, step=1, error.prob=0.01)
bcp.scanone <- scanone(bcp, method="hk", model="binary")
bcp.perm <- scanone(bcp, method="hk", model="binary", n.perm=500, n.cluster = 7)

save(bcp, bcp.scanone, bcp.perm, file = "data/scanone_infertility_cat.rdata")

summary(bcp.perm)
plot(bcp.scanone, ylim = c(0, max(bcp.scanone$lod, summary(bcp.perm)[[1]])), main="Infertility (cat.)")
abline(h=summary(bcp.perm)[[1]], col="red", lty=2)
abline(h=summary(bcp.perm)[[2]], col="blue", lty=2)


