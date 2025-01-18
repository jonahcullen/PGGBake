#!/usr/bin/env Rscript

# modified from https://github.com/AnimalGenomicsETH/bovine-graphs

### Given mash distance construct phylogenetic tree from it

#loading library 
library("tidyverse")
library("ape")
library("ragg")

ref <- snakemake@params[["ref"]]
outmat <- snakemake@output[[1]]
outviz  <- snakemake@output[[2]]
disfile  <- snakemake@input[[1]]

datdis  <- read.table(disfile,header=FALSE, stringsAsFactors =FALSE)

#rename the header
colnames(datdis)  <- c("anim1","anim2","distr","comp4","comp5")

# give correct assembly name
pattern <- "(?<=/)[^/]+(?=\\.clean\\.fa)"
datdis$anim1c  <- str_extract(pattern = pattern, string = datdis$anim1)
datdis$anim2c  <- str_extract(pattern = pattern, string = datdis$anim2)

datwide <- datdis %>%
    dplyr::select(anim1c, anim2c, distr) %>%
    pivot_wider(names_from = anim2c, values_from = distr)

datmat <- as.matrix(datwide %>% select(-anim1c))
rownames(datmat) <- datwide$anim1c
print(rownames(datmat))
write.table(datmat, file = outmat, sep = "\t", quote = FALSE)

# select outgroup from the distance matrix
outgroup <- which.max(datmat[ref,]) %>% names()
print(paste("Outgroup:", outgroup))

# apply neighbor joining 
tr  <- nj(datmat)

# visualize the tree
agg_tiff(outviz, width = 600, height = 500, units="px")

plot.phylo(root(tr, outgroup = outgroup), cex=2, edge.width=1)
axisPhylo(backward=FALSE, cex.axis=2)

dev.off()

