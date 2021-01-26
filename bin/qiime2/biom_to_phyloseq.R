library(biomformat)
library(phyloseq)
library(ggplot2)
library(ape)
library(tidyr)


source("/SAN/Susanas_den/gitProj/HMHZ/bin/1_data_prep.R")

otu <- read.table("tmp/qiime2/phyloseq_pooled/otu_table.txt", header=TRUE)
tax <- read.table("tmp/qiime2/phyloseq_pooled/taxonomy.tsv", header=TRUE, sep="\t")
colnames(otu)
rownames(otu) <- otu[,1]
otu <- otu[,-1]
tax[1,1]
rownames(tax) <- tax[,1]
tax <- tax[,-1]
tax[1,]
tax <- tax %>%
    separate(Taxon, c("domain", "phyla", "class", "order", "family", "genus", "species"),"; ")
head(tax)
OTU <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))
META <- sample_data(sample.data)

sample_names(OTU)
sample_names(META)
taxa_names(OTU)
taxa_names(TAX)

sample_names(OTU) <- gsub("^.*?\\.", "", sample_names(OTU))
sample_names(OTU) <- gsub("\\.", "_", sample_names(OTU))

sample_names(OTU)
sample_names(META)

ps_d <- phyloseq(OTU, TAX, META)
ps_d

saveRDS(ps_d, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/PSqiime.RDS")

## again for the short reads dataset

otu <- read.table("tmp/qiime2/phyloseq_pooled200/otu_table.txt", header=TRUE)
tax <- read.table("tmp/qiime2/phyloseq_pooled200/taxonomy.tsv", header=TRUE, sep="\t")
colnames(otu)
rownames(otu) <- otu[,1]
otu <- otu[,-1]
tax[1,1]
rownames(tax) <- tax[,1]
tax <- tax[,-1]
tax[1,]
tax <- tax %>%
    separate(Taxon, c("domain", "phyla", "class", "order", "family", "genus", "species"),"; ")
head(tax)

OTU <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))
META <- sample_data(sample.data)

sample_names(OTU)

sample_names(META)

taxa_names(OTU)

taxa_names(TAX)

sample_names(OTU) <- gsub("^.*?\\.", "", sample_names(OTU))
sample_names(OTU) <- gsub("\\.", "_", sample_names(OTU))

sample_names(OTU)
sample_names(META)

ps_d <- phyloseq(OTU, TAX, META)
ps_d

saveRDS(ps_d, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/PSqiime200.RDS")
