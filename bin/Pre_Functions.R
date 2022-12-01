fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # abundance filtering at 0.005%
    keepTaxa = (x / sum(x) > 0.00005)
#    keepTaxa = (x / sum(x) > 0.0005)
    summary(keepTaxa)
    ps = phyloseq::prune_taxa(keepTaxa, ps)
# plus prevalnce filter at 1%
    KeepTaxap <- microbiome::prevalence(ps)>0.01
    ps <- phyloseq::prune_taxa(KeepTaxap, ps)
# subset samples based on total read count (100 reads)
#ps <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 100)
    ps <- phyloseq::prune_samples(sample_sums(ps)>100, ps)
    ps
}

library(ggplot2)
library(reshape)
#library(phyloseq, lib="/usr/local/lib/R/site-library/")
library(data.table)
#library(parallel)
library(microbiome)
#library("pheatmap")
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
library(vegan)
library(tidyr)
library(tidyverse)
#library(geosphere)
library(ggpubr)
## using the devel
devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")
library("ape")
library(magrittr)
library(Biostrings)
library(DECIPHER)
library(ShortRead)
library(RColorBrewer)
library(phyloseq)

#### preprocessing
PS <- readRDS(file = "tmp/interData/PhyloSeqCombi_HMHZ_All.Rds")
PS.l <- readRDS(file = "tmp/interData/PhyloSeqList_HMHZ_All.Rds")

## filtering MA by amplicon
fPS.l <- list()
for (i in 1:38) {
    try(fPS.l[[i]] <- fil(PS.l[[i]]), silent=TRUE)
}

fPS.l

fPS <- fPS.l[[1]]

for (i in 2:38){
    fPS <- try(merge_phyloseq(fPS,fPS.l[[i]]))
#    print(fPS)
    }

# filter pooled for comparison
fPS.pool <- fil(PS)

# How many Eimeria asv's do we have?
subset_taxa(PS, Genus%in%"g__Eimeria") # 225 before filtering
#subset_taxa(PS, genus%in%"Eimeria") # 144 before filtering

subset_taxa(fPS, Genus%in%"g__Eimeria") # 225 before filtering

#get_taxa_unique(fPS.l[[2]], "family")

get_taxa_unique(fPS, "Family")

get_taxa_unique(fPS, "Phylum")

#how many primers amplify Apicomplexa and which families?
for (i in 1:44) {
#    print(names(all.PS.l)[[i]])
    try(p <- subset_taxa(fPS.l[[i]],Phylum=="p__Apicomplexa"), silent=TRUE)
#    try(get_taxa_unique(p, "family"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "Family")
        print(paste(i, "- ", names(fPS.l[i]), ": ", length(a), sep=""))
        print(a)
}
    rm(p)
}

##### and how many amplicons have eimeria?
for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(PS.l[[i]],Genus=="g__Eimeria"), silent=TRUE)
#    try(get_taxa_unique(p, "genus"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "Genus")
        print(paste(i, "- ", names(PS.l[i]), ": ", nrow(p@tax_table), sep=""))
        print(a)
}
    rm(p)
}


####
#fPS.l <- list()
#for (i in 1:48) {
#    try(fPS.l[[i]] <- fil(PS.l[[i]]), silent=TRUE)
#}

for (i in 1:48) {
    try(names(fPS.l[i]) <- names(PS.l[i]), silent=TRUE)
}


### and what happens when we filter?
nmEim <- list()
names18S <- list()

for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(fPS.l[[i]],Genus=="g__Eimeria"), silent=TRUE)
#    try(get_taxa_unique(p, "genus"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "Genus")
        print(paste(i, "- ", names(PS.l[i]), ": ", nrow(p@tax_table), sep=""))
#        print(a)
        nmEim[i] <- names(PS.l[i])
        names18S[[i]] <- paste(rep(names(PS.l[i]), nrow(p@tax_table)), "ASV", seq(1, nrow(p@tax_table), 1), sep="_") 
}
    rm(p)
}

nmEim <- unlist(nmEim)
names18S <- unlist(names18S)

names18S

## 19 amplicons have 59 EImeira ASV's
# save ASV's from 18S amplicons

# first we need to know which genes are we targeting
primerL <- read.csv("/SAN/Susanas_den/HMHZ/data/primerInputUnique.csv")
#quick fix here
primerL$Primer_name[122] <- "27M_F_98_F.Klin0341_CR_18_R"

#targets
primerL$Gen[which(primerL$Primer_name%in%nmEim)]
primerL$Primer_name[which(primerL$Primer_name%in%nmEim)]

# these amplicons are not in the list but they target 18S, 18S and 28S
nmEim[which(!nmEim%in%primerL$Primer_name)]

# amplicons targeting 28S: "D3A_5Mod_46_F.D3B_5Mod_46_R" and "NLF184cw_74_F.NL818cw_74_R"
# we need to ignore the 28S Eimeria ASV's
#Eim.28S <- merge_phyloseq(subset_taxa(fPS.l[[12]], genus%in%"Eimeria"), subset_taxa(fPS.l[[28]], genus%in%"Eimeria"))
#ASV.28S <- rownames(Eim.28S@tax_table)

fPS.TSS <- transform_sample_counts(fPS, function(x) x / sum(x))
T.PS <- fPS.TSS
T.PS@otu_table <- T.PS@otu_table*T.PS@sam_data$Concentration

Eim <- subset_taxa(fPS, Genus%in%"g__Eimeria")
Eim.TSS <- subset_taxa(fPS.TSS, Genus%in%"g__Eimeria")
Eim.T <- subset_taxa(T.PS, Genus%in%"g__Eimeria")
#Eim <- prune_taxa(!rownames(Eim@tax_table)%in%ASV.28S, Eim)

Eim <- prune_taxa(taxa_sums(Eim)>0, Eim)
Eim.TSS <- prune_taxa(taxa_sums(Eim.TSS)>0, Eim.TSS)
Eim.T <- prune_taxa(taxa_sums(Eim.T)>0, Eim.T)



# separating 18S and 28S
Eim.ASV <- DNAStringSet(taxa_names(Eim))

names18S

names(Eim.ASV) <- names18S

Eim.ASV

Seq18 <- Eim.ASV[which(!grepl("D3A", names(Eim.ASV)))]
Seq28 <- Eim.ASV[which(grepl("D3A", names(Eim.ASV)))]

writeFasta(Eim.ASV, "tmp/EimeriaASV.fasta")
writeFasta(Seq18, "tmp/Eimeria18S.fasta")
writeFasta(Seq28, "tmp/Eimeria28S.fasta")



############## old preproc

##Eliminate Unassigned to superkingdom level
#PS <- subset_taxa(PS, !is.na(superkingdom) & !superkingdom %in% c("", "uncharacterized"))

# subset samples based on total read count (1000 reads)
#median(phyloseq::sample_sums(PS))

#hist(phyloseq::sample_sums(PS))

#PS <- phyloseq::subset_samples(PS, phyloseq::sample_sums(PS) > 1000)

# prevalence filtering at 5%
#pPPS=phyloseq_filter_prevalence(pPS, prev.trh=0.05)

###A lot of Mus :(
## Host read numbers
#sum(otu_table(subset_taxa(PS, genus%in%"Mus")))/sum(otu_table(PS))
###Eliminate reads assigned as "Mus"
#PS <- subset_taxa(PS, !genus %in% "Mus") ##Eliminate reads :S

# Eliminate samples with no reads
#PS <- prune_samples(sample_sums(PS)>0, PS)

# abundance filtering to 0.01%? Or keep prevalence filtering?
#x = taxa_sums(PS)
#keepTaxa = (x / sum(x) > 0.0001)
#summary(keepTaxa)
#pPS = prune_taxa(keepTaxa, PS)

# Eliminate samples with no reads
#pPS <- prune_samples(sample_sums(pPS)>0, pPS)
