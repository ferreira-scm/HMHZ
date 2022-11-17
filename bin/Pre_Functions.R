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
# subset samples based on total read count (500 reads)
#ps <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 500)
    ps <- phyloseq::prune_samples(sample_sums(ps)>500, ps)
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

fPS <- fPS.l[[1]]

for (i in 2:36){
    fPS <- try(merge_phyloseq(fPS,fPS.l[[i]]))
#    print(fPS)
    }

# filter pooled for comparison
fPS.pool <- fil(PS)

# How many Eimeria asv's do we have?
subset_taxa(PS, genus%in%"Eimeria") # 144 before filtering

#how many primers amplify Apicomplexa and which families?
for (i in 1:44) {
#    print(names(all.PS.l)[[i]])
    try(p <- subset_taxa(fPS.l[[i]],phylum=="Apicomplexa"), silent=TRUE)
#    try(get_taxa_unique(p, "family"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "family")
        print(paste(i, "- ", names(fPS.l[i]), ": ", length(a), sep=""))
        print(a)
}
    rm(p)
}

##### and how many amplicons have eimeria?
for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(PS.l[[i]],genus=="Eimeria"), silent=TRUE)
#    try(get_taxa_unique(p, "genus"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "genus")
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
    try(p <- subset_taxa(fPS.l[[i]],genus=="Eimeria"), silent=TRUE)
#    try(get_taxa_unique(p, "genus"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "genus")
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
Eim <- subset_taxa(fPS, genus%in%"Eimeria")
Eim.TSS <- subset_taxa(fPS.TSS, genus%in%"Eimeria")
#Eim <- prune_taxa(!rownames(Eim@tax_table)%in%ASV.28S, Eim)


# let's do some co-infection analyses. First we need to do some preparations
amplicon <- data.frame(names18S)
amplicon$species <- "sp"

# sanity check
rownames(Eim@tax_table)==rownames(Eim.TSS@tax_table)
Eim.TSS@tax_table[,6] <- names18S
Eim@tax_table[,6] <- names18S
names18S

## import tree
library(ape)
library(ggtree)

# I need set 3 groups of reference sequences (ferrisi, vermiformis, falciformis)
fer <- phyloT$tip.label[grep("ferrisi",phyloT$tip.label)]
fal <- c(phyloT$tip.label[grep("falciformis",phyloT$tip.label)], "KU174464_Eimeria_apionodes", "AF311644_Eimeria_sevilletensis", "JQ993657_Eimeria_sp._ex_Apodemus_agrarius")
ver <- c(phyloT$tip.label[grep("vermiformis",phyloT$tip.label)], "KU174465_Eimeria_apionodes", "KU174467_Eimeria_apionodes")

species_ASV <- function(fer, name) {
    phyloT <- ape::read.tree("tmp/siteCF.cf.tree")
    clade.nr <- MRCA(phyloT, fer)
    ferT <- extract.clade(phyloT, clade.nr)
    ferrisiASV <- which(names18S%in%ferT$tip.label)
    amplicon$species[ferrisiASV] <- name
    amplicon$species
}

amplicon$species <- species_ASV(fer, "ferrisi")

amplicon$species <- species_ASV(ver, "vermiformis")

amplicon$species <- species_ASV(fal, "falciformis")

amplicon$species[grep("D3A", amplicon$names18S)] <- "28S"
                 
amplicon

### manual labeling based on phylogenetic tree
#amp5 <- c("ferrisi", "ferrisi", "Sp", "falciformis", "vermiformis")
#amp13 <- c("28S", "28S", "28S", "28S", "28S", "28S", "28S", "28S")
#amp16 <- c("ferrisi", "falciformis", "vermiformis")
#amp23 <- c("ferrisi", "Sp", "ferrisi", "ferrisi", "ferrisi", "Sp", "Sp")
#amp28 <- c("ferrisi", "ferrisi", "Sp", "ferrisi", "ferrisi")
#amp30 <- c("ferrisi", "falciformis")
#amp31 <- "ferrisi"
#amp33 <- c("ferrisi", "falciformis")
#amp35 <- c("ferrisi", "falciformis")

#falciformis
#which(names18S%in%c("18S_0067a_deg_3Mod_53_F.NSR399_3Mod_53_R_ASV_4", "Prot1702_32_F.wang1624CR6S_16_R_ASV_2", "Proti15_25_F.Proti440R_28_R_ASV_2", "18S_0067a_deg_5Mod_52_F.NSR399_5Mod_52_R_ASV_2", "Nem_0425_6CR_14_F.EukB_MH_24_R_ASV_2"))

#vermiformis
#which(names18S%in%c("18S_0067a_deg_3Mod_53_F.NSR399_3Mod_53_R_ASV_5", "18S_0067a_deg_5Mod_52_F.NSR399_5Mod_52_R_ASV_3"))

#sp
#which(names18S%in%c("18S_0067a_deg_3Mod_53_F.NSR399_3Mod_53_R_ASV_3", "MarkN_10_F.Euk360_CR_21_R_ASV_7", "MarkN_10_F.Euk360_CR_21_R_ASV_2", "MarkN_10_F.Euk360_CR_21_R_ASV_6", "MarkN_10_F.Proti440R_28_R_ASV_3"))

#ferrisi
#which(names18S%in%c("Proti15_25_F.Proti440R_28_R_ASV_1", "MarkN_10_F.Proti440R_28_R_ASV_1", "MarkN_10_F.Euk360_CR_21_R_ASV_5", "wang1141_13_F.Nem_0425_6_3_R_ASV_1", "Prot1702_32_F.wang1624CR6S_16_R_ASV_1", "Nem_0425_6CR_14_F.EukB_MH_24_R_ASV_1", "18S_0067a_deg_5Mod_52_F.NSR399_5Mod_52_R_ASV_1", "18S_0067a_deg_3Mod_53_F.NSR399_3Mod_53_R_ASV_1", "MarkN_10_F.Proti440R_28_R_ASV_5", "MarkN_10_F.Euk360_CR_21_R_ASV_3", "MarkN_10_F.Proti440R_28_R_ASV_4", "MarkN_10_F.Euk360_CR_21_R_ASV_1", "MarkN_10_F.Proti440R_28_R_ASV_2", "MarkN_10_F.Euk360_CR_21_R_ASV_4", "18S_0067a_deg_3Mod_53_F.NSR399_3Mod_53_R_ASV_2"))

#amplicon$spec <- c(amp5, amp13, amp16, amp23, amp28, amp30, amp31, amp33, amp35)

#amp.8 <- c("ferrisi", "ferrisi", "ferrisi", "ferrisi", "sp.", "sp.", "vermiformis", "vermiformis")
#amp.12 <- c("ferrisi", "ferrisi", "ferrisi", "ferrisi")
#amp.22 <- rep("28S", 10)
#amp.25 <- c("ferrisi", "sp.", "vermiformis")
#amp.27 <- c("ferrisi", "falciformis")
#amp.33 <- c("ferrisi", "ferrisi", "ferrisi", "ferrisi", "sp.", "sp.")
#amp.35 <- c("ferrisi", "ferrisi")
#amp.41 <- c("ferrisi", "ferrisi", "ferrisi", "ferrisi", "ferrisi", "sp.")
#amp.43 <- c("falciformis", "ferrisi")
#amp.45 <- c("ferrisi", "ferrisi")
#amp.47 <- c("28S", "28S")
#amp.50 <- c("ferrisi", "falciformis", "ferrisi")
#amp.52 <- c("ferrisi", "falciformis")
#amp.55 <- c("ferrisi", "falciformis", "vermiformis")
#amp.59 <- c("ferrisi", "falciformis", "ferrisi", "falciformis")

amplicon


#amplicon$spec <- c(amp.8,amp.12, amp.22, amp.25, amp.27, amp.33, amp.35, amp.41, amp.43, amp.45, amp.47, amp.50, amp.52, amp.55, amp.59)

#amplicon[amplicon$names18S=="Proti643_31_F.Proti1580R_32_R_ASV_2",]

amplicon$names18S[amplicon$spec=="ferrisi"]

amplicon$names18S[amplicon$spec=="falciformis"]

amplicon$names18S[amplicon$spec=="sp"]

amplicon$names18S[amplicon$spec=="vermiformis"]

#sanity check
Eim@tax_table[,7] <- amplicon$spec
Eim.TSS@tax_table[,7] <- amplicon$spec

get_taxa_unique(Eim, "genus")

amp_names <- gsub("_ASV.*", "", names18S)
Eim@tax_table[,6] <- amp_names
Eim.TSS@tax_table[,6] <- amp_names

Eim.m <- Eim
Eim.TSS.m <- Eim.TSS

#Eim.m@tax_table[,6] <- amplicon$names18S
#Eim.TSS.m@tax_table[,6] <- amplicon$names18S

# removing empty samples
Eim <- phyloseq::prune_samples(sample_sums(Eim)>0, Eim)
Eim.TSS <- phyloseq::prune_samples(sample_sums(Eim.TSS)>0, Eim.TSS)

# separating 18S and 28S

amp_names

Eim.ASV <- DNAStringSet(taxa_names(Eim))

names(Eim.ASV) <- amplicon$names18S

Eim.ASV

Seq18 <- Eim.ASV[which(!grepl("D3A", names(Eim.ASV)))]
Seq28 <- Eim.ASV[which(grepl("D3A", names(Eim.ASV)))]

writeFasta(Eim.ASV, "tmp/EimeriaASV.fasta")
writeFasta(Seq18, "tmp/Eimeria18S.fasta")
writeFasta(Seq28, "tmp/Eimeria28S.fasta")


Eim18 <- Eim
Eim.TSS18 <- Eim.TSS
Eim28 <- Eim
Eim.TSS28 <- Eim.TSS
Eim18 <- subset_taxa(Eim18, !species=="28S")
Eim.TSS18 <- subset_taxa(Eim.TSS18, !species=="28S")
Eim28 <- subset_taxa(Eim28, species=="28S")
Eim.TSS28 <- subset_taxa(Eim.TSS28, species=="28S")

#Eim28@tax_table[10,6]

get_taxa_unique(Eim18, "species")

# now this is for plotting co-infections
Eim18 <- phyloseq::prune_samples(sample_sums(Eim18)>0, Eim18)
Eim.TSS18 <- phyloseq::prune_samples(sample_sums(Eim.TSS18)>0, Eim.TSS18)
#eim18.TSS <- transform_sample_counts(Eim18, function(x) x / sum(x)) 
Eim18_sp <- Eim.TSS18

Eim18_sp@tax_table[,6] <- Eim18_sp@tax_table[,5]

Eim18_sp <- tax_glom(Eim18_sp, taxrank="species")

Eim18_sp

#sanity check
colnames(Eim18_sp@otu_table)==rownames(Eim18_sp@tax_table)
#colnames(Eim18_sp@otu_table) <- Eim18_sp@tax_table[,5]




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
