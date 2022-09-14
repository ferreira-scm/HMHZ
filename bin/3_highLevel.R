# high level analysis of the "wild microbiome"

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

PS <- readRDS(file = "tmp/interData/PhyloSeqCombi_HMHZ_All.Rds")
PS.l <- readRDS(file = "tmp/interData/PhyloSeqList_HMHZ_All.Rds")

## filtering MA by amplicon
fPS.l <- list()
for (i in 1:48) {
    try(fPS.l[[i]] <- fil(PS.l[[i]]), silent=TRUE)
}

fPS <- fPS.l[[1]]
for (i in 2:44){
    fPS <- try(merge_phyloseq(fPS,fPS.l[[i]]))
#    print(fPS)
    }

# filter pooled for comparison
fPS.pool <- fil(PS)

##Eliminate Unassigned to superkingdom level
#PS <- subset_taxa(PS, !is.na(superkingdom) & !superkingdom %in% c("", "uncharacterized"))

# subset samples based on total read count (1000 reads)
median(phyloseq::sample_sums(PS))

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

# How many Eimeria asv's do we have?
subset_taxa(PS, genus%in%"Eimeria") # 218 before filtering

PSeim <- subset_taxa(fPS, genus%in%"Eimeria")

PSeim

PSeim@tax_table[,6]

PS.l

fPS.l

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
fPS.l <- list()
for (i in 1:48) {
    try(fPS.l[[i]] <- fil(PS.l[[i]]), silent=TRUE)
}

names(fPS.l) <- names(PS.l)

### and what happens when we filter?
nmEim <- list()
names18S <- list()

for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(fPS.l[[i]],genus=="Eimeria"), silent=TRUE)
#    try(get_taxa_unique(p, "genus"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "genus")
        print(paste(i, "- ", names(fPS.l[i]), ": ", nrow(p@tax_table), sep=""))
                                        #        print(a)
        nmEim[i] <- names(fPS.l[i])
        names18S[[i]] <- paste(rep(names(fPS.l[i]), nrow(p@tax_table)), "ASV", seq(1, nrow(p@tax_table), 1), sep="_") 
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
Eim <- subset_taxa(fPS, genus%in%"Eimeria")
#Eim <- prune_taxa(!rownames(Eim@tax_table)%in%ASV.28S, Eim)

Eim

library(ShortRead)
## I also want to store the amplicon info to know which ASV's came from where
EimASV <- DNAStringSet(rownames(Eim@tax_table))

names(EimASV) <- names18S

writeFasta(EimASV, "tmp/EimeriaASV.fasta")

EimASV <- EimASV[which(!grepl("D3A", names(EimASV)))]
EimASV <- EimASV[which(!grepl("NLF", names(EimASV)))]

writeFasta(EimASV, "tmp/Eimeria18S.fasta")

EimASV

#PS <- subset_taxa(PS, !genus %in% "Mus") ##Eliminate reads :S

#names(fPS.l)[which(!names(fPS.l)%in%primerL$Primer_name)]



PScri <- subset_taxa(PS, family%in%"Cryptosporidiidae")
#PSsar <- subset_taxa(PS, family%in%"Sarcocystidae")
# for nematoda
PSoxy <- subset_taxa(PS, family%in%"Oxyuridae")
PStri <- subset_taxa(PS, family%in%"Trichuridae")
PShex <- subset_taxa(PS, family%in%"Heteroxynematidae")
#PSstr <- subset_taxa(PS, family%in%"Strongyloididae")
PSasc <- subset_taxa(PS, family%in%"Ascaridiidae")
PSspi <- subset_taxa(PS, family%in%"Spirocercidae")
PShek <- subset_taxa(PS, family%in%"Heterakidae")
# for platyhelminth
PShym <- subset_taxa(PS, family%in%"Hymenolepididae")

sample_data(PS)$eimASV <- sample_sums(PSeim)
sample_data(PS)$criASV <- sample_sums(PScri)
sample_data(PS)$oxyASV <- sample_sums(PSoxy)
sample_data(PS)$triASV <- sample_sums(PStri)
sample_data(PS)$hexASV <- sample_sums(PShex)
sample_data(PS)$ascASV <- sample_sums(PSasc)
sample_data(PS)$hekASV <- sample_sums(PShek)
sample_data(PS)$spiASV <- sample_sums(PSspi)
sample_data(PS)$hymASV <- sample_sums(PShym)


### now I subset based on parasites for downstream analyses
# first at phyla level
PSnem <- subset_taxa(pPS, phylum%in%"Nematoda")
PSapi <- subset_taxa(pPS, phylum%in%"Apicomplexa")
PSpla <- subset_taxa(pPS, phylum%in%"Platyhelminthes")
#for apicomplexa
PSeim <- subset_taxa(pPS, family%in%"Eimeriidae")
PScri <- subset_taxa(pPS, family%in%"Cryptosporidiidae")
PSsar <- subset_taxa(PS, family%in%"Sarcocystidae")
# for nematoda
PSoxy <- subset_taxa(pPS, family%in%"Oxyuridae")
PStri <- subset_taxa(pPS, family%in%"Trichuridae")
PShex <- subset_taxa(pPS, family%in%"Heteroxynematidae")
PSstr <- subset_taxa(PS, family%in%"Strongyloididae")
PSasc <- subset_taxa(pPS, family%in%"Ascaridiidae")
PSspi <- subset_taxa(pPS, family%in%"Spirocercidae")
PShek <- subset_taxa(pPS, family%in%"Heterakidae")
# for platyhelminth
PShym <- subset_taxa(pPS, family%in%"Hymenolepididae")
PSano <- subset_taxa(PS, family%in%"Anoplocephalidae")

TruthN <- sample_sums(PSnem)>=1
TruthA <- sample_sums(PSapi)>=1
Trutheim <- sample_sums(PSeim)>=1 # 475
Truthcri <- sample_sums(PScri)>=1 # 54 samples
#Truthsar <- sample_sums(PSsar)>=1 # 5 samples
Truthoxy <- sample_sums(PSoxy)>=1 # 352 samples
Truthtri <- sample_sums(PStri)>=1 # 241 samples
Truthhex <- sample_sums(PShex)>=1 # 356 samples
Truthstr <- sample_sums(PSstr)>=1 # 2 samples
Truthasc <- sample_sums(PSasc)>=1 # 23
Truthspi <- sample_sums(PSspi)>=1 #12 samples
Truthhek <- sample_sums(PShek)>=1 #24 samples
Truthhym <- sample_sums(PShym)>=1 #8 samples
Truthano <- sample_sums(PSano)>=1 #1 sample


# make a new variable called parasite (yes/no) if nematode or apicomplexa are present
for (i in 1:nsamples(pPS))
{ if (TruthN[i]==TRUE)
  {sample_data(pPS)$NemASV[i]="yes"}
  else {sample_data(pPS)$NemASV[i]="no"}
}
for (i in 1:nsamples(pPS))
{ if (TruthA[i]==TRUE)
  {sample_data(pPS)$ApiASV[i]="yes"}
  else {sample_data(pPS)$ApiASV[i]="no"}
}
# apicomplexa
for (i in 1:nsamples(pPS))
{ if (Trutheim[i]==TRUE)
  {sample_data(pPS)$peimASV[i]="yes"}
  else {sample_data(pPS)$peimASV[i]="no"}
}
for (i in 1:nsamples(pPS))
{ if (Truthcri[i]==TRUE)
  {sample_data(pPS)$pcriASV[i]="yes"}
  else {sample_data(pPS)$pcriASV[i]="no"}
}
#nematoda
for (i in 1:nsamples(pPS))
{ if (Truthoxy[i]==TRUE)
  {sample_data(pPS)$poxyASV[i]="yes"}
  else {sample_data(pPS)$poxyASV[i]="no"}
}
for (i in 1:nsamples(pPS))
{ if (Truthtri[i]==TRUE)
  {sample_data(pPS)$ptriASV[i]="yes"}
  else {sample_data(pPS)$ptriASV[i]="no"}
}
for (i in 1:nsamples(pPS))
{ if (Truthhex[i]==TRUE)
  {sample_data(pPS)$phexASV[i]="yes"}
  else {sample_data(pPS)$phexASV[i]="no"}
}
for (i in 1:nsamples(pPS))
{ if (Truthasc[i]==TRUE)
  {sample_data(pPS)$pascASV[i]="yes"}
  else {sample_data(pPS)$pascASV[i]="no"}
}
for (i in 1:nsamples(pPS))
{ if (Truthhek[i]==TRUE)
  {sample_data(pPS)$phekASV[i]="yes"}
  else {sample_data(pPS)$phekASV[i]="no"}
}
#Platyhelminths
for (i in 1:nsamples(pPS))
{ if (Truthhym[i]==TRUE)
  {sample_data(pPS)$phymASV[i]="yes"}
  else {sample_data(pPS)$phymASV[i]="no"}
}

# let's consider positive for Eimeria when it has Eimeria reads and oocysts
#sample_data(pPS)$peimASV[na.omit(sample_data(pPS)$OPG>0)]<- "yes"
#sample_data(pPS)$pcriASV[na.omit(sample_data(pPS)$Crypto_Oocyst_calculated>0)] <- "yes"
# same for nematodes
#sample_data(pPS)$poxyASV[na.omit(sample_data(pPS)$Pinworms>0)] <- "yes"
#sample_data(pPS)$pasc[na.omit(sample_data(pPS)$Aspiculuris_tetraptera>0)] <- "yes"
#sample_data(pPS)$poxyASV[na.omit(sample_data(pPS)$Syphacia_obvelata>0)] <- "yes"
#sample_data(pPS)$phekASV[na.omit(sample_data(pPS)$Heterakis_spumosa>0)] <- "yes"
#sample_data(pPS)$phymASV[na.omit(sample_data(pPS)$Hymenolepis_microstoma>0)] <- "yes"
#sample_data(pPS)$ptriASV[na.omit(sample_data(pPS)$Trichuris_muris>0)] <- "yes"

# at genus level for Pintworms?
#summary(sample_data(PS)$Taenia_martis>0)
#summary(sample_data(PS)$Catenotaenia_pusilla>0)
#summary(sample_data(PS)$Hymenolepis_diminuta>0)
#summary(sample_data(PS)$Hymenolepis_microstoma>0)
#summary(sample_data(PS)$Mastophorus>0)

###### categorize intensity level
Trutheim <- sample_sums(PSeim)>mean(sample_sums(PSeim))
Truthcri <- sample_sums(PScri)>mean(sample_sums(PScri))
#Truthsar <- sample_sums(PSsar)>=1 # 5 samples
Truthoxy <- sample_sums(PSoxy)>mean(sample_sums(PSoxy))
Truthtri <- sample_sums(PStri)>mean(sample_sums(PStri))
Truthhex <- sample_sums(PStri)>mean(sample_sums(PStri))
#Truthstr <- sample_sums(PSstr)>=1 # 2 samples
Truthasc <- sample_sums(PSasc)>mean(sample_sums(PSasc))
#Truthspi <-sample_sums(PStri)>mean(sample_sums(PStri))
Truthhek <- sample_sums(PShek)>mean(sample_sums(PShek))
Truthhym <- sample_sums(PShym)>mean(sample_sums(PShym))
# make a new variable called parasite (yes/no) if nematode or apicomplexa are present
# apicomplexa
for (i in 1:nsamples(pPS))
{ if (Trutheim[i]==TRUE)
  {sample_data(pPS)$iEimASV[i]="high"}
  else {sample_data(pPS)$iEimASV[i]="low"}
}
for (i in 1:nsamples(pPS))
{ if (Truthcri[i]==TRUE)
  {sample_data(pPS)$icriASV[i]="high"}
  else {sample_data(pPS)$icriASV[i]="low"}
}
#nematoda
for (i in 1:nsamples(pPS))
{ if (Truthoxy[i]==TRUE)
  {sample_data(pPS)$ioxyASV[i]="high"}
  else {sample_data(pPS)$ioxyASV[i]="low"}
}
for (i in 1:nsamples(pPS))
{ if (Truthtri[i]==TRUE)
  {sample_data(pPS)$itriASV[i]="high"}
  else {sample_data(pPS)$itriASV[i]="low"}
}
for (i in 1:nsamples(pPS))
{ if (Truthhex[i]==TRUE)
  {sample_data(pPS)$ihexASV[i]="high"}
  else {sample_data(pPS)$ihexASV[i]="low"}
}
for (i in 1:nsamples(pPS))
{ if (Truthasc[i]==TRUE)
  {sample_data(pPS)$iascASV[i]="high"}
  else {sample_data(pPS)$iascASV[i]="low"}
}
for (i in 1:nsamples(pPS))
{ if (Truthhek[i]==TRUE)
  {sample_data(pPS)$ihekASV[i]="high"}
  else {sample_data(pPS)$ihekASV[i]="low"}
}
#Platyhelminths
for (i in 1:nsamples(pPS))
{ if (Truthhym[i]==TRUE)
  {sample_data(pPS)$ihymASV[i]="high"}
  else {sample_data(pPS)$ihymASV[i]="low"}
}

sample_data(pPS)$eimASV <- sample_sums(PSeim)
sample_data(pPS)$criASV <- sample_sums(PScri)
sample_data(pPS)$oxyASV <- sample_sums(PSoxy)
sample_data(pPS)$triASV <- sample_sums(PStri)
sample_data(pPS)$hexASV <- sample_sums(PShex)
sample_data(pPS)$ascASV <- sample_sums(PSasc)
sample_data(pPS)$hekASV <- sample_sums(PShek)
sample_data(pPS)$hymASV <- sample_sums(PShym)
sample_data(pPS)$spiASV <- sample_sums(PSspi)

# transform and agglomerate

sample_data(pPS)$LibSize <- sample_sums(pPS)

sample_data(PS)$LibSize <- sample_sums(PS)

sample_sums(pPS)

                                        # save
saveRDS(PS, file="tmp/PSpre.R")
saveRDS(pPS, file="tmp/pPSpre.R")

## checking apicomplexa
PSapi <- subset_taxa(pPS, phylum%in%"Apicomplexa")
bra_nmds <- ordinate(PSapi, "NMDS", "bray")
plot_ordination(PSapi, ordinate(PSapi, "NMDS", "bray"))
