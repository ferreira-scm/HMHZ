# high level analysis of the "wild microbiome"

library(ggplot2)
library(reshape)
library(phyloseq, lib="/usr/local/lib/R/site-library/")
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

# merge all the runs
PS11 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi11.Rds")
PS12 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi12.Rds")
PS21 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi21.Rds")
PS22 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi22.Rds")
PS <- merge_phyloseq(PS11, PS12, PS21, PS22)

##Eliminate Unassigned to superkingdom level
PS <- subset_taxa(PS, !is.na(superkingdom) & !superkingdom %in% c("", "uncharacterized"))

#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/PSCombi.Rds")

# subset samples based on total read count
sort(phyloseq::sample_sums(PS))
PS <- phyloseq::subset_samples(PS, phyloseq::sample_sums(PS) > 1000)

# abundance filtering to 0.01%? Or keep prevalence filtering?

# prevalence filtering at 0.05%
PS=phyloseq_filter_prevalence(PS, prev.trh=0.005)

###A lot of Mus :(
## Host read numbers
sum(otu_table(subset_taxa(PS, genus%in%"Mus")))/sum(otu_table(PS))
###Eliminate reads assigned as "Mus"
PS <- subset_taxa(PS, !genus %in% "Mus") ##Eliminate reads :S

# Eliminate samples with no reads
PS <- prune_samples(sample_sums(PS)>0, PS)


# transform and agglomerate
PSphy <-  tax_glom(PS, "phylum", NArm = FALSE)
sample_data(PS)$depth <- sample_sums(PS)
PSTSS = transform_sample_counts(PS, function(x){x / sum(x)})
PSTSS1 = transform_sample_counts(PS, function(x){(x / sum(x))*10000})

### now I subset based on parasites for downstream analyses
# first at phyla level
PSnem <- subset_taxa(PS, phylum%in%"Nematoda")
PSapi <- subset_taxa(PS, phylum%in%"Apicomplexa")
PSpla <- subset_taxa(PS, phylum%in%"Platyhelminthes")
#for apicomplexa
PSeim <- subset_taxa(PS, family%in%"Eimeriidae")
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
#PSano <- subset_taxa(PS, family%in%"Anoplocephalidae")

TruthN <- sample_sums(PSnem)>=1
TruthA <- sample_sums(PSapi)>=1
Trutheim <- sample_sums(PSeim)>=1 # 475
Truthcri <- sample_sums(PScri)>=1 # 54 samples
#Truthsar <- sample_sums(PSsar)>=1 # 5 samples
Truthoxy <- sample_sums(PSoxy)>=1 # 352 samples
Truthtri <- sample_sums(PStri)>=1 # 241 samples
Truthhex <- sample_sums(PShex)>=1 # 356 samples
#Truthstr <- sample_sums(PSstr)>=1 # 2 samples
Truthasc <- sample_sums(PSasc)>=1 # 23
#Truthspi <- sample_sums(PSspi)>=1 #12 samples
Truthhek <- sample_sums(PShek)>=1 #24 samples
Truthhym <- sample_sums(PShym)>=1 #8 samples
#Truthano <- sample_sums(PSano)>=1 1 sample

# make a new variable called parasite (yes/no) if nematode or apicomplexa are present
for (i in 1:nsamples(PS))
{ if (TruthN[i]==TRUE)
  {sample_data(PS)$NemASV[i]="yes"}
  else {sample_data(PS)$NemASV[i]="no"}
}
for (i in 1:nsamples(PS))
{ if (TruthA[i]==TRUE)
  {sample_data(PS)$ApiASV[i]="yes"}
  else {sample_data(PS)$ApiASV[i]="no"}
}
# apicomplexa
for (i in 1:nsamples(PS))
{ if (Trutheim[i]==TRUE)
  {sample_data(PS)$EimASV[i]="yes"}
  else {sample_data(PS)$EimASV[i]="no"}
}
for (i in 1:nsamples(PS))
{ if (Truthcri[i]==TRUE)
  {sample_data(PS)$criASV[i]="yes"}
  else {sample_data(PS)$criASV[i]="no"}
}
#nematoda
for (i in 1:nsamples(PS))
{ if (Truthoxy[i]==TRUE)
  {sample_data(PS)$oxyASV[i]="yes"}
  else {sample_data(PS)$oxyASV[i]="no"}
}
for (i in 1:nsamples(PS))
{ if (Truthtri[i]==TRUE)
  {sample_data(PS)$triASV[i]="yes"}
  else {sample_data(PS)$triASV[i]="no"}
}
for (i in 1:nsamples(PS))
{ if (Truthhex[i]==TRUE)
  {sample_data(PS)$hexASV[i]="yes"}
  else {sample_data(PS)$hexASV[i]="no"}
}
for (i in 1:nsamples(PS))
{ if (Truthasc[i]==TRUE)
  {sample_data(PS)$ascASV[i]="yes"}
  else {sample_data(PS)$ascASV[i]="no"}
}
for (i in 1:nsamples(PS))
{ if (Truthhek[i]==TRUE)
  {sample_data(PS)$hekASV[i]="yes"}
  else {sample_data(PS)$hekASV[i]="no"}
}
#Platyhelminths
for (i in 1:nsamples(PS))
{ if (Truthhym[i]==TRUE)
  {sample_data(PS)$hymASV[i]="yes"}
  else {sample_data(PS)$hymASV[i]="no"}
}

# let's consider positive for Eimeria when it has Eimeria reads and oocysts
sample_data(PS)$EimASV[na.omit(sample_data(PS)$OPG>0)] <- "yes"
sample_data(PS)$criASV[na.omit(sample_data(PS)$Crypto_Oocyst_calculated>0)] <- "yes"

# same for nematodes
sample_data(PS)$NemASV[na.omit(sample_data(PS)$nema>0)] <- "yes"
sample_data(PS)$oxyASV[na.omit(sample_data(PS)$Pinworms>0)] <- "yes"
sample_data(PS)$oxyASV[na.omit(sample_data(PS)$Aspiculuris_tetraptera>0)] <- "yes"
sample_data(PS)$oxyASV[na.omit(sample_data(PS)$Syphacia_obvelata>0)] <- "yes"
sample_data(PS)$hekASV[na.omit(sample_data(PS)$Heterakis_spumosa>0)] <- "yes"
sample_data(PS)$oxyASV[na.omit(sample_data(PS)$Syphacia_obvelata>0)] <- "yes"
sample_data(PS)$hymASV[na.omit(sample_data(PS)$Hymenolepis_microstoma>0)] <- "yes"
sample_data(PS)$triASV[na.omit(sample_data(PS)$Trichuris_muris>0)] <- "yes"

# at genus level for Pintworms?

#summary(sample_data(PS)$Taenia_martis>0)
#summary(sample_data(PS)$Catenotaenia_pusilla>0)
#summary(sample_data(PS)$Hymenolepis_diminuta>0)
#summary(sample_data(PS)$Hymenolepis_microstoma>0)
#summary(sample_data(PS)$Mastophorus>0)

###### categorize intensity level

Trutheim <- sample_sums(PSeim)>mean(sample_sums(PSeim))

Truthcri <- sample_sums(PScri)>mean(sample_sums(PScri))


Truthcri

#Truthsar <- sample_sums(PSsar)>=1 # 5 samples
Truthoxy <- sample_sums(PSoxy)>mean(sample_sums(PSoxy))
Truthtri <- sample_sums(PStri)>mean(sample_sums(PStri))
Truthhex <- sample_sums(PStri)>mean(sample_sums(PStri))
#Truthstr <- sample_sums(PSstr)>=1 # 2 samples
Truthasc <- sample_sums(PSasc)>mean(sample_sums(PSasc))
#Truthspi <-sample_sums(PStri)>mean(sample_sums(PStri))
Truthhek <- sample_sums(PShek)>mean(sample_sums(PShek))
Truthhym <- sample_sums(PShym)>mean(sample_sums(PShym))
#Truthano <-

# make a new variable called parasite (yes/no) if nematode or apicomplexa are present
# apicomplexa
for (i in 1:nsamples(PS))
{ if (Trutheim[i]==TRUE)
  {sample_data(PS)$iEimASV[i]="high"}
  else {sample_data(PS)$iEimASV[i]="low"}
}

for (i in 1:nsamples(PS))
{ if (Truthcri[i]==TRUE)
  {sample_data(PS)$icriASV[i]="high"}
  else {sample_data(PS)$icriASV[i]="low"}
}

                                        #nematoda
for (i in 1:nsamples(PS))
{ if (Truthoxy[i]==TRUE)
  {sample_data(PS)$ioxyASV[i]="high"}
  else {sample_data(PS)$ioxyASV[i]="low"}
}
for (i in 1:nsamples(PS))
{ if (Truthtri[i]==TRUE)
  {sample_data(PS)$itriASV[i]="high"}
  else {sample_data(PS)$itriASV[i]="low"}
}
for (i in 1:nsamples(PS))
{ if (Truthhex[i]==TRUE)
  {sample_data(PS)$ihexASV[i]="high"}
  else {sample_data(PS)$ihexASV[i]="low"}
}
for (i in 1:nsamples(PS))
{ if (Truthasc[i]==TRUE)
  {sample_data(PS)$iascASV[i]="high"}
  else {sample_data(PS)$iascASV[i]="low"}
}
for (i in 1:nsamples(PS))
{ if (Truthhek[i]==TRUE)
  {sample_data(PS)$ihekASV[i]="high"}
  else {sample_data(PS)$ihekASV[i]="low"}
}
#Platyhelminths
for (i in 1:nsamples(PS))
{ if (Truthhym[i]==TRUE)
  {sample_data(PS)$ihymASV[i]="high"}
  else {sample_data(PS)$ihymASV[i]="low"}
}

sample_data(PS)$rEim[i]=sample_sums(PSeim)

sample_data(PS)$rcri[i]=sample_sums(PScri)

sample_data(PS)$rtri[i]=sample_sums(PStri)

sample_data(PS)$rhex[i]=sample_sums(PShex)

sample_data(PS)$rasc[i]=sample_sums(PSasc)

sample_data(PS)$rhek[i]=sample_sums(PShek)

sample_data(PS)$rhym[i]=sample_sums(PShym)

saveRDS(PS, file="tmp/PSpre.R")
