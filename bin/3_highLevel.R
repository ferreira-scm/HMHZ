# high level analysis of the "wild microbiome"
# preprocessing

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

PSl11 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist11.Rds")
PSl12 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist12.Rds")
PSl21 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist21.Rds")
PSl22 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist22.Rds")

###Eliminate the empty primers first and then merge the phyloseq lists... empty list make the next function bug
PSl11[sapply(PSl11, is.null)]<- NULL
PSl12[sapply(PSl12, is.null)]<- NULL
PSl21[sapply(PSl21, is.null)]<- NULL
PSl22[sapply(PSl22, is.null)]<- NULL
 ##Merge all the information from both experiments
along<- names(PSl11)
PSl <- lapply(along, function(i) merge_phyloseq(PSl11[[i]], PSl12[[i]],PSl21[[i]], PSl22[[i]]))
names(PSl) <- names(PSl11)

# adjusting amplicon names for some reason
x<- names(PSl)
x[6]<-"BGf_132_F.BGr_132_R"
x[22]<-"LSU_Fwd_2_3Mod_55_F.LSU_Rev_4_54_R"
x[28]<-"NLF184cw_74_F.NL818cw_74_R"
names(PSl)<- x

saveRDS(PSl, file="tmp/PSl.Rds")
saveRDS(PS, file="tmp/PSCombi.Rds")

PSl <- readRDS(file="tmp/PSl.Rds")

# load qiime-derived PS

PSq <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/PSqiime200.RDS")
PSq
PS

rm(PS11)
rm(PS12)
rm(PS21)
rm(PS22)
rm(PSl11)
rm(PSl12)
rm(PSl21)
rm(PSl22)
rm(along)

##Eliminate Unassigned to superkingdom level
PS <- subset_taxa(PS, !is.na(superkingdom) & !superkingdom %in% c("", "uncharacterized"))

PSl <- lapply(PSl, function(x){
    subset_taxa(x, !is.na(superkingdom) & !superkingdom %in% c("", "uncharacterized"))})

    
# subset samples based on total read count (1000 reads)
median(phyloseq::sample_sums(PS))

PS <- phyloseq::subset_samples(PS, phyloseq::sample_sums(PS) > 1000)

# can't do it for the PSl
#lapply(PSl, function(x){
#    subset_samples(x, phyloseq::sample_sums(x)>50
#})

# prevalence filtering at 0.05%
#PS=phyloseq_filter_prevalence(PS, prev.trh=0.005)

###A lot of Mus :(
## Host read numbers
sum(otu_table(subset_taxa(PS, genus%in%"Mus")))/sum(otu_table(PS))
###Eliminate reads assigned as "Mus"
PS <- subset_taxa(PS, !genus %in% "Mus") ##Eliminate reads :S

PSl <- lapply(PSl, function(x){
     subset_taxa(x, !genus %in% "Mus")
})


# Eliminate samples with no reads
PS <- prune_samples(sample_sums(PS)>0, PS)

# abundance filtering to 0.01%? Or keep prevalence filtering?
x = taxa_sums(PS)
keepTaxa = (x / sum(x) > 0.0001)
summary(keepTaxa)
pPS = prune_taxa(keepTaxa, PS)

# Eliminate samples with no reads
pPS <- prune_samples(sample_sums(pPS)>0, pPS)


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
#PSsar <- subset_taxa(PS, family%in%"Sarcocystidae")
# for nematoda
PSoxy <- subset_taxa(pPS, family%in%"Oxyuridae")
PStri <- subset_taxa(pPS, family%in%"Trichuridae")
PShex <- subset_taxa(pPS, family%in%"Heteroxynematidae")
#PSstr <- subset_taxa(PS, family%in%"Strongyloididae")
PSasc <- subset_taxa(pPS, family%in%"Ascaridiidae")
PSspi <- subset_taxa(pPS, family%in%"Spirocercidae")
PShek <- subset_taxa(pPS, family%in%"Heterakidae")
# for platyhelminth
PShym <- subset_taxa(pPS, family%in%"Hymenolepididae")
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
