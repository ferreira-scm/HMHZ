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

source("bin/Pre_Functions.R")

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

# How many Eimeria asv's do we have?
subset_taxa(PS, genus%in%"Eimeria") # 218 before filtering

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



#Eim <- prune_taxa(!rownames(Eim@tax_table)%in%ASV.28S, Eim)

## Now the alignment
library(ShortRead)
## I also want to store the amplicon info to know which ASV's came from where
EimASV <- DNAStringSet(rownames(Eim@tax_table))
names(EimASV) <- names18S
writeFasta(EimASV, "tmp/EimeriaASV.fasta")
EimASV <- EimASV[which(!grepl("D3A", names(EimASV)))]
EimASV <- EimASV[which(!grepl("NLF", names(EimASV)))]
writeFasta(EimASV, "tmp/Eimeria18S.fasta")

library(DECIPHER)
EimSeqs <- readFasta("../LabMicrobiome/tmp/Eimeria_seqs.fa")
name_seqs <- EimSeqs@id
EimSeqs <- sread(EimSeqs)
names(EimSeqs) <- name_seqs
allEim <- c(EimSeqs, EimASV)
Eim.align <- AlignSeqs(allEim, anchor=NA, iterations=10, refinements=10, processors=90)
Eim.align <-AdjustAlignment(Eim.align)
#names(Eim.align)[c(1:205)] <- name_seqs

writeFasta(Eim.align, "tmp/Eimeria18S_wild_lab_ref.fasta")
writeFasta(allEim, "tmp/Eimeria18S_wild_lab_ref_notAligned.fasta")

# and this is how I did the tree outside R
#~/iqtree-2.2.0-Linux/bin/iqtree2 -s tmp/Eimeria18S_wild_lab_ref.fasta -m MFP
#~/iqtree-2.2.0-Linux/bin/iqtree2 -s  tmp/Eimeria18S_wild_lab_ref.fasta -m TN+F+R2 -B 5000 -redo 

### Now we need  to analyse which ASV's are likely from the same Eimeira.
## ASV's within the same sample from different amplicons that correlate well.

# let's start by selecting a few samples with lots of Eimeria

histogram(phyloseq::sample_sums(Eim)) 

Eim@tax_table[,6] <- names18S

Eim.lots <- prune_samples(phyloseq::sample_sums(Eim)>0.2, Eim)
# removing empty ASV's



eim <- psmelt(Eim.lots)

for (i in rownames(Eim.lots@otu_table)){
eim %>%
    filter(Sample==i) %>%
    select(Abundance, species) %>%
    filter(Abundance>0) -> x
print(i)
print(x)
}

amplicon <- data.frame(names18S)

names18S

### manual labeling based on phylogenetic tree
amp.8 <- c("ferrisi", "ferrisi", "ferrisi", "ferrisi", "na_vermiformis", "na_vermiformis", "vermiformis", "vermiformis")
amp.12 <- c("chimera", "chimera", "chimera", "chimera")
amp.22 <- rep("28S", 10)
amp.25 <- c("ferrisi", "na_vermiformis", "vermiformis")
amp.27 <- c("ferrisi", "falciformis")
amp.33 <- c("ferrisi", "ferrisi", "ferrisi", "ferrisi", "na_vermiformis", "ferrisi")
amp.35 <- c("ferrisi", "ferrisi")
amp.41 <- c("ferrisi", "ferrisi", "na_vermiformis", "ferrisi", "ferrisi", "ferrisi")
amp.43 <- c("falciformis", "ferrisi")
amp.45 <- c("chimera", "chimera")
amp.47 <- c("28S", "28S")
amp.50 <- c("ferrisi", "falciformis", "ferrisi")
amp.52 <- c("ferrisi", "falciformis")
amp.55 <- c("ferrisi", "vermiformis", "vermiformis")
amp.59 <- c("ferrisi", "falciformis", "ferrisi", "falciformis")


amplicon$spec <- c(amp.8,amp.12, amp.22, amp.25, amp.27, amp.33, amp.35, amp.41, amp.43, amp.45, amp.47, amp.50, amp.52, amp.55, amp.59)

#sanity check
Eim@tax_table[,6]==amplicon$names18S

Eim@tax_table[,5] <- amplicon$spec

# removing empty samples
Eim <- phyloseq::prune_samples(sample_sums(Eim)>0, Eim)

eim.TSS <- transform_sample_counts(Eim, function(x) x / sum(x)) 

Eim_sp <- tax_glom(eim.TSS, "genus")

colnames(Eim_sp@otu_table)==rownames(Eim_sp@tax_table)
colnames(Eim_sp@otu_table) <- Eim_sp@tax_table[,5]

eim_sp <- psmelt(Eim_sp)
eim <- psmelt(eim.TSS)
eim2 <- psmelt(Eim)

eim$genus <- as.factor(eim$genus)
                                        # relevel
dist_bc <- (vegdist(Eim.TSS@otu_table, method="bray"))

dist_bc2 <- (vegdist(Eim_sp@otu_table, method="bray"))
res2 <- pcoa(dist_bc2)
res <- pcoa(dist_bc)

EH_sort <- names(sort(res$vectors[,1]))

EH_sort2 <- names(sort(res2$vectors[,1]))

eim$Sample <- factor(eim$Sample, levels= EH_sort2)
#eim_sp$Sample <- factor(eim_sp$Sample, levels=EH_sort2)

EH_sort2==levels(eim$Sample)

com_plot <- ggplot(eim, aes(x=Sample, y=Abundance, fill=genus))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=c("#999999", "lightgray", "forestgreen", "dodgerblue4", "darksalmon", "darkred"))+
    labs(fill="Eimeria ASV species", x="Sample", y="Eimeria relative abundance")+
    theme()+
    theme(axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
      strip.background = element_rect(colour="black", fill="white"))
#    coord_flip()

library(viridis)
library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")

Eim_heat <- ggplot(eim_sp, aes(Sample, genus, fill=Abundance))+
    geom_tile()+
      labs(y="Eimeria ASV species", x="Sample", fill="Eimeria relative abundance")+
    scale_fill_gradientn(colours = pal)+
    theme()+
        theme(axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
      strip.background = element_rect(colour="black", fill="white"))
#    scale_fill_distiller(palette = "RdPu")
#    scale_fill_gradient2(low="#00AFBB", mid="#E7B800", high="#FC4E07", midpoint=0.5)
#    scale_fill_viridis(discrete=FALSE)

library(cowplot)

Eim_species_comp <- plot_grid(com_plot, Eim_heat, ncol=1, align="v", labels="auto")

ggsave("fig/Eimeria_heat_composition.pdf", Eim_species_comp, height=5, width=10, dpi=600)

phyloseq::plot_bar(eim.TSS, fill="genus")+
    geom_bar(stat="identity", position="stack")+
    scale_fill_manual(values=c("white", "darkgray", "forestgreen", "dodgerblue4", "deepskyblue", "darksalmon", "darkred"))

### calculate dissimilarity distances
library(vegan)

library(ape)

Eim.TSS <- eim.TSS


# sanity check
colnames(eim.TSS@otu_table)==rownames(Eim.TSS@tax_table)

asvN <- seq(1,59,1)
asvN <- paste("ASV", asvN, sep="")
colnames(Eim.TSS@otu_table) <- paste(asvN, Eim.TSS@tax_table[,5], sep="_")
colnames(Eim.TSS@otu_table)


biplot(res2, Eim_sp@otu_table)


otu.st <- apply(Eim.TSS@otu_table, 2, scale, center=TRUE, scale=TRUE)
otu.st2 <- apply(Eim_sp@otu_table, 2, scale, center=TRUE, scale=TRUE)

biplot(res2, otu.st2,
       xlabs=rep("o", 362)
       )

biplot(res, otu.st,
       xlabs=rep("o", 362),
       cex=c(0.7, 0.8))


sample_sums(Eim)

fPS

#phyloseq::plot_heatmap(fPS)

plot_net(Eim, "bray")

phyloseq::plot_richness(Eim)



pcoa_bc <- ordinate(eim.TSS, "PCoA", "bray")

pcoa.p <- plot_ordination(eim.TSS, pcoa_bc)

pcoa.p

pcoa.p$data[,2]

### 
# correlation matrix
library(Hmisc)
library(Matrix)

library(SpiecEasi)

library(igraph)

# change species name to track amplicon

Eim@tax_table[,6] <- names18S

eim <- (Eim@otu_table)
tax <- data.frame(Eim@tax_table) 

head(eim)

head(tax)

otu.cor <- rcorr(as.matrix(eim), type="spearman")

otu.pval <- forceSymmetric(otu.cor$P)

cor.p <- p.adjust(otu.pval, method="BH") 

# consider only significant and strong correlations
#cor.r1[which(!cor.r1 > 0.8 | !cor.r1 < -0.8)]=NA
#cor.p[which(cor.p>0.01)]=NA
otu.pval@x <- cor.p

sel.tax <- tax[rownames(otu.pval),,drop=FALSE]

#sanity check
all.equal(rownames(sel.tax), rownames(otu.pval))

p.yes <- otu.pval<0.05

p.yes

r.val = otu.cor$r # select all the correlation values

p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion

p.yes.r

## select asv based on rho
p.yes.r <- abs(p.yes.r)>0.25 # output is logical vector
p.yes.rr <- p.yes.r*r.val # use logical vector for subscripting.

p.yes.rr

adjm <- as.matrix(p.yes.rr)

adjm



colnames(adjm) <- names18S
rownames(adjm) <- names18S

net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE) 

V(net.grph)

plot(net.grph)

### spiec easi

# parallel multicores
pargs <- list(rep.num=1000, seed=10010, ncores=90, thresh=0.05)
## mb
#se.net <- spiec.easi(Eim, method="mb", pulsar.params=pargs)

saveRDS(se.net, "tmp/se.Eimnet.RDS")
#se.net <- readRDS("tmp/se.Eimnet.RDS")  

# looking at lambda path
se.net$select$stars$summary

# coding ASV/ARG nodes
net.ids <- names18S

#### now we improve our visualization
## I want weighted edges
bm=symBeta(getOptBeta(se.net), mode="maxabs")
diag(bm) <- 0
#weights <- Matrix::summary(t(bm))[,3] # includes negative weights
weights <- (1-Matrix::summary(t(bm))[,3])/2 # ort

netEim <- adj2igraph(Matrix::drop0(getRefit(se.net)),
                     edge.attr=list(weight=weights),
                     vertex.attr = list(name=net.ids))


#betaMat=as.matrix(symBeta(getOptBeta(se.10net)))
# we want positive edges to be green and negative to be red
edges <- E(netEim)
edge.colors=c()
for (e.index in 1:length(edges)){
    adj.nodes=ends(netEim, edges[e.index])
    xindex=which(net.ids==adj.nodes[1])
    yindex=which(net.ids==adj.nodes[2])
    beta=bm[xindex, yindex]
    if (beta>0){
        edge.colors=append(edge.colors, "#1B7837")
    }else if(beta<0){
        edge.colors=append(edge.colors, "#762A83")
    }
}
E(netEim)$color=edge.colors 

### defining attributes



plot(netEim,
#     layout=layout_with_fr(netARG),
     vertex.label=net.ids,
     frame.col="grey")
#     vertex.color=adjustcolor(colorb,0.8),

###################################################################
############################################# parasite analysis
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
