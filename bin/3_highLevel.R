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
library("ape")
library(magrittr)
library(Biostrings)
library(DECIPHER)
library(ShortRead)
library(RColorBrewer)
library(ape)
library(ggtree)

source("bin/Pre_Functions.R")

## Now the alignment
## Eimeria from 18S
EimASV <- readFasta("tmp/Eimeria18S.fasta")
name_EimASV <- EimASV@id
EimASV <- sread(EimASV)
names(EimASV) <- name_EimASV

EimSeqs <- readFasta("../LabMicrobiome/tmp/Eimeria_seqs.fa")
name_seqs <- EimSeqs@id
EimSeqs <- sread(EimSeqs)
names(EimSeqs) <- name_seqs
allEim <- c(EimSeqs, EimASV)

names(EimSeqs)

allEim


allEim <- OrientNucleotides(allEim)

set.seed(12345)
Eim.align <- AlignSeqs(allEim, anchor=NA, iterations=20, refinements=20, processors=90)

Eim.align <- AdjustAlignment(Eim.align)

writeFasta(Eim.align, "tmp/Eimeria18S_wild_lab_ref.fasta")
#writeFasta(allEim, "tmp/Eimeria18S_wild_lab_ref_notAligned.fasta")
#writeFasta(Wild.align, "tmp/Eimeria18S_wild.fasta")
writeFasta(allEim, "tmp/Eimeria18S.fasta")

## removing long branches
#run_treeshrink.py -t tmp/Eimeria18S_wild_lab_ref.fasta.contree -m per-gene --force

#longB <- read.table("tmp/Eimeria18S_wild_lab_ref.fasta_treeshrink/output.txt", header=TRUE, sep="\t")

#longB <- names(longB)[1:6] 

#FinalEimeria.lb <- FinalEimeria[-which(names(FinalEimeria)%in%longB)]

#Eim.alignlb <- AlignSeqs(FinalEimeria.lb, anchor=NA, iterations=20, refinements=20, processors=90)

#Eim.alignlb <- AdjustAlignment(Eim.alignlb)

#writeFasta(Eim.align, "tmp/Eimeria18S_lb.fasta")

# I need set 3 groups of reference sequences (ferrisi, vermiformis, falciformis)
fer <- c(names(EimSeqs)[grep("ferrisi",names(EimSeqs))],"AF246717_Eimeria_telekii", "AF311643_Eimeria_separata")
fal <- names(EimSeqs)[grep("falciformis|JQ993650|KU174466|KU174449|KU174463|KU174481|JQ993657|JQ993666|KU174452|KU174470|KU174472|KU174473|KU174471|KU174460|KU174450|KU174455|KU174451|KU174456|KU174474|KU174461|KU174454|KU174458|KU174464|AF311644|KU174457|KU174453|KU174469|KU174476|KU174486|KU174478|KU174484|KU192956|JQ993665",names(EimSeqs))]
ver <- names(EimSeqs)[grep("vermiformis|KU192931|KU192936|KU174467|KU192916|KU174465",names(EimSeqs))]
#sp1 <-names(EimSeqs)[grep("U40263|KU174483|KU192965|KU192958|KU192961|JQ993660|JQ993649", names(EimSeqs))]
#sp2 <-names(EimSeqs)[grep("JQ993662|JQ993663|KU174485|KU174468|KU174462|JQ993658|KU174479|KU174480|KU174487|KU174459", names(EimSeqs))]

amp_names <- gsub("_ASV.*", "", names18S)
amp_names

## removing lab
LabEim <- EimSeqs[grep("Lab_single", names(EimSeqs))]

names(LabEim) <- paste("wang1141_13_F.Nem_0425_6_3_R", names(LabEim), sep="_")

names(LabEim)

refSeq <- EimSeqs[-grep("Lab_single", names(EimSeqs))]

amp_levels <- levels(as.factor(amp_names))

# removing 28S primer
amp_levels <-amp_levels[!amp_levels=="D3A_5Mod_46_F.D3B_5Mod_46_R"]

# alignming ASV per amplicon to reference
amp_al <- list()
for (i in seq(1:length(amp_levels))){
    amp_al[i] <- AlignSeqs(c(EimASV[grep(amp_levels[i], names(EimASV))], refSeq), anchor=NA, iterations=20, refinements=20, processors=90)
}

### quick fix here to add lab ASVs
#amp_al[9] <- AlignSeqs(c(EimASV[grep(amp_levels[9], names(EimASV))], refSeq, LabEim), anchor=NA, iterations=20, refinements=20, processors=90)
#amp_al[[9]] <- Lab_al

#save fastas
for (i in seq(1:length(amp_levels))){
writeFasta(amp_al[[i]], paste("tmp/amplicon_alignments/amplicon", i, sep=""))
}


#cd tmp
#for f in amplicon_alignments/* ; do ~/iqtree-2.2.0-Linux/bin/iqtree2 -s $f -m TN+F+R2 -B 5000 -T AUTO; done

#~/iqtree-2.2.0-Linux/bin/iqtree2 -s tmp/amplicon_alignments/amplicon9 -m TN+F+R2 -B 5000 -T AUTO -redo



########################################################################
### for 28S ASV alone, they seem aligned, but just in case there are gaps...

ASV28 <- readFasta("tmp/Eimeria28S.fasta")
nameASV28 <- ASV28@id
ASV28 <- sread(ASV28)
names(ASV28) <- nameASV28
ASV28.align <- AlignSeqs(ASV28, anchor=NA, iterations=10, refinements=10, processors=10)

writeFasta(ASV28.align, "tmp/Eimeira_ASV28S_aligned.fasta")


# and this is how I did the tree outside R
#~/iqtree-2.2.0-Linux/bin/iqtree2 -s tmp/Eimeria18S_wild_lab_ref.fasta -m MFP

#~/iqtree-2.2.0-Linux/bin/iqtree2 -s  tmp/Eimeria18S_wild_lab_ref.fasta -m TN+F+R2 -B 5000 -T AUTO -redo

#~/iqtree-2.2.0-Linux/bin/iqtree2 -t tmp/Eimeria18S_wild_lab_ref.fasta.contree -s tmp/Eimeria18S_wild_lab_ref.fasta --scf 100 --prefix siteCF

#~/iqtree-2.2.0-Linux/bin/iqtree2 -s  tmp/Eimeria18S_wild.fasta -m TN+F+R2 -B 5000 -redo

#~/iqtree-2.2.0-Linux/bin/iqtree2 -s  tmp/Eimeria18S_lb.fasta -m TN+F+R2 -B 5000

#~/iqtree-2.2.0-Linux/bin/iqtree2 -s  tmp/Eimeira28S_aligned.fasta -m MFP -B 1000 -T AUTO


### Now we need  to analyse which ASV's are likely from the same Eimeira.
## ASV's within the same sample from different amplicons that correlate well.

####################################################################
############################################################################



# let's do some co-infection analyses. First we need to do some preparations
names_ASVs <- c(names18S, names(LabEim))
amplicon <- data.frame(names_ASVs)
amplicon$species <- "sp"
amplicon$BS <- NA

# sanity check
rownames(Eim@tax_table)==rownames(Eim.T@tax_table)
Eim.T@tax_table[,6] <- names18S
Eim@tax_table[,6] <- names18S

OutGroup <- names(refSeq)[1:3]

species_ASV <- function(fer, name, phyloT, OutGroup) {
    phyloT <- root(phyloT, OutGroup)
    clade.nr <- ggtree::MRCA(phyloT, fer)
    ferT <- extract.clade(phyloT, clade.nr)
    ferrisiASV <- which(names_ASVs%in%ferT$tip.label)
    amplicon$species[ferrisiASV] <- name
    amplicon$species
}

## I want to extract the bootstrap support value in a separate function
BS_ASV <- function(fer, phyloT, OutGroup) {
    phyloT <- root(phyloT, OutGroup)
    clade.nr <- ggtree::MRCA(phyloT, fer)
    ferT <- extract.clade(phyloT, clade.nr)
    ferrisiASV <- which(names_ASVs%in%ferT$tip.label)
    amplicon$BS[ferrisiASV] <- phyloT$node.label[clade.nr-length(phyloT$tip.label)]
    amplicon$BS
}





## import tree
#phyloT <- read.tree("tmp/Eimeria18S_wild_lab_ref.fasta.contree")
phyloT <- list()
for (i in (1:9)){
    phyloT[[i]] <- ggtree::read.tree(paste("tmp/amplicon_alignments/amplicon", i, ".contree", sep=""))
}

for (i in (1:9)){
    amplicon$species <- species_ASV(fer, "ferrisi", phyloT[[i]], OutGroup)
    amplicon$species <- species_ASV(ver, "vermiformis", phyloT[[i]], OutGroup)
    amplicon$species <- species_ASV(fal, "falciformis", phyloT[[i]], OutGroup)
    amplicon$BS <- BS_ASV(fer, phyloT[[i]], OutGroup)
    amplicon$BS <- BS_ASV(ver, phyloT[[i]], OutGroup)
    amplicon$BS <- BS_ASV(fal, phyloT[[i]], OutGroup)
}

amplicon$species[grep("D3A_5", amplicon$names_ASVs)] <- "28S"
amplicon$BS[grep("D3A_5", amplicon$names_ASVs)] <- 1

# inspection
amplicon$names_ASVs[amplicon$spec=="ferrisi"]
amplicon$names_ASVs[amplicon$spec=="falciformis"]

amplicon$names_ASVs[amplicon$spec=="sp"]

amplicon$names_ASVs[amplicon$spec=="vermiformis"]

#### manual visualization from the tree for the 2 unassigned  asv'S

amplicon$species[amplicon$names_ASVs=="MarkN_10_F.Proti440R_28_R_ASV_5"] <- "falciformis"
amplicon$BS[amplicon$names_ASVs=="MarkN_10_F.Proti440R_28_R_ASV_5"] <- 96

amplicon$BS[amplicon$names_ASVs=="18S_0067a_deg_3Mod_53_F.NSR399_3Mod_53_R_ASV_4"] <- 59

amplicon$species[amplicon$names_ASVs=="wang1141_13_F.Nem_0425_6_3_R_ASV_2"] <- "falciformis"
amplicon$BS[amplicon$names_ASVs=="wang1141_13_F.Nem_0425_6_3_R_ASV_2"] <- 43

amplicon$species[amplicon$names_ASVs=="Prot1702_32_F.wang1624CR6S_16_R_ASV_2"] <- "falciformis"
amplicon$BS[amplicon$names_ASVs=="Prot1702_32_F.wang1624CR6S_16_R_ASV_2"] <- 62

## now we save this in a new table and remove the lab ASV's
amplicon2 <- amplicon

amplicon <- amplicon[-c(36, 37, 38, 39),]

#sanity check
Eim@tax_table[,7] <- amplicon$spec
Eim.TSS@tax_table[,7] <- amplicon$spec
Eim.T@tax_table[,7] <- amplicon$spec

get_taxa_unique(Eim, "Genus")
################################################################
######### Preparing phyloseq objects for plotting and so on
amp_names <- gsub("_ASV.*", "", names18S)
Eim@tax_table[,6] <- amp_names
Eim.TSS@tax_table[,6] <- amp_names
Eim.T@tax_table[,6] <- amp_names

Eim.m <- Eim
Eim.TSS.m <- Eim.TSS
Eim.T.m <- Eim.T

#Eim.m@tax_table[,6] <- amplicon$names18S
#Eim.TSS.m@tax_table[,6] <- amplicon$names18S

# removing empty samples
Eim <- phyloseq::prune_samples(sample_sums(Eim)>0, Eim)
Eim.TSS <- phyloseq::prune_samples(sample_sums(Eim.TSS)>0, Eim.TSS)
Eim.T <- phyloseq::prune_samples(sample_sums(Eim.T)>0, Eim.T)


Eim18 <- Eim
Eim.TSS18 <- Eim.TSS
Eim.T18 <- Eim.T
Eim28 <- Eim
Eim.TSS28 <- Eim.TSS
Eim.T28 <- Eim.T
Eim18 <- subset_taxa(Eim18, !Species=="28S")
Eim.TSS18 <- subset_taxa(Eim.TSS18, !Species=="28S")
Eim.T18 <- subset_taxa(Eim.T18, !Species=="28S")
Eim28 <- subset_taxa(Eim28, Species=="28S")
Eim.TSS28 <- subset_taxa(Eim.TSS28, Species=="28S")
Eim.T28 <- subset_taxa(Eim.T28, Species=="28S")

#Eim28@tax_table[10,6]

get_taxa_unique(Eim18, "Species")

# now this is for plotting co-infections
Eim18 <- phyloseq::prune_samples(sample_sums(Eim18)>0, Eim18)
Eim.TSS18 <- phyloseq::prune_samples(sample_sums(Eim.TSS18)>0, Eim.TSS18)
Eim.T18 <- phyloseq::prune_samples(sample_sums(Eim.T18)>0, Eim.T18)
#eim18.TSS <- transform_sample_counts(Eim18, function(x) x / sum(x)) 
Eim18_sp <- Eim.T18

Eim18_sp@tax_table[,6] <- Eim18_sp@tax_table[,5]

Eim18_sp <- tax_glom(Eim18_sp, taxrank="Species")

Eim18_sp

#sanity check
colnames(Eim18_sp@otu_table)==rownames(Eim18_sp@tax_table)
#colnames(Eim18_sp@otu_table) <- Eim18_sp@tax_table[,5]


#########################################################
# correlation matrix
library(Hmisc)
library(Matrix)
library(SpiecEasi)
library(igraph)

eim <- (Eim.T@otu_table)
tax <- data.frame(Eim.T@tax_table) 

head(eim)

head(tax)

taxa_sums(Eim.T)

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
p.yes.r.str <- abs(p.yes.r)>0.5 # output is logical vector

p.yes.rr <- p.yes.r.str*r.val # use logical vector for subscripting.

p.yes.rr

adjm <- as.matrix(p.yes.r)

adjm

colnames(adjm) <- names18S
rownames(adjm) <- names18S

net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE) 

V(net.grph)


summary(adjm<0)

E(net.grph)$weight

E(net.grph)$color <- "dodgerblue4"

E(net.grph)$color[which(E(net.grph)$weight<0)] <- "#FF00FF"

E(net.grph)$color


# we also want the node color to code for amplicon
amp <- as.factor(gsub("_ASV_[0-9]", "", names18S))
nb.col <- length(levels(amp))
coul <- colorRampPalette(brewer.pal(8, "Accent"))(nb.col)
mc <- coul[as.numeric(amp)]
mc2 <- mc[-28]

# sanity check
names(V(net.grph))==amplicon$names_ASVs

amplicon$BS <- as.numeric(amplicon$BS)

pdf("fig/Eimeria_ASVs_Network.pdf",
                width =15, height = 15)
set.seed(1235)
plot(net.grph,
     vertex.label=amplicon$species,
#     edge.width=as.integer(cut(E(net.grph)$weight, breaks=6))/2,
     edge.width=abs(E(net.grph)$weight)*1.5,
     vertex.color=adjustcolor(mc, 0.8),
     vertex.size=amplicon$BS/10,
     frame.col="grey")
dev.off()

modules =cluster_fast_greedy(net.grph, weights=abs(E(net.grph)$weight))
modularity(modules)
sizes(modules)

#V(net.grph2)$cluster <- modules$membership
plot_dendrogram(modules)

cluster_ID <- modules$membership



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
