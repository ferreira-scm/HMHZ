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


## We make a tree per amplicon
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

##### there will be some sections here commented out because they take time to run and they don't need to be ran all the time.
## align
#allEim <- OrientNucleotides(allEim)
#set.seed(12345)
#Eim.align <- AlignSeqs(allEim, anchor=NA, iterations=20, refinements=20, processors=90)
#Eim.align <- AdjustAlignment(Eim.align)

#writeFasta(Eim.align, "tmp/Eimeria18S_wild_lab_ref.fasta")
#writeFasta(allEim, "tmp/Eimeria18S_wild_lab_ref_notAligned.fasta")

################ getting amplicon names
amp_names <- gsub("_ASV.*", "", names18S)
amp_names

# selecting lab asv's
LabEim <- EimSeqs[grep("Lab_single", names(EimSeqs))]
names(LabEim) <- paste("wang1141_13_F.Nem_0425_6_3_R", names(LabEim), sep="_")
## removing lab asv's
refSeq <- EimSeqs[-grep("Lab_single", names(EimSeqs))]
amp_levels <- levels(as.factor(amp_names))

# removing 28S primer
amp_levels <-amp_levels[!amp_levels=="D3A_5Mod_46_F.D3B_5Mod_46_R"]

# alignming ASV per amplicon to reference
#amp_al <- list()
#for (i in seq(1:length(amp_levels))){
#    amp_al[i] <- AlignSeqs(c(EimASV[grep(amp_levels[i], names(EimASV))], refSeq), anchor=NA, iterations=20, refinements=20, processors=90)
#}

### quick fix here to add lab ASVs
#amp_al[9] <- AlignSeqs(c(EimASV[grep(amp_levels[9], names(EimASV))], refSeq, LabEim), anchor=NA, iterations=20, refinements=20, processors=90)
#amp_al[[9]] <- Lab_al

#save fastas
#for (i in seq(1:length(amp_levels))){
#writeFasta(amp_al[[i]], paste("tmp/amplicon_alignments/amplicon", i, sep=""))
#}


############## run this outside R
#cd tmp
#for f in amplicon_alignments/* ; do ~/iqtree-2.2.0-Linux/bin/iqtree2 -s $f -m TN+F+R2 -B 5000 -T AUTO; done

#~/iqtree-2.2.0-Linux/bin/iqtree2 -s tmp/amplicon_alignments/amplicon9 -m TN+F+R2 -B 5000 -T AUTO -redo

########################################################################

### for 28S ASV alone, they seem aligned, but just in case there are gaps...
#ASV28 <- readFasta("tmp/Eimeria28S.fasta")
#nameASV28 <- ASV28@id
#ASV28 <- sread(ASV28)
#names(ASV28) <- nameASV28
#ASV28.align <- AlignSeqs(ASV28, anchor=NA, iterations=10, refinements=10, processors=10)
#writeFasta(ASV28.align, "tmp/Eimeira_ASV28S_aligned.fasta")

# trees outside R
#~/iqtree-2.2.0-Linux/bin/iqtree2 -s tmp/Eimeria18S_wild_lab_ref.fasta -m MFP

#~/iqtree-2.2.0-Linux/bin/iqtree2 -s  tmp/Eimeria18S_wild_lab_ref.fasta -m TN+F+R2 -B 5000 -T AUTO -redo

#~/iqtree-2.2.0-Linux/bin/iqtree2 -t tmp/Eimeria18S_wild_lab_ref.fasta.contree -s tmp/Eimeria18S_wild_lab_ref.fasta --scf 100 --prefix siteCF

#~/iqtree-2.2.0-Linux/bin/iqtree2 -s  tmp/Eimeira28S_aligned.fasta -m MFP -B 1000 -T AUTO


### Now we need  to analyse which ASV's are likely from the same Eimeira.
## ASV's within the same sample from different amplicons that correlate well.
####################################################################
############################################################################
# I need set 3 groups of reference sequences (ferrisi, vermiformis, falciformis)
fer <- c(names(EimSeqs)[grep("ferrisi",names(EimSeqs))],"AF246717_Eimeria_telekii", "AF311643_Eimeria_separata")
fal <- names(EimSeqs)[grep("falciformis|JQ993650|KU174466|KU174449|KU174463|KU174481|JQ993657|JQ993666|KU174452|KU174470|KU174472|KU174473|KU174471|KU174460|KU174450|KU174455|KU174451|KU174456|KU174474|KU174461|KU174454|KU174458|KU174464|AF311644|KU174457|KU174453|KU174469|KU174476|KU174486|KU174478|KU174484|KU192956|JQ993665",names(EimSeqs))]
ver <- names(EimSeqs)[grep("vermiformis|KU192931|KU192936|KU174467|KU192916|KU174465",names(EimSeqs))]
#sp1 <-names(EimSeqs)[grep("U40263|KU174483|KU192965|KU192958|KU192961|JQ993660|JQ993649", names(EimSeqs))]
#sp2 <-names(EimSeqs)[grep("JQ993662|JQ993663|KU174485|KU174468|KU174462|JQ993658|KU174479|KU174480|KU174487|KU174459", names(EimSeqs))]

# let's do some co-infection analyses. First we need to do some preparations


names_ASVs <- c(names18S, names(LabEim))
amplicon <- data.frame(names_ASVs)
amplicon$species <- "sp"
amplicon$BS <- NA

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

## import trees
#phyloT <- read.tree("tmp/Eimeria18S_wild_lab_ref.fasta.contree")
phyloT <- list()
for (i in (1:9)){
    phyloT[[i]] <- ggtree::read.tree(paste("tmp/amplicon_alignments/amplicon", i, ".contree", sep=""))
}

# now assign species based on trees
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

#### manual visualization from the tree for the 4 unassigned  asv'S
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

get_taxa_unique(Eim, "Species")

#########################################################
# correlation matrix and network
library(Hmisc)
library(Matrix)
library(igraph)

# removing empty samples
Eim <- phyloseq::prune_samples(sample_sums(Eim)>0, Eim)
Eim.TSS <- phyloseq::prune_samples(sample_sums(Eim.TSS)>0, Eim.TSS)
Eim.T <- phyloseq::prune_samples(sample_sums(Eim.T)>0, Eim.T)

eim <- (Eim.T@otu_table)
tax <- data.frame(Eim.T@tax_table) 

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
r.val = otu.cor$r # select all the correlation values
p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion
## select asv based on rho
p.yes.r.str <- abs(p.yes.r)>0.5 # output is logical vector
p.yes.rr <- p.yes.r.str*r.val # use logical vector for subscripting.
adjm <- as.matrix(p.yes.r)
colnames(adjm) <- names18S
rownames(adjm) <- names18S
net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE) 

V(net.grph)
### negative correlations
summary(adjm<0)
### colour negative edges
E(net.grph)$weight
E(net.grph)$color <- "dodgerblue4"
E(net.grph)$color[which(E(net.grph)$weight<0)] <- "#FF00FF"
E(net.grph)$color

# we also want the node color to code for amplicon
amp <- as.factor(gsub("_ASV_[0-9]", "", names18S))
nb.col <- length(levels(amp))
coul <- colorRampPalette(brewer.pal(8, "Accent"))(nb.col)
mc <- coul[as.numeric(amp)]

# sanity check
names(V(net.grph))==amplicon$names_ASVs
amplicon$BS <- as.numeric(amplicon$BS)

# now plotting
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

# the modules now
modules =cluster_fast_greedy(net.grph, weights=abs(E(net.grph)$weight))
modularity(modules)
sizes(modules)

#V(net.grph2)$cluster <- modules$membership
#plot_dendrogram(modules)

amplicon$species_FINAL[modules$membership==1] <- "ferrisi"
amplicon$species_FINAL[modules$membership==2] <- "falciformis"
amplicon$species_FINAL[amplicon$species=="vermiformis"] <- "vermiformis"

# assigning species
Eim@tax_table[,7] <- amplicon$species_FINAL
Eim.T@tax_table[,7] <- amplicon$species_FINAL
Eim.TSS@tax_table[,7] <- amplicon$species_FINAL

############################ we need to normalise between amplicons, before we merge
# z-score normalization

# substract the mean Eimeria (species level) from the Eimeria in each sample divided by the standard deviation

# first, divide 18S from 28S

#Eim.T18

#Eim.T28

# then agglomerate species

#Eim.T18@tax_table[,6] <- "Eimeria"

#Eim.T28@tax_table[,6] <- "Eimeria"

#Eim.T18_sp <- tax_glom(Eim.T18, "Species")
#Eim.T28_sp <- tax_glom(Eim.T28, "Species")

# do normalization

#mean18 <- mean(sample_sums(Eim.T18_sp))
#sd18 <- sd(sample_sums(Eim.T18_sp))
#Z.Eim18 <- Eim.T18_sp
#Z.Eim18@otu_table <- ((Eim.T18_sp@otu_table-mean18)/sd18)+

#mean28 <- mean(sample_sums(Eim.T28_sp))
#sd28 <- sd(sample_sums(Eim.T28_sp))
#Z.Eim28 <- Eim.T28_sp
#Z.Eim28@otu_table <- ((Eim.T28_sp@otu_table-mean28)/sd28)

#Z.Eim28@otu_table


# join

##### Do this but for the sum of all ASVs per amplicon


#Eim.T@tax_table[,6] <- "Eimeria"

#Eim.T_sp <- tax_glom(Eim.T, "Species")


#taxa_sums(Eim.T_sp)

#Eim_s <- tax_glom(Eim, "Species")

#Eim_s

#fPS

#161/619*100
