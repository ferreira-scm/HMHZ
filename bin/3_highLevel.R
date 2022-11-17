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

source("bin/Pre_Functions.R")

##########################################
### composition plots
tPS <- subset_taxa(fPS, !phylum%in%NA)
tPSTSS <- transform_sample_counts(tPS, function(x) x / sum(x)) 
tps <- psmelt(tPSTSS)

mycolor=colorRampPalette(brewer.pal(8, "Paired"))(22)

comp <- ggplot(tps, aes(x=Sample, y=Abundance, fill=fct_reorder(phylum, Abundance)))+
    geom_bar(position="stack", stat="identity")+
#    scale_x_discrete(labels=tps$EH_ID, breaks=tps$Sample)+
    scale_fill_manual("genus", values=mycolor)+
#    facet_grid(~dpi, scales="free")+
    labs(fill="Phylum", y="Proportion within total sequence")+
    theme_bw(base_size=14)+
    guides(fill=guide_legend(ncol=4))+
    theme(axis.text.y = element_text(colour = 'black', size = 14),
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         legend.key = element_blank(),
#        strip.background = element_rect(colour="black", fill="white"),
         legend.text = element_text(colour = 'black', size = 12),
         legend.position="bottom") 

## making other tables for plotting
tPS@sam_data$one <- "one"
ave.ps <- merge_samples(tPS, "one")
ave.tss <- transform_sample_counts(ave.ps, function(x) x / sum(x)) 
ave.tss <- psmelt(ave.tss)

comp.ave <- ggplot(ave.tss, aes(x=Sample, y=Abundance, fill=phylum))+
    geom_bar(position="stack", stat="identity")+
#    scale_x_discrete(labels=tps$EH_ID, breaks=tps$Sample)+
    scale_fill_manual("genus", values=mycolor)+
#    facet_grid(~dpi, scales="free")+
    labs(y="Proportion within total sequence")+
    theme_bw(base_size=14)+
    theme(axis.text.y = element_text(colour = 'black', size = 14),
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         legend.key = element_blank(),
         strip.background = element_rect(colour="black", fill="white"),
#        legend.text = element_text(colour = 'black', size = 12, face = 'italic'),
         legend.position="none") 

com_g <- plot_grid(comp.ave, comp, ncol=2, align="h", axis="b", rel_widths=c(1, 5))
ggsave("fig/Composition_barplot.pdf", com_g, height=7, width=10, dpi=400)
ggsave("fig/Composition_barplot.png", com_g, height=7, width=10, dpi=400)

#####################################################################
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

Eim.fg <- readFasta("../LabMicrobiome/tmp/Eimeria_falciformis_NODE_3780_18S.fa")
Eimfg <- sread(Eim.fg)
names(Eimfg)  <- "Eimeira_falciformis_NODE3780"
Eimfg

Eim.vg <- readFasta("../LabMicrobiome/tmp/Eimeria_vermiformis_NODE_17893_18S_plus.fa")
Eimvg <- sread(Eim.vg)
names(Eimvg)  <- "Eimeira_vermiformis_NODE17893"
Eimvg

FinalEimeria <- c(Eimfg, Eimvg, EimSeqs, EimASV)

#### I want to align only with wild Eimeiras
refEim <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_ref.fa")  
names(refEim) <- gsub("(\\s)", "_", names(refEim))
Wild.Eim <- c(EimASV, refEim)
# ops forgout the outgroup
T.out <- c("JQ993669", "JQ993670", "JQ993671")
T.out <- read.GenBank(T.out)
T.out.ID <- paste(names(T.out), attr(T.out, "species"), sep="_") 
T.O <- T.out %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
names(T.O) <- T.out.ID
# and I want the lab EImerias too
LabASV <- readFasta("../LabMicrobiome/tmp/Eimeria_lab_ASV.fa")
nmLab <- LabASV@id
nmLab
LabASV <- sread(LabASV)
names(LabASV) <- nmLab
#Lab <- LabASV@sread
WildEim <- c(Wild.Eim, T.O, LabASV)
#Wild.align <- AlignSeqs(WildEim, anchor=NA, iterations=20, refinements=20, processors=90)
#Wild.align <- AdjustAlignment(Wild.align)
names(allEim[c("ASV1", "ASV2", "ASV3", "ASV4", "ASV5")]) <- nmLab
########################################################################

library(DECIPHER)
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








#########################################################
### calculate dissimilarity distancelibrary(vegan)
library(ape)

Eim18.TSS <- eim18.TSS
# sanity check
colnames(eim18.TSS@otu_table)==rownames(Eim18.TSS@tax_table)

plot_ordination(eim18.TSS, res18)

asvN <- seq(1,47,1)
asvN <- paste("ASV", asvN, sep="")
colnames(Eim18.TSS@otu_table) <- paste(asvN, Eim18.TSS@tax_table[,5], sep="_")

colnames(Eim18.TSS@otu_table)

Eim18_sp@otu_table

biplot(res18.2, Eim18_sp@otu_table)

biplot(res18, Eim18.TSS@otu_table)

otu.st <- apply(Eim18.TSS@otu_table, 2, scale, center=TRUE, scale=TRUE)
otu.st2 <- apply(Eim18_sp@otu_table, 2, scale, center=TRUE, scale=TRUE)

biplot(res18.2, otu.st2,
       xlabs=rep("o", 362)
       )

biplot(res18, otu.st,
       xlabs=rep("o", 362),
       cex=c(0.7, 0.8))


#phyloseq::plot_heatmap(fPS)


plot_net(Eim18, "bray")

phyloseq::plot_richness(Eim18)



pcoa_bc <- ordinate(eim18.TSS, "PCoA", "bray")

pcoa.p <- plot_ordination(eim18.TSS, pcoa_bc)

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
