setwd("/SAN/Susanas_den/gitProj/HMHZ/")

library(ggplot2)
library(reshape)
library(phyloseq, lib="/usr/local/lib/R/site-library/")
library(data.table)
library(microbiome)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
library(vegan)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(ggeffects, lib="/usr/local/lib/R/site-library/")

PS=readRDS(file="tmp/PSpre.R")

pPS=readRDS(file="tmp/pPSpre.R")

pPS

# let's include host and ecological variables
pPS <- subset_samples(pPS, !Sex=="<NA>")
sample_data(pPS)$Sex

pPS <- subset_samples(pPS, !BMI=="NA")

sample_data(pPS)$BMI

pPS

#PS <- subset_samples(PS, !Sex=="<NA>")
#PS <- subset_samples(PS, !BMI=="NA")


### beta diversity (Bray-Curtis distances), PERMANOVA and visualization with NMDS
# we will do it 3 times:
# 1) No transformation, and accounting for library size (position 1) on the PERMANOVA
#### permanovas on beta diversisty

PCdis=phyloseq::distance((pPS), method="bray", type="samples")

sdata=sample_data(pPS)

permaPS=adonis2(PCdis~
            sdata$eimASV+
            sdata$oxyASV +
            sdata$triASV+
            sdata$hexASV +
            sdata$ascASV +
            sdata$criASV +
            sdata$hekASV+
            sdata$hymASV,
        permutations = 10000, method = "bray")

permaPS

# 2) Using variance stabilized transformation
library(DESeq2)
pPSvs = pPS
otu_table(pPSvs) <- otu_table(pPSvs)+1
pPSds = phyloseq_to_deseq2(pPSvs, ~1)
pPSds = estimateSizeFactors(pPSds)
pPSds = estimateDispersions(pPSds, fitType = "local")
pPSvst = pPS

otu_table(pPSvst) <- otu_table(getVarianceStabilizedData(pPSds), taxa_are_rows = TRUE)
# replace negatives with zero
otu_table(pPSvst)[otu_table(pPSvst)<0] <- 0

sample_sums(pPSvst)

# adding parasite info
vstPSeim <- subset_taxa(pPSvst, family%in%"Eimeriidae")
vstPScri <- subset_taxa(pPSvst, family%in%"Cryptosporidiidae")
#PSsar <- subset_taxa(PS, family%in%"Sarcocystidae")
# for nematoda
vstPSoxy <- subset_taxa(pPSvst, family%in%"Oxyuridae")
vstPStri <- subset_taxa(pPSvst, family%in%"Trichuridae")
vstPShex <- subset_taxa(pPSvst, family%in%"Heteroxynematidae")
#PSstr <- subset_taxa(PS, family%in%"Strongyloididae")
vstPSasc <- subset_taxa(pPSvst, family%in%"Ascaridiidae")
vstPSspi <- subset_taxa(pPSvst, family%in%"Spirocercidae")
vstPShek <- subset_taxa(pPSvst, family%in%"Heterakidae")
# for platyhelminth
vstPShym <- subset_taxa(pPSvst, family%in%"Hymenolepididae")

sample_data(pPSvst)$eimASV <- sample_sums(vstPSeim)
sample_data(pPSvst)$criASV <- sample_sums(vstPScri)
sample_data(pPSvst)$oxyASV <- sample_sums(vstPSoxy)
sample_data(pPSvst)$triASV <- sample_sums(vstPStri)
sample_data(pPSvst)$hexASV <- sample_sums(vstPShex)
sample_data(pPSvst)$ascASV <- sample_sums(vstPSasc)
sample_data(pPSvst)$hekASV <- sample_sums(vstPShek)
sample_data(pPSvst)$hymASV <- sample_sums(vstPShym)
sample_data(pPSvst)$spiASV <- sample_sums(vstPSspi)


pPSvst <- subset_taxa(pPSvst, !family %in% "Eimeriidae")
pPSvst <- subset_taxa(pPSvst, !family %in% "Cryptosporidiidae")
pPSvst <- subset_taxa(pPSvst, !family %in% "Oxyuridae")
pPSvst <- subset_taxa(pPSvst, !family %in% "Trichuridae")
pPSvst <- subset_taxa(pPSvst, !family %in% "Heteroxynematidae")
pPSvst <- subset_taxa(pPSvst, !family %in% "Ascaridiidae")
pPSvst <- subset_taxa(pPSvst, !family %in% "Heterakidae") 
pPSvst <- subset_taxa(pPSvst, !family %in% "Hymenolepididae")
pPSvst <- subset_taxa(pPSvst, !family %in% "Spirocercidae")

# does it make a difference when I normalize after removing parasite counts?
ppPS <- pPS
ppPS <- subset_taxa(ppPS, !family %in% "Eimeriidae")
ppPS <- subset_taxa(ppPS, !family %in% "Cryptosporidiidae")
ppPS <- subset_taxa(ppPS, !family %in% "Oxyuridae")
ppPS <- subset_taxa(ppPS, !family %in% "Trichuridae")
ppPS <- subset_taxa(ppPS, !family %in% "Heteroxynematidae")
ppPS <- subset_taxa(ppPS, !family %in% "Ascaridiidae")
ppPS <- subset_taxa(ppPS, !family %in% "Heterakidae") 
ppPS <- subset_taxa(ppPS, !family %in% "Hymenolepididae")
ppPS <- subset_taxa(ppPS, !family %in% "Spirocercidae")

ppPSvs = ppPS

otu_table(ppPSvs) <- otu_table(ppPSvs)+1
ppPSds = phyloseq_to_deseq2(ppPSvs, ~1)
ppPSds = estimateSizeFactors(ppPSds)
ppPSds = estimateDispersions(ppPSds, fitType = "local")

ppPSvst = ppPS

otu_table(ppPSvst) <- otu_table(getVarianceStabilizedData(ppPSds), taxa_are_rows = TRUE)
# replace negatives with zero
otu_table(ppPSvst)[otu_table(ppPSvst)<0] <- 0

sample_data(ppPSvst)$eimASV <- sample_data(pPSvst)$eimASV
sample_data(ppPSvst)$criASV <- sample_data(pPSvst)$criASV
sample_data(ppPSvst)$oxyASV <- sample_data(pPSvst)$oxyASV
sample_data(ppPSvst)$triASV <- sample_data(pPSvst)$triASV
sample_data(ppPSvst)$hexASV <- sample_data(pPSvst)$hexASV
sample_data(ppPSvst)$ascASV <- sample_data(pPSvst)$ascASV
sample_data(ppPSvst)$hekASV <- sample_data(pPSvst)$hekASV
sample_data(ppPSvst)$hymASV <- sample_data(pPSvst)$hymASV
sample_data(ppPSvst)$spiASV <- sample_data(pPSvst)$spiASV

#### finaly we are ready to do the permanovas
vstdis=phyloseq::distance((pPSvst), method="bray", type="samples")

sdata=sample_data(pPSvst)

permaPS=adonis2(vstdis~
            sdata$Year+
            sdata$Transect+
            sdata$Sex+
            sdata$BMI+
            sdata$HI+
            sdata$eimASV+
            sdata$oxyASV +
            sdata$triASV+
            sdata$hexASV +
            sdata$ascASV +
            sdata$criASV +
            sdata$hekASV+
            sdata$spiASV+
            sdata$hymASV,
        permutations = 10000, method = "bray")
permaPS
write.table(permaPS, "tmp/permapPS_vst_prem.txt")

vvstdis=phyloseq::distance((ppPSvst), method="bray", type="samples")
vvstjac <- phyloseq::distance(ppPSvst, method="jaccard")

sdata=sample_data(ppPSvst)

ppermaPS=adonis2(vvstdis~
            sdata$Year+
            sdata$Transect+
            sdata$Sex+
            sdata$BMI+
            sdata$HI+
            sdata$eimASV+
            sdata$oxyASV +
            sdata$triASV+
            sdata$hexASV +
            sdata$ascASV +
            sdata$criASV +
            sdata$hekASV+
            sdata$spiASV+
            sdata$hymASV,
        permutations = 10000, method = "bray")

ppermaPS

write.table(ppermaPS, "tmp/permapPS_prem_vst.txt")

jacpermaPS=adonis2(vstjac~
            sdata$Year+
            sdata$Transect+
            sdata$Sex+
            sdata$BMI+
            sdata$HI+
            sdata$eimASV+
            sdata$oxyASV +
            sdata$triASV+
            sdata$hexASV +
            sdata$ascASV +
            sdata$criASV +
            sdata$hekASV+
            sdata$spiASV+
            sdata$hymASV,
        permutations = 10000, method = "jaccard")

jacpermaPS


#5) on seperate domains
### Average relative OTU abundances per domain
euPS <- subset_taxa(pPS, superkingdom%in%"Eukaryota")
baPS<- subset_taxa(pPS, superkingdom%in%"Bacteria")

euvs = euPS
otu_table(euvs) <- otu_table(euvs)+1
euds = phyloseq_to_deseq2(euvs, ~1)
euds = estimateSizeFactors(euds)
euds = estimateDispersions(euds, fitType = "local")
euvst = euPS

otu_table(euvst) <- otu_table(getVarianceStabilizedData(euds), taxa_are_rows = TRUE)
# replace negatives with zero
otu_table(euvst)[otu_table(euvst)<0] <- 0

bavs = baPS
otu_table(bavs) <- otu_table(bavs)+1
bads = phyloseq_to_deseq2(bavs, ~1)
bads = estimateSizeFactors(bads)
bads = estimateDispersions(bads, fitType = "local")
bavst = baPS

otu_table(bavst) <- otu_table(getVarianceStabilizedData(bads), taxa_are_rows = TRUE)
# replace negatives with zero
otu_table(bavst)[otu_table(bavst)<0] <- 0




# adding parasite info
vstPSeim <- subset_taxa(euvst, family%in%"Eimeriidae")
vstPScri <- subset_taxa(euvst, family%in%"Cryptosporidiidae")
#PSsar <- subset_taxa(PS, family%in%"Sarcocystidae")
# for nematoda
vstPSoxy <- subset_taxa(euvst, family%in%"Oxyuridae")
vstPStri <- subset_taxa(euvst, family%in%"Trichuridae")
vstPShex <- subset_taxa(euvst, family%in%"Heteroxynematidae")
#PSstr <- subset_taxa(PS, family%in%"Strongyloididae")
vstPSasc <- subset_taxa(euvst, family%in%"Ascaridiidae")
vstPSspi <- subset_taxa(euvst, family%in%"Spirocercidae")
vstPShek <- subset_taxa(euvst, family%in%"Heterakidae")
# for platyhelminth
vstPShym <- subset_taxa(euvst, family%in%"Hymenolepididae")

sample_data(euvst)$eimASV <- sample_sums(vstPSeim)
sample_data(euvst)$criASV <- sample_sums(vstPScri)
sample_data(euvst)$oxyASV <- sample_sums(vstPSoxy)
sample_data(euvst)$triASV <- sample_sums(vstPStri)
sample_data(euvst)$hexASV <- sample_sums(vstPShex)
sample_data(euvst)$ascASV <- sample_sums(vstPSasc)
sample_data(euvst)$hekASV <- sample_sums(vstPShek)
sample_data(euvst)$hymASV <- sample_sums(vstPShym)
sample_data(euvst)$spiASV <- sample_sums(vstPSspi)

sample_data(bavst)$eimASV <- sample_sums(vstPSeim)
sample_data(bavst)$criASV <- sample_sums(vstPScri)
sample_data(bavst)$oxyASV <- sample_sums(vstPSoxy)
sample_data(bavst)$triASV <- sample_sums(vstPStri)
sample_data(bavst)$hexASV <- sample_sums(vstPShex)
sample_data(bavst)$ascASV <- sample_sums(vstPSasc)
sample_data(bavst)$hekASV <- sample_sums(vstPShek)
sample_data(bavst)$hymASV <- sample_sums(vstPShym)
sample_data(bavst)$spiASV <- sample_sums(vstPSspi)


euvst <- subset_taxa(euvst, !family %in% "Eimeriidae")
euvst <- subset_taxa(euvst, !family %in% "Cryptosporidiidae")
euvst <- subset_taxa(euvst, !family %in% "Oxyuridae")
euvst <- subset_taxa(euvst, !family %in% "Trichuridae")
euvst <- subset_taxa(euvst, !family %in% "Heteroxynematidae")
euvst <- subset_taxa(euvst, !family %in% "Ascaridiidae")
euvst <- subset_taxa(euvst, !family %in% "Heterakidae") 
euvst <- subset_taxa(euvst, !family %in% "Hymenolepididae")
euvst <- subset_taxa(euvst, !family %in% "Spirocercidae")


# permanova
evstdis=phyloseq::distance((euvst), method="bray", type="samples")
sdata=sample_data(euvst)

epermaPS=adonis2(evstdis~
            sdata$Year+
            sdata$Transect+
            sdata$Sex+
            sdata$BMI+
            sdata$HI+
            sdata$eimASV+
            sdata$oxyASV +
            sdata$triASV+
            sdata$hexASV +
            sdata$ascASV +
            sdata$criASV +
            sdata$hekASV +
            sdata$spi +
            sdata$hymASV,
        permutations = 10000, method = "bray")

epermaPS

write.table(epermaPS, "tmp/euka_permapPS.txt")

bvstdis=phyloseq::distance((bavst), method="bray", type="samples")
sdata=sample_data(bavst)

bpermaPS=adonis2(bvstdis~
            sdata$Year+
            sdata$Transect+
            sdata$Sex+
            sdata$BMI+
            sdata$HI+
            sdata$eimASV+
            sdata$oxyASV +
            sdata$triASV+
            sdata$hexASV +
            sdata$ascASV +
            sdata$criASV +
            sdata$hekASV +
            sdata$spi +
            sdata$hymASV,
        permutations = 10000, method = "bray")

bpermaPS
write.table(bpermaPS, "tmp/bac_permapPS.txt")

plot_ordination(pPSvst, ordinate(pPSvst, "NMDS", "bray"), color = "Year")

# 3) relative abundance
PSTSS = transform_sample_counts(pPS, function(x) x / sum(x) )
PSeim <- subset_taxa(PSTSS, family%in%"Eimeriidae")
PScri <- subset_taxa(PSTSS, family%in%"Cryptosporidiidae")
#PSsar <- subset_taxa(PS, family%in%"Sarcocystidae")
# for nematoda
PSoxy <- subset_taxa(PSTSS, family%in%"Oxyuridae")
PStri <- subset_taxa(PSTSS, family%in%"Trichuridae")
PShex <- subset_taxa(PSTSS, family%in%"Heteroxynematidae")
#PSstr <- subset_taxa(PS, family%in%"Strongyloididae")
PSasc <- subset_taxa(PSTSS, family%in%"Ascaridiidae")
PSspi <- subset_taxa(PSTSS, family%in%"Spirocercidae")
PShek <- subset_taxa(PSTSS, family%in%"Heterakidae")
# for platyhelminth
PShym <- subset_taxa(PSTSS, family%in%"Hymenolepididae")

sample_data(PSTSS)$eimASV <- sample_sums(PSeim)
sample_data(PSTSS)$criASV <- sample_sums(PScri)
sample_data(PSTSS)$oxyASV <- sample_sums(PSoxy)
sample_data(PSTSS)$triASV <- sample_sums(PStri)
sample_data(PSTSS)$hexASV <- sample_sums(PShex)
sample_data(PSTSS)$ascASV <- sample_sums(PSasc)
sample_data(PSTSS)$hekASV <- sample_sums(PShek)
sample_data(PSTSS)$hymASV <- sample_sums(PShym)
sample_data(PSTSS)$spiASV <- sample_sums(PSspi)

PSTSS <- subset_taxa(PSTSS, !family %in% "Eimeriidae")
PSTSS <- subset_taxa(PSTSS, !family %in% "Cryptosporidiidae")
PSTSS <- subset_taxa(PSTSS, !family %in% "Oxyuridae")
PSTSS <- subset_taxa(PSTSS, !family %in% "Trichuridae")
PSTSS <- subset_taxa(PSTSS, !family %in% "Heteroxynematidae")
PSTSS <- subset_taxa(PSTSS, !family %in% "Ascaridiidae")
PSTSS <- subset_taxa(PSTSS, !family %in% "Heterakidae") 
PSTSS <- subset_taxa(PSTSS, !family %in% "Hymenolepididae")
PSTSS <- subset_taxa(PSTSS, !family %in% "Spirocercidae")


TSSdis=phyloseq::distance((PSTSS), method="bray", type="samples")

sdata=sample_data(PSTSS)

tpermaPS=adonis2(TSSdis~
            sdata$Year+
            sdata$Transect+
            sdata$Sex+
            sdata$BMI+
            sdata$HI+
            sdata$eimASV+
            sdata$oxyASV +
            sdata$triASV+
            sdata$hexASV +
            sdata$ascASV +
            sdata$criASV +
            sdata$hekASV+
            sdata$spiASV+
            sdata$hymASV,
        permutations = 10000, method = "bray")

tpermaPS

# 4) rarefy
library("metagMisc")
library("iNEXT")

#ipPS <- prepare_inext(pPS)
#a <- iNEXT(ipPS, q=0)
#b <- ggiNEXT(a, type=1)
#c <- b+ scale_x_continuous(breaks=seq(0, 80000, 10000)) + theme(legend.position = "none")

#png(filename = "fig/rarec.png",
#        width =10, height = 5, units = "in", res= 200)
#c
#dev.off()

set.seed(1234)
PSrare <- rarefy_even_depth(pPS, sample.size=5000, replace=F)

PSrare
PSeim <- subset_taxa(PSrare, family%in%"Eimeriidae")
PScri <- subset_taxa(PSrare, family%in%"Cryptosporidiidae")
#PSsar <- subset_taxa(PS, family%in%"Sarcocystidae")
# for nematoda
PSoxy <- subset_taxa(PSrare, family%in%"Oxyuridae")
PStri <- subset_taxa(PSrare, family%in%"Trichuridae")
PShex <- subset_taxa(PSrare, family%in%"Heteroxynematidae")
#PSstr <- subset_taxa(PS, family%in%"Strongyloididae")
PSasc <- subset_taxa(PSrare, family%in%"Ascaridiidae")
PSspi <- subset_taxa(PSrare, family%in%"Spirocercidae")
PShek <- subset_taxa(PSrare, family%in%"Heterakidae")
# for platyhelminth
PShym <- subset_taxa(PSrare, family%in%"Hymenolepididae")

sample_data(PSrare)$eimASV <- sample_sums(PSeim)
sample_data(PSrare)$criASV <- sample_sums(PScri)
sample_data(PSrare)$oxyASV <- sample_sums(PSoxy)
sample_data(PSrare)$triASV <- sample_sums(PStri)
sample_data(PSrare)$hexASV <- sample_sums(PShex)
sample_data(PSrare)$ascASV <- sample_sums(PSasc)
sample_data(PSrare)$hekASV <- sample_sums(PShek)
sample_data(PSrare)$hymASV <- sample_sums(PShym)
sample_data(PSrare)$spiASV <- sample_sums(PSspi)

PSrare <- subset_taxa(PSrare, !family %in% "Eimeriidae")
PSrare <- subset_taxa(PSrare, !family %in% "Cryptosporidiidae")
PSrare <- subset_taxa(PSrare, !family %in% "Oxyuridae")
PSrare <- subset_taxa(PSrare, !family %in% "Trichuridae")
PSrare <- subset_taxa(PSrare, !family %in% "Heteroxynematidae")
PSrare <- subset_taxa(PSrare, !family %in% "Ascaridiidae")
PSrare <- subset_taxa(PSrare, !family %in% "Heterakidae") 
PSrare <- subset_taxa(PSrare, !family %in% "Hymenolepididae")
PSrare <- subset_taxa(PSrare, !family %in% "Spirocercidae")


Bra_MDSRA <- ordinate(PSrare, "NMDS", "bray")

BC_ieim <-plot_ordination(PSrare, Bra_MDSRA)

BC_ieim <- ggplot(BC_ieim$data, BC_ieim$mapping)+ geom_point(size=1, alpha=0.5)+
    labs(x="NMDS1", y="NMDS2")+
    theme_minimal()

BC_ieim

PSraredis=phyloseq::distance((PSrare), method="bray", type="samples")

sdata=sample_data(PSrare)

permaPSrare=adonis2(PSraredis~
            sdata$Year+
            sdata$Transect+
            sdata$Sex+
            sdata$BMI+
            sdata$HI+
            sdata$eimASV+
            sdata$oxyASV +
            sdata$triASV+
            sdata$hexASV +
            sdata$ascASV +
            sdata$criASV +
            sdata$hekASV +
            sdata$spiASV +
            sdata$hymASV,
        permutations = 10000, method = "bray")
permaPSrare

write.table(permaPSrare, "tmp/permaPSrare.txt")

PSraredisJAC=phyloseq::distance((PSrare), method="jaccard", type="samples")

JpermaPSrare=adonis2(PSraredisJAC~
            sdata$Year+
            sdata$Transect+
            sdata$Sex+
            sdata$BMI+
            sdata$HI+
            sdata$eimASV+
            sdata$oxyASV +
            sdata$triASV+
            sdata$hexASV +
            sdata$ascASV +
            sdata$criASV +
            sdata$hekASV +
            sdata$spiASV +
            sdata$hymASV,
        permutations = 10000, method = "jaccard")
JpermaPSrare
write.table(JpermaPSrare, "tmp/JpermaPSrare.txt")

#for separate domains
erPS <- subset_taxa(PSrare, superkingdom%in%"Eukaryota")
brPS<- subset_taxa(PSrare, superkingdom%in%"Bacteria")

ePSraredisJAC=phyloseq::distance((erPS), method="jaccard", type="samples") 
bPSraredisJAC=phyloseq::distance((brPS), method="jaccard", type="samples") 
sdata=sample_data(PSrare)

eJpermaPSrare=adonis2(ePSraredisJAC~
            sdata$Year+
            sdata$Transect+
            sdata$Sex+
            sdata$BMI+
            sdata$HI+
            sdata$eimASV+
            sdata$oxyASV +
            sdata$triASV+
            sdata$hexASV +
            sdata$ascASV +
            sdata$criASV +
            sdata$hekASV +
            sdata$spiASV +
            sdata$hymASV,
        permutations = 10000, method = "jaccard")
eJpermaPSrare
write.table(eJpermaPSrare, "tmp/EUK_Jac_permaPSrare.txt")
bJpermaPSrare=adonis2(bPSraredisJAC~
            sdata$Year+
            sdata$Transect+
            sdata$Sex+
            sdata$BMI+
            sdata$HI+
            sdata$eimASV+
            sdata$oxyASV +
            sdata$triASV+
            sdata$hexASV +
            sdata$ascASV +
            sdata$criASV +
            sdata$hekASV +
            sdata$spiASV +
            sdata$hymASV,
        permutations = 10000, method = "jaccard")
bJpermaPSrare
write.table(bJpermaPSrare, "tmp/BAC_Jac_permaPSrare.txt")



library(breakaway)
library(DivNet)

#### richness test
library(vegan)
PS <- subset_taxa(PS, !family %in% "Eimeriidae")
PS <- subset_taxa(PS, !family %in% "Cryptosporidiidae")
PS <- subset_taxa(PS, !family %in% "Oxyuridae")
PS <- subset_taxa(PS, !family %in% "Trichuridae")
PS <- subset_taxa(PS, !family %in% "Heteroxynematidae")
PS <- subset_taxa(PS, !family %in% "Ascaridiidae")
PS <- subset_taxa(PS, !family %in% "Heterakidae") 
PS <- subset_taxa(PS, !family %in% "Hymenolepididae")
PS <- subset_taxa(PS, !family %in% "Spirocercidae")

ePS <- subset_taxa(PS, superkingdom%in%"Eukaryota")
bPS<- subset_taxa(PS, superkingdom%in%"Bacteria")

rich = estimate_richness(PS)
rownames(rich) == rownames(sample_data(PS))
sample_data(PS)$richObs <- rich$Observed
sample_data(PS)$richChao <- rich$Chao1
sample_data(PS)$richShan <- rich$Shannon


richdf <- data.frame(
            richChao=sample_data(PS)$richObs,
            richShan=sample_data(PS)$richShan,
            richObs=sample_data(PS)$richObs,
            Year=sample_data(PS)$Year,
            Transect=sample_data(PS)$Transect,
            Sex=sample_data(PS)$Sex,
            BMI=sample_data(PS)$BMI,
            eimASV=sample_data(PS)$eimASV,
            oxyASV=sample_data(PS)$oxyASV,
            triASV=sample_data(PS)$triASV,
            hexASV=sample_data(PS)$hexASV,
            ascASV=sample_data(PS)$ascASV,
            criASV=sample_data(PS)$criASV,
            hekASV=sample_data(PS)$hekASV,
            spiASV=sample_data(PS)$spiASV,
            hymASV=sample_data(PS)$hymASV,
            LibSize=sample_data(PS)$LibSize)

library(MASS)
# Chao1
#find optimal lambda for Box-Cox transformation
bc1 <- boxcox(sample_data(PS)$richChao~
            sample_data(PS)$Year+
            sample_data(PS)$Transect+
            sample_data(PS)$Sex+
            sample_data(PS)$BMI+
            sample_data(PS)$eimASV +
            sample_data(PS)$oxyASV +
            sample_data(PS)$triASV +
            sample_data(PS)$hexASV +
            sample_data(PS)$ascASV +
            sample_data(PS)$criASV +
            sample_data(PS)$hekASV +
            sample_data(PS)$spiASV +
            sample_data(PS)$hymASV +
            sample_data(PS)$LibSize)

lambda1 <- bc1$x[which.max(bc1$y)]
lambda1

bcchaoLM <- lm(((sample_data(PS)$richChao^lambda-1)/lambda)~
            sample_data(PS)$Year+
            sample_data(PS)$Transect+
            sample_data(PS)$Sex+
            sample_data(PS)$BMI+
            sample_data(PS)$eimASV +
            sample_data(PS)$oxyASV +
            sample_data(PS)$triASV +
            sample_data(PS)$hexASV +
            sample_data(PS)$ascASV +
            sample_data(PS)$criASV +
            sample_data(PS)$hekASV +
            sample_data(PS)$spiASV +
            sample_data(PS)$hymASV +
            sample_data(PS)$LibSize)

chaoLM <- lm(richChao~
            Year+
            Transect+
            Sex+
            BMI+
            eimASV +
            oxyASV +
            triASV +
            hexASV +
            ascASV +
            criASV +
            hekASV +
            spiASV +
            hymASV +
            LibSize, data=richdf)


summary(bcchaoLM)
summary(chaoLM)
#plot(chaoLM)

chaoLMsumm <- summary(chaoLM)

sink("tmp/chaoLM.txt")
chaoLMsumm
sink()

chaoLRT <- drop1(chaoLM, test = "Chisq")

sink("tmp/chaoLM_LRT.txt")
chaoLRT
sink()


library(relaimpo)
calc.relimp(chaoLM)
calc.relimp(chaoLM, rela=TRUE)

# Shannon's diversity
shanLM <- lm(richShan~
            Year+
            Transect+
            Sex+
            BMI+
            eimASV +
            oxyASV +
            triASV +
            hexASV +
            ascASV +
            criASV +
            hekASV +
            spiASV +
            hymASV +
            LibSize, data=richdf)
#plot(shanLM)

sink("tmp/shanLM.txt")
summary(shanLM)
sink()

sink("tmp/shanLM_LRT.txt")
drop1(shanLM, test="Chisq")
sink()

#observed richness
obsLM <- lm(richObs~
            Year+
            Transect+
            Sex+
            BMI+
            eimASV +
            oxyASV +
            triASV +
            hexASV +
            ascASV +
            criASV +
            hekASV +
            spiASV +
            hymASV +
            LibSize, data=richdf)

summary(obsLM)

sink("tmp/obsLM.txt")
summary(obsLM)
sink()

sink("tmp/obsLM_LRT.txt")
drop1(obsLM, test="Chisq")
sink()

library(ggeffects)

png(filename = "fig/RichOBS_eim.png",
        width =4, height = 3, units = "in", res= 300)
plot(ggpredict(obsLM, terms = "eimASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Eimeriidae reads", y="ASV richness")
dev.off()

png(filename = "fig/RichOBS_tri.png",
        width =4, height = 3, units = "in", res= 300)
plot(ggpredict(obsLM, terms = "triASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Trichuriidae reads", y="ASV richness")
dev.off()

png(filename = "fig/RichOBS_hek.png",
        width =4, height = 3, units = "in", res= 300)
plot(ggpredict(obsLM, terms = "hekASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Heterakidae reads", y="ASV richness")
dev.off()


## for bacteria and eukarya
erich = estimate_richness(ePS)
brich = estimate_richness(bPS)

Richdf <- data.frame(
            erichChao=erich$Observed,
            erichShan=erich$Shannon,
            erichObs=erich$Chao1,
            brichChao=brich$Observed,
            brichShan=brich$Shannon,
            brichObs=brich$Chao1,
            Year=sample_data(PS)$Year,
            HI=sample_data(PS)$HI,
            Transect=sample_data(PS)$Transect,
            Sex=sample_data(PS)$Sex,
            BMI=sample_data(PS)$BMI,
            eimASV=sample_data(PS)$eimASV,
            oxyASV=sample_data(PS)$oxyASV,
            triASV=sample_data(PS)$triASV,
            hexASV=sample_data(PS)$hexASV,
            ascASV=sample_data(PS)$ascASV,
            criASV=sample_data(PS)$criASV,
            hekASV=sample_data(PS)$hekASV,
            spiASV=sample_data(PS)$spiASV,
            hymASV=sample_data(PS)$hymASV,
            LibSize=sample_data(PS)$LibSize)

eObsLM <- lm(erichObs~
            Year+
            Transect+
            Sex+
            BMI+
            eimASV +
            oxyASV +
            triASV +
            hexASV +
            ascASV +
            criASV +
            hekASV +
            spiASV +
            hymASV +
            LibSize, data=Richdf)

bObsLM <- lm(brichObs~
            Year+
            Transect+
            Sex+
            BMI+
            eimASV +
            oxyASV +
            triASV +
            hexASV +
            ascASV +
            criASV +
            hekASV +
            spiASV +
            hymASV +
            LibSize, data=Richdf)

eShanLM <- lm(erichShan~
            Year+
            Transect+
            Sex+
            BMI+
            eimASV +
            oxyASV +
            triASV +
            hexASV +
            ascASV +
            criASV +
            hekASV +
            spiASV +
            hymASV +
            LibSize, data=Richdf)

bShanLM <- lm(brichShan~
            Year+
            Transect+
            Sex+
            BMI+
            eimASV +
            oxyASV +
            triASV +
            hexASV +
            ascASV +
            criASV +
            hekASV +
            spiASV +
            hymASV +
            LibSize, data=Richdf)

eChaoLM <- lm(erichChao~
            Year+
            Transect+
            Sex+
            BMI+
            eimASV +
            oxyASV +
            triASV +
            hexASV +
            ascASV +
            criASV +
            hekASV +
            spiASV +
            hymASV +
            LibSize, data=Richdf)

bChaoLM <- lm(brichChao~
            Year+
            Transect+
            Sex+
            BMI+
            eimASV +
            oxyASV +
            triASV +
            hexASV +
            ascASV +
            criASV +
            hekASV +
            spiASV +
            hymASV +
            LibSize, data=Richdf)

sink("tmp/bObsLM_LRT.txt")
drop1(bObsLM, test="Chisq")
sink()

sink("tmp/bObsLM.txt")
summary(bObsLM)
sink()

sink("tmp/eObsLM_LRT.txt")
drop1(eObsLM, test="Chisq")
sink()

sink("tmp/eObsLM.txt")
summary(eObsLM)
sink()

sink("tmp/eChaoLM.txt")
summary(eChaoLM)
sink()

sink("tmp/eChaoLM_LRT.txt")
drop1(eChaoLM, test="Chisq")
sink()

sink("tmp/bChaoLM.txt")
summary(bChaoLM)
sink()

sink("tmp/bChaoLM_LRT.txt")
drop1(bChaoLM, test="Chisq")
sink()

sink("tmp/eShanLM.txt")
summary(eShanLM)
sink()

sink("tmp/eShanLM_LRT.txt")
drop1(eShanLM, test="Chisq")
sink()

sink("tmp/bShanLM.txt")
summary(bShanLM)
sink()

sink("tmp/bShanLM_LRT.txt")
drop1(bShanLM, test="Chisq")
sink()


calc.relimp(eChaoLM)

calc.relimp(bChaoLM)

summary(bChaoLM)


##Obs
pdf("fig/figures4man/EukRichOBS_year.pdf",
        width =3, height = 5)
plot(ggpredict(eObsLM, terms = "Year"), add.data = TRUE, show.title=FALSE)+
labs(x="Year of sampling", y="ASV richness")
dev.off()

pdf("fig/figures4man/BacRichOBS_transect.pdf",
        width =3, height = 5)
plot(ggpredict(bObsLM, terms = "Transect"), add.data = TRUE, show.title=FALSE)+
labs(x="Transect", y="ASV richness")
dev.off()

pdf("fig/figures4man/EukRichOBS_eim.pdf",
        width =3, height = 5)
plot(ggpredict(eObsLM, terms = "eimASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Eimeriidae reads", y="ASV richness")
dev.off()

pdf("fig/figures4man/EukRichOBS_oxy.pdf",
        width =3, height = 5)
plot(ggpredict(eObsLM, terms = "oxyASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Oxyuridae reads", y="ASV richness")
dev.off()

pdf("fig/figures4man/EukRichOBS_tri.pdf",
        width =3, height = 5)
plot(ggpredict(eObsLM, terms = "triASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Trichuridae reads", y="ASV richness")
dev.off()

pdf("fig/figures4man/EukRichOBS_hex.pdf",
        width =3, height = 5)
plot(ggpredict(eObsLM, terms = "hexASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Heteroxynematidae reads", y="ASV richness")
dev.off()

pdf("fig/figures4man/BacRichOBS_eim.pdf",
        width =3, height = 5)
plot(ggpredict(bObsLM, terms = "eimASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Eimeriidae reads", y="ASV richness")
dev.off()

pdf("fig/figures4man/BacRichOBS_tri.pdf",
        width =3, height = 5)
plot(ggpredict(bObsLM, terms = "triASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Trichuriidae reads", y="ASV richness")
dev.off()

## Shan
png(filename = "fig/BacRichShan_eim.png",
        width =3, height = 5, units = "in", res= 300)
plot(ggpredict(bShanLM, terms = "eimASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Eimeriidae reads", y="ASV richness")
dev.off()

png(filename = "fig/BacRichShan_tri.png",
        width =3, height = 5, units = "in", res= 300)
plot(ggpredict(bShanLM, terms = "triASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Trichuridae reads", y="ASV richness")
dev.off()

png(filename = "fig/EukRichShan_hex.png",
        width =3, height = 5, units = "in", res= 300)
plot(ggpredict(eShanLM, terms = "hexASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Heteroxynematidae reads", y="ASV richness")
dev.off()

png(filename = "fig/EukRichShan_tri.png",
        width =3, height = 5, units = "in", res= 300)
plot(ggpredict(eShanLM, terms = "triASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Trichuridae reads", y="ASV richness")
dev.off()

                                        #Chao
png(filename = "fig/BacRichChao_eim.png",
        width =3, height = 3, units = "in", res= 300)
plot(ggpredict(bChaoLM, terms = "eimASV"), add.data = TRUE, show.title=FALSE)+
    labs(x="Eimeridae reads", y="ASV richness")+ ggtitle("Bacterial community")
dev.off()

png(filename = "fig/BacRichChao_tri.png",
        width =3, height = 3, units = "in", res= 300)
plot(ggpredict(bChaoLM, terms = "triASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Trichuridae reads", y="ASV richness")
dev.off()

png(filename = "fig/EukRichChao_eim.png",
        width =3, height = 3, units = "in", res= 300)
plot(ggpredict(eChaoLM, terms = "eimASV"), add.data = TRUE, show.title=FALSE)+
    labs(x="Eimeriidae reads", y="ASV richness")+ ggtitle("Eukaryotic community")
dev.off()

png(filename = "fig/EukRichChao_tri.png",
        width =3, height = 5, units = "in", res= 300)
plot(ggpredict(eChaoLM, terms = "triASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Trichuridae reads", y="ASV richness")
dev.off()

png(filename = "fig/EukRichChao_oxy.png",
        width =3, height = 5, units = "in", res= 300)
plot(ggpredict(eChaoLM, terms = "oxyASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Oxyuridae reads", y="ASV richness")
dev.off()

png(filename = "fig/EukRichChao_hex.png",
        width =3, height = 5, units = "in", res= 300)
plot(ggpredict(eChaoLM, terms = "hexASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Heteroxynematidae reads", y="ASV richness")
dev.off()

png(filename = "fig/EukRichChao_hym.png",
        width =3, height = 5, units = "in", res= 300)
plot(ggpredict(eChaoLM, terms = "hymASV"), add.data = TRUE, show.title=FALSE)+
labs(x="Hymenolepididae reads", y="ASV richness")
dev.off()


## plotting DA
# Differencial abundance analysis
#remotes::install <- github("FrederickHuangLin/ANCOMBC")

################ abundance analysis


### Ancom is very cool for dichotomous variables
library(ANCOMBC)
library(DT)

library(RColorBrewer)
library(ggpubr)

# Plot differencial analysis with ANCOMBC
plotDA= function(PS, name) {
out= ancombc(phyloseq=PS, formula=name, lib_cut=1)
res=out$res
res=as.data.frame(res)
names(res)=c("log.changes", "se", "W", "p.val", "p-adj", "diff.ab")
phyla = data.frame(phylum=PS@tax_table[,c(2, 5)])
res=merge(phyla, res, by=0)
#str(res)
alpha.adj=0.05/nrow(res)
critic.val=qnorm(1-alpha.adj/2)
res$ci.lo.adj=res$log.changes - critic.val*res$se
res$ci.up.adj=res$log.changes + critic.val*res$se
resS=res[res$diff.ab==TRUE,]
nb.col <- nrow(resS)
mycolor <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.col)
oxyDA=ggbarplot(resS, x="Row.names", y="log.changes",
          fill= "phylum.genus",
          palette="mycolor",
          sor.by.groups=FALSE,
          sort.val="asc",
          xlab="ASV",
          legend.title="genus",
          rotate=TRUE
          ) +
    rremove("y.text")+
    geom_errorbar(aes(ymin=ci.lo.adj, ymax=ci.up.adj), width=0.2)+
    theme(legend.title=element_text(size=8),
          legend.text=element_text(size=6))+
    guides(fill = guide_legend(ncol=4, bycol=TRUE))
oxyDA
}

eimDA=plotDA(PS, "iEimASV")
oxyDA=plotDA(PS, "ioxyASV")
criDA=plotDA(PS, "icriASV")
triDA=plotDA(PS, "itriASV")
hexDA=plotDA(PS, "ihexASV")
ascDA=plotDA(PS, "iascASV")
hekDA=plotDA(PS, "ihekASV")
hymDA=plotDA(PS, "ihymASV")

png(filename = "fig/DAeim.png",
        width =4, height = 6, units = "in", res= 400)
eimDA
dev.off()

png(filename = "fig/DAoxy.png",
        width =4, height = 6, units = "in", res= 400)
oxyDA
dev.off()

png(filename = "fig/DAtri.png",
        width =4, height = 6, units = "in", res= 400)
triDA
dev.off()

png(filename = "fig/DAhex.png",
        width =4, height = 6, units = "in", res= 400)
hexDA
dev.off()

png(filename = "fig/DAasc.png",
        width =4, height = 6, units = "in", res= 400)
ascDA
dev.off()

png(filename = "fig/DAhek.png",
        width =4, height = 2, units = "in", res= 400)
hekDA
dev.off()

png(filename = "fig/DAhym.png",
        width =4, height = 3, units = "in", res= 400)
hymDA
dev.off()

png(filename = "fig/DAnem.png",
        width =8, height = 4, units = "in", res= 400)
ggarrange(oxyDA, triDA, hexDA, ascDA, ncol=4,
          labels=c("Oxyuriidae", "Trichuridae", "Heteroxynematidae", "Ascaridiidae"))
dev.off()

png(filename = "fig/DAplat.png",
        width =4, height = 4, units = "in", res= 400)
hymDA
dev.off()

#### Calculate normalized stochasticity ration(NST)
library("NST")
tNST(otups, sample_data(PS)$EimASV)

