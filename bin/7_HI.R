library(ggplot2)
#library(reshape)
library(phyloseq, lib="/usr/local/lib/R/site-library/")
#library(data.table)
#library(microbiome)
#library(plyr)
#library(dplyr)
#library(gridExtra)
#library(grid)
#library(vegan)
#library(tidyr)
#library(tidyverse)
#library(ggpubr)
#library(ggeffects, lib="/usr/local/lib/R/site-library/")
library(reshape)
library(data.table)
library(parallel)
library(microbiome)
library("pheatmap")
library(dplyr)
library(gridExtra)
library(grid)
library("fitdistrplus")
library("optimx")
#library(FSA, lib="/usr/local/lib/R/site-library/")

#devtools::install_github("alicebalard/parasiteLoad")

library(parasiteLoad)

### Functions
source("bin/TestDistributions.R")
source("bin/bananaplotNoCI.R")

PS=readRDS(file="tmp/PSpre.R")

pPS=readRDS(file="tmp/pPSpre.R")

PSTSS = transform_sample_counts(pPS, function(x) x / sum(x) )

#sample_data(pPS)$HI

set.seed(1234)
PSrare <- rarefy_even_depth(PS, sample.size=5000, replace=F)

erPS <- subset_taxa(PSrare, superkingdom%in%"Eukaryota")
brPS<- subset_taxa(PSrare, superkingdom%in%"Bacteria")

ePS <- subset_taxa(PS, superkingdom%in%"Eukaryota")
bPS<- subset_taxa(PS, superkingdom%in%"Bacteria")

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

sample_data(PSrare)$all <- "one"

source("bin/Functions.R")

summarize_taxa(PSeim, "genus", "all")
summarize_taxa(PScri, "genus", "all")
summarize_taxa(PSoxy, "genus", "all")
summarize_taxa(PStri, "genus", "all")
summarize_taxa(PShex, "genus", "all")
summarize_taxa(PShek, "genus", "all")
summarize_taxa(PSasc, "genus", "all")
summarize_taxa(PShym, "genus", "all")
summarize_taxa(PSspi, "genus", "all")

sample_data(PSrare)$eimASV <- sample_sums(PSeim)
sample_data(PSrare)$criASV <- sample_sums(PScri)
sample_data(PSrare)$oxyASV <- sample_sums(PSoxy)
sample_data(PSrare)$triASV <- sample_sums(PStri)
sample_data(PSrare)$hexASV <- sample_sums(PShex)
sample_data(PSrare)$ascASV <- sample_sums(PSasc)
sample_data(PSrare)$hekASV <- sample_sums(PShek)
sample_data(PSrare)$hymASV <- sample_sums(PShym)
sample_data(PSrare)$spiASV <- sample_sums(PSspi)

PSeimg <- subset_taxa(pPS, genus%in%"Eimeria")
sample_data(pPS)$eimASVg <- sample_sums(PSeimg)

#library(microbiome)
#PSTSS = transform(pPS, "compositional")

#PSeim <- subset_taxa(PSTSS, family%in%"Eimeriidae")
#sample_data(PSTSS)$eimASV <- sample_sums(PSeim)
#e2PS <- subset_taxa(PSTSS, superkingdom%in%"Eukaryota")
#b2PS<- subset_taxa(PSTSS, superkingdom%in%"Bacteria")

erich = estimate_richness(ePS)
brich = estimate_richness(bPS)

errich = estimate_richness(erPS)
brrich = estimate_richness(brPS)

Richdf <- data.frame(
            erichChao=round(erich$Observed*100),
            erichShan=round(erich$Shannon*100),
            erichObs=round(erich$Chao1*100),
            brichChao=round(brich$Observed*100),
            brichShan=round(brich$Shannon*100),
            brichObs=round(brich$Chao1*100),
            Year=sample_data(PS)$Year,
            Transect=sample_data(PS)$Transect,
            Sex=sample_data(PS)$Sex,
            HI=sample_data(PS)$HI,
            BMI=sample_data(PS)$BMI,
    eimASV=sample_data(PS)$eimASV,
    eimOPG=sample_data(PS)$OPG_Eimeria,
            logeimASV=log10(sample_data(PS)$eimASV+1),
            oxyASV=sample_data(PS)$oxyASV,
            triASV=sample_data(PS)$triASV,
            hexASV=sample_data(PS)$hexASV,
            ascASV=sample_data(PS)$ascASV,
            criASV=sample_data(PS)$criASV,
            hekASV=sample_data(PS)$hekASV,
            spiASV=sample_data(PS)$spiASV,
            hymASV=sample_data(PS)$hymASV,
            LibSize=sample_data(PS)$LibSize)


sample_data(PS)$OPG_Eimeria

pPSdf <- data.frame(
            Year=sample_data(pPS)$Year,
            Transect=sample_data(pPS)$Transect,
            Sex=sample_data(pPS)$Sex,
            HI=sample_data(pPS)$HI,
            BMI=sample_data(pPS)$BMI,
    eimASV=sample_data(pPS)$eimASV,
    eimASVg=sample_data(pPS)$eimASVg,
            oxyASV=sample_data(pPS)$oxyASV,
            triASV=sample_data(pPS)$triASV,
            hexASV=sample_data(pPS)$hexASV,
            ascASV=sample_data(pPS)$ascASV,
            criASV=sample_data(pPS)$criASV,
            hekASV=sample_data(pPS)$hekASV,
            spiASV=sample_data(pPS)$spiASV,
            hymASV=sample_data(pPS)$hymASV,
    LibSize=sample_data(pPS)$LibSize)


relAbdf <- data.frame(
            Year=sample_data(PSTSS)$Year,
            Transect=sample_data(PSTSS)$Transect,
            Sex=sample_data(PSTSS)$Sex,
            HI=sample_data(PSTSS)$HI,
    BMI=sample_data(PSTSS)$BMI,
    logeimASV=round(log10(sample_data(PSTSS)$eimASV+1)*10000),
    logeimgASV=round(log10(sample_data(PSTSS)$eimASVg+1)*10000),
    logoxyASV=round(log10(sample_data(PSTSS)$oxyASV+1)*10000),
    logtriASV=round(log10(sample_data(PSTSS)$triASV+1)*10000),
    loghexASV=round(log10(sample_data(PSTSS)$hexASV+1)*10000),
    logascASV=round(log10(sample_data(PSTSS)$ascASV+1)*10000),
    logcriASV=round(log10(sample_data(PSTSS)$criASV+1)*10000),
    loghekASV=round(log10(sample_data(PSTSS)$hekASV+1)*10000),
    logspiASV=round(log10(sample_data(PSTSS)$spiASV+1)*10000),
    loghymASV=round(log10(sample_data(PSTSS)$hymASV+1)*10000),
        eimASV=round((sample_data(PSTSS)$eimASV)*10000),
    eimgASV=round((sample_data(PSTSS)$eimASVg)*10000),
    oxyASV=round((sample_data(PSTSS)$oxyASV)*10000),
    triASV=round((sample_data(PSTSS)$triASV)*10000),
    hexASV=round((sample_data(PSTSS)$hexASV)*10000),
    ascASV=round((sample_data(PSTSS)$ascASV)*10000),
    criASV=round((sample_data(PSTSS)$criASV)*10000),
    hekASV=round((sample_data(PSTSS)$hekASV)*10000),
    spiASV=round((sample_data(PSTSS)$spiASV)*10000),
    hymASV=round((sample_data(PSTSS)$hymASV)*10000),
    LibSize=sample_data(PSTSS)$LibSize)

Richrdf <- data.frame(
            erichChao=round(errich$Observed*100),
            erichShan=round(errich$Shannon*100),
            erichObs=round(errich$Chao1*100),
            brichChao=round(brrich$Observed*100),
            brichShan=round(brrich$Shannon*100),
            brichObs=round(brrich$Chao1*100),
            Year=sample_data(PSrare)$Year,
            Transect=sample_data(PSrare)$Transect,
            Sex=sample_data(PSrare)$Sex,
            HI=sample_data(PSrare)$HI,
            BMI=sample_data(PSrare)$BMI,
            eimASV=sample_data(PSrare)$eimASV,
            logeimASV=log10(sample_data(PSrare)$eimASV+1),
            oxyASV=sample_data(PSrare)$oxyASV,
            triASV=sample_data(PSrare)$triASV,
            hexASV=sample_data(PSrare)$hexASV,
            ascASV=sample_data(PSrare)$ascASV,
            criASV=sample_data(PSrare)$criASV,
            hekASV=sample_data(PSrare)$hekASV,
            spiASV=sample_data(PSrare)$spiASV,
            hymASV=sample_data(PSrare)$hymASV,
            LibSize=sample_data(PSrare)$LibSize)


hist(log10(Richrdf$eimASV+1))

hist(relAbdf$eimASV)

## remove zeros
#eim <- paradf[paradf$eimASV>0,]
eim <- Richdf[Richdf$eimASV>0,]
oxy <- Richdf[Richdf$oxyASV>0,]
tri <- Richdf[Richdf$triASV>0,]
hex <- Richdf[Richdf$hexASV>0,]
asc <- Richdf[Richdf$ascASV>0,]
cri <- Richdf[Richdf$criASV>0,]
hek <- Richdf[Richdf$hekASV>0,]
hym <- Richdf[Richdf$hymASV>0,]
spi <- Richdf[Richdf$spiASV>0,]
eimOo <-  na.omit(Richdf[Richdf$eimOPG>0,])

eimOo$eimOPG <- round(eimOo$eimOPG)

eimOo$eimOPG

eimr <- Richrdf[Richrdf$eimASV>0,]
oxyr <- Richrdf[Richrdf$oxyASV>0,]
trir <- Richrdf[Richrdf$triASV>0,]
hexr <- Richrdf[Richrdf$hexASV>0,]
ascr <- Richrdf[Richrdf$ascASV>0,]
crir <- Richrdf[Richrdf$criASV>0,]
hekr <- Richrdf[Richrdf$hekASV>0,]
hymr <- Richrdf[Richrdf$hymASV>0,]
spir <- Richrdf[Richrdf$spiASV>0,]

eimtss <- relAbdf[relAbdf$logeimASV>0,]
oxytss <- relAbdf[relAbdf$logoxyASV>0,]
tritss <- relAbdf[relAbdf$logtriASV>0,]
hextss <- relAbdf[relAbdf$loghexASV>0,]
asctss <- relAbdf[relAbdf$logascASV>0,]
critss <- relAbdf[relAbdf$logcriASV>0,]
hektss <- relAbdf[relAbdf$loghekASV>0,]
hymtss <- relAbdf[relAbdf$loghymASV>0,]
spitss <- relAbdf[relAbdf$logspiASV>0,]

eimp <- pPSdf[pPSdf$eimASV>0,]
oxyp <- pPSdf[pPSdf$oxyASV>0,]
trip <- pPSdf[pPSdf$triASV>0,]
hexp <- pPSdf[pPSdf$hexASV>0,]
ascp <- pPSdf[pPSdf$ascASV>0,]
crip <- pPSdf[pPSdf$criASV>0,]
hekp <- pPSdf[pPSdf$hekASV>0,]
hymp <- pPSdf[pPSdf$hymASV>0,]
spip <- pPSdf[pPSdf$spiASV>0,]



findGoodDist(eimtss$logeimASV,distribs = c("normal", "negative binomial", "poisson"),
                          distribs2 = c("norm", "nbinom", "pois"))

head(eimtss)

# for eimeria (test)
fitp <- parasiteLoad::analyse(data = eimp, response = "eimASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

fittss <- parasiteLoad::analyse(data = eimtss, response = "eimASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

fit <- parasiteLoad::analyse(data = eim, response = "eimASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

fitopg <- parasiteLoad::analyse(data = eimOo, response = "eimOPG", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

fitr <- parasiteLoad::analyse(data = eimr, response = "eimASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

fitOPG <- parasiteLoad::analyse(data = PS, response = "OPG_Eimeira", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

names(sample_data(PS))

# for oxy (test)
ofitp <- parasiteLoad::analyse(data = oxyp, response = "oxyASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

ofittss <- parasiteLoad::analyse(data = oxytss, response = "oxyASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

ofit <- parasiteLoad::analyse(data = oxy, response = "oxyASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

ofitr <- parasiteLoad::analyse(data = oxyr, response = "oxyASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")


# for tri (test)
tfitp <- parasiteLoad::analyse(data = trip, response = "triASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

tfittss <- parasiteLoad::analyse(data = tritss, response = "triASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

tfit <- parasiteLoad::analyse(data = tri, response = "triASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

tfitr <- parasiteLoad::analyse(data = trir, response = "triASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

#####
eimtss.plot <- parasiteLoad::bananaPlot(mod = fittss$H2,
                                       data = eimtss,
                                       response = "eimASV",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = T,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
pdf("fig/figures4man/Eimtss_banana.pdf",
        width =5, height = 5)
eimtss.plot
dev.off()

fit <- parasiteLoad::analyse(data = oxy, response = "oxyASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")
oxy.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = oxy,
                                       response = "oxyASV",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))

pdf("fig/figures4man/Oxy_banana.pdf",
        width =5, height = 5)
oxy.plot
dev.off()

fit <- parasiteLoad::analyse(data = tri, response = "triASV", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")
tri.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = tri,
                                       response = "triASV",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))

pdf("fig/figures4man/Tri_banana.pdf",
        width =5, height = 5)
tri.plot
dev.off()

## Observed richness
# for eukaryotes
fit <- parasiteLoad::analyse(data = Richdf, response = "erichObs", model = "weibull", group = "Sex", hybridIndex = "HI", myparamBounds="default")

fit <- parasiteLoad::analyse(data = Richrdf, response = "erichObs", model = "weibull", group = "Sex", hybridIndex = "HI", myparamBounds="default")


erichObs.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = Richdf,
                                       response = "erichObs",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))

erichObs.plot


pdf("fig/figures4man/EukRichOBS_banana.pdf",
        width =5, height = 5)
erichObs.plot
dev.off()

# for bacteria
fit <- parasiteLoad::analyse(data = Richdf, response = "brichObs", model = "weibull", group = "Sex", hybridIndex = "HI", myparamBounds="default")

brichObs.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = Richdf,
                                       response = "brichObs",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
brichObs.plot

pdf("fig/figures4man/BacRichOBS_banana.pdf",
        width =5, height = 5)
brichObs.plot
dev.off()

## Shannon richness
# for eukaryotes
fit <- parasiteLoad::analyse(data = Richdf, response = "erichShan", model = "weibull", group = "Sex", hybridIndex = "HI", myparamBounds="default")

fit <- parasiteLoad::analyse(data = Richrdf, response = "erichShan", model = "weibull", group = "Sex", hybridIndex = "HI", myparamBounds="default")


erichShan.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = Richdf,
                                       response = "erichShan",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
erichShan.plot


pdf("fig/figures4man/EukRichShan_banana.pdf",
        width =5, height = 5)
erichShan.plot
dev.off()

# for bacteria
fit <- parasiteLoad::analyse(data = Richdf, response = "brichShan", model = "weibull", group = "Sex", hybridIndex = "HI", myparamBounds="default")

brichShan.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = Richdf,
                                       response = "brichShan",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
brichShan.plot

pdf("fig/figures4man/BacRichShan_banana.pdf",
        width =3, height = 3)
brichShan.plot
dev.off()

## Chao1 richness
# for eukaryotes
fit <- parasiteLoad::analyse(data = Richdf, response = "erichChao", model = "weibull", group = "Sex", hybridIndex = "HI", myparamBounds="default")

erichChao.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = Richdf,
                                       response = "erichChao",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
erichChao.plot


pdf("fig/figures4man/EukRichChao_banana.pdf",
        width =3, height = 3)
erichChao.plot
dev.off()

# for bacteria
fit <- parasiteLoad::analyse(data = Richdf, response = "brichChao", model = "weibull", group = "Sex", hybridIndex = "HI", myparamBounds="default")

brichChao.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = Richdf,
                                       response = "brichChao",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
brichChao.plot

pdf("fig/figures4man/BacRichChao_banana.pdf",
        width =3, height = 3)
brichChao.plot
dev.off()

##### for overall diversity
rich = estimate_richness(PS)

sample_data(PS)$richObs <- rich$Observed
sample_data(PS)$richChao <- rich$Chao1
sample_data(PS)$richShan <- rich$Shannon

rich1df <- data.frame(
            richChao=sample_data(PS)$richObs,
            richShan=sample_data(PS)$richShan,
            richObs=sample_data(PS)$richObs,
            Sex=sample_data(PS)$Sex,
            HI=sample_data(PS)$HI)

fit <- parasiteLoad::analyse(data = rich1df, response = "richChao", model = "weibull", group = "Sex", hybridIndex = "HI", myparamBounds="default")
richChao.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = rich1df,
                                       response = "richChao",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
richChao.plot


fit <- parasiteLoad::analyse(data = rich1df, response = "richShan", model = "weibull", group = "Sex", hybridIndex = "HI", myparamBounds="default")
richShan.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = rich1df,
                                       response = "richShan",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
richShan.plot

#### Taxa from mva
# eimeria:
eim1 <- "GGATGAACGCTAGCTACAGGCTTAACACATGCAAGTCGAGGGGCAGCATGGCCTATCTTTCGGGATGGGCCGATGGCGACCGGCGCACGGGTGAGTAACGCGTATCCAACCTTCCCTTTACTGGGGTCCAGCCCGTCGAAAGGCGGATTAATCCCCCATGTTCTCCGTCCCGGACATCTGTGTCGGAGCAAAGATTTATCGGTAAAGGATGGGGATGCGTCCGATTAGCTTGTTGGCGGGGTAACGGCCCACCAAGGCATCGATCGGTAGGGGTTCTGAGAGGAAGGTCCCCCACACTGGAACTGAGACACGGTCCAGA"

eim2 <- "GTGCATGTGTAAGTATGAACTAATTCAGACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTATCTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATACGTGCAACAAACCCCGACTTCTGGAAGGGATGCATTTATTAGATAAAAGGTCAACACAGGNNNNNNNNNNGCATTTATTAGATAAAAGGTCAACACAGGCTCTGCCTGTTGCTTTGATGATTCATGATAACTCGTCGGATCGCACGGCCTTTGTGCCGGCGACGCATCATTCAAATTTCTGCCCTATCAACTTTCGATGGTAGGATAGTGGCCTACCATGGTGGTGACGGGTGACGGAGAATTAGGGTTCGATT"

HI1 <- "GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCTTGCCTAGATGATTTTAGTGCTTGCACTAAATGAAACTAGATACAAGCGAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCCAAGAGACTGGGATAACACCTGGAAACAGATGCTAATACCGGATAACAACACTAGACGCATGTCTAGAGTTTAAAAGATGGTTCTGCTATCACTCTTGGATGGACCTGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCC"

HI2 <- "TGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATTGTTAGTTGCCAGCATTAAGTTGGGCACTCTAGCAAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACGGTACAACGANNNNNNNNNNCTGCAACTCGCCTACATGAAGTCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGAGAGTTTGTAACACCCAAAGCCGGTGGGGTAACCTTTTGGAGCCAGCCGTCTAAGGTGGGACAGATGATTAGGGTGA"

tax_table(pPS)[eim1, 6] # Prevotella NA

tax_table(pPS)[eim2, 6] # Glycine NA

tax_table(pPS)[HI1, 6] # Lactobacillus NA

tax_table(pPS)[HI2, 6] # Lactobacillus NA

tax_table(pPS)[eim6, 6] # 

### selected taxa from Victor
PSbla <- subset_taxa(pPS, genus%in%"Blautia")
PSbact <- subset_taxa(pPS, genus%in%"Bacteroides")
PSMur <- subset_taxa(pPS, genus%in%"Muribaculum")
PSpar <- subset_taxa(pPS, genus%in%"Parasutterella")
PSpar <- subset_taxa(pPS, genus%in%"Parasutterella")

PSeimG <- subset_taxa(pPS, genus%in%"Eimeria")

PSPrevo <- prune_taxa(eim1, pPS)
PSGlyc <- prune_taxa(eim2, pPS)
PSLac1 <- prune_taxa(HI1, pPS)
PSLac2 <- prune_taxa(HI2, pPS)


Vdf <- data.frame(bla=sample_sums(PSbla),
                  bact=sample_sums(PSbact),
                  eimG=sample_sums(PSeimG),
                  mur=sample_sums(PSMur),
                  par=sample_sums(PSpar),
                  prevo= sample_sums(PSPrevo),
                  gly= sample_sums(PSGlyc),
                  lac1=sample_sums(PSLac1),
                  lac2=sample_sums(PSLac2),
                   HI=sample_data(pPS)$HI,
                   Sex=sample_data(pPS)$Sex)


Vdfbac <- Vdf[Vdf$bact>0,]

Vdfbla <- Vdf[Vdf$bla>0,]

VdfeimG <- Vdf[Vdf$eimG>0,]

Vdfmur <- Vdf[Vdf$mur>0,]
Vdfpar <- Vdf[Vdf$par>0,]
Vdfprevo <- Vdf[Vdf$prevo>0,]
Vdfgly <- Vdf[Vdf$gly>0,]
Vdflac1 <- Vdf[Vdf$lac1>0,]
Vdflac2 <- Vdf[Vdf$lac2>0,]



fit <- parasiteLoad::analyse(data = Vdfbac, response = "bact", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Vbact.plot <- parasiteLoad::bananaPlot(mod = fit$H1,
                                       data = Vdfbac,
                                       response = "bact",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = T,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
Vbact.plot


fit <- parasiteLoad::analyse(data = Vdfmur, response = "mur", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Vmur.plot <- parasiteLoad::bananaPlot(mod = fit$H1,
                                       data = Vdfmur,
                                       response = "mur",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = T,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))

Vmur.plot

pdf("fig/figures4man/Muribaculum_banana.pdf",
        width =5, height = 4)
Vmur.plot
dev.off()

fit <- parasiteLoad::analyse(data = Vdfpar, response = "par", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Vpar.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = Vdfpar,
                                       response = "par",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = T,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
Vpar.plot

pdf("fig/figures4man/Parasutterella_banana.pdf",
        width =5, height = 4)
Vpar.plot
dev.off()


fit <- parasiteLoad::analyse(data = Vdfprevo, response = "prevo", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Vprevo.plot <- parasiteLoad::bananaPlot(mod = fit$H1,
                                       data = Vdfprevo,
                                       response = "prevo",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
Vprevo.plot

fit <- parasiteLoad::analyse(data = Vdfgly, response = "gly", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Vgly.plot <- parasiteLoad::bananaPlot(mod = fit$H1,
                                       data = Vdfgly,
                                       response = "gly",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
Vgly.plot

fit <- parasiteLoad::analyse(data = Vdflac1, response = "lac1", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Vlac1.plot <- parasiteLoad::bananaPlot(mod = fit$H1,
                                       data = Vdflac1,
                                       response = "lac1",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
Vlac1.plot

fit <- parasiteLoad::analyse(data = Vdflac2, response = "lac2", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Vlac2.plot <- parasiteLoad::bananaPlot(mod = fit$H1,
                                       data = Vdflac2,
                                       response = "lac2",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
Vlac2.plot


fit <- parasiteLoad::analyse(data = VdfeimG, response = "eimG", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

VeimG.plot <- parasiteLoad::bananaPlot(mod = fit$H1,
                                       data = VdfeimG,
                                       response = "eimG",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = T,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
VeimG.plot



fit <- parasiteLoad::analyse(data = Vdfbla, response = "bla", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Vdf$bla

Vbla.plot <- parasiteLoad::bananaPlot(mod = fit$H1,
                                       data = Vdfbla,
                                       response = "bla",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = T,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))

Vbla.plot


pdf("fig/figures4man/PrevotellaNA_banana.pdf",
        width =5, height = 4)
Vprevo.plot
dev.off()
pdf("fig/figures4man/GlycineNA_banana.pdf",
        width =5, height = 4)
Vgly.plot
dev.off()
pdf("fig/figures4man/LactobacillusNA1_banana.pdf",
        width =5, height = 4)
Vlac1.plot
dev.off()
pdf("fig/figures4man/LactobacillusNA2_banana.pdf",
        width =5, height = 4)
Vlac2.plot
dev.off()


##### Selected taxa from cluster analysis

cluster4 <- readRDS(file="tmp/cluster4.RDS")

PSc4 <- prune_taxa(cluster4, pPS)
PS4eim <- subset_taxa(PSc4, family%in%"Eimeriidae")
PS4sacc <- subset_taxa(PSc4, family%in%"Saccharomycetaceae")
c4df <- data.frame(eim=sample_sums(PS4eim),
                   sacc=sample_sums(PS4sacc),
                   HI=sample_data(PS4eim)$HI,
                   Sex=sample_data(PS4eim)$Sex)

c4dfeim <- c4df[c4df$eim>0,]

c4dfsacc <- c4df[c4df$sacc>0,]

fit <- parasiteLoad::analyse(data = c4dfeim, response = "eim", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

eim4.plot <- parasiteLoad::bananaPlot(mod = fit$H1,
                                       data = c4dfeim,
                                       response = "eim",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))

eim4.plot

pdf("fig/figures4man/Eimeria4_bananawZ.pdf",
        width =5, height = 4)
eim4.plot
dev.off()

fit <- parasiteLoad::analyse(data = c4dfsacc, response = "sacc", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

sacc4.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = c4dfsacc,
                                       response = "sacc",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))

pdf("fig/figures4man/Saccha4_banana.pdf",
        width =5, height = 4)
sacc4.plot
dev.off()


### cluster 9


cluster9 <- readRDS(file="tmp/cluster9.RDS")
PSc9 <- prune_taxa(cluster9, pPS)
PS9eim <- subset_taxa(PSc9, family%in%"Eimeriidae")
PS9sacc <- subset_taxa(PSc9, family%in%"Saccharomycetaceae")
c9df <- data.frame(eim=sample_sums(PS9eim),
                   sacc=sample_sums(PS9sacc),
                   HI=sample_data(PS9eim)$HI,
                   Sex=sample_data(PS9eim)$Sex)

c9dfeim <- c9df[c9df$eim>0,]

c9dfsacc <- c9df[c9df$sacc>0,]

fit <- parasiteLoad::analyse(data = c9dfeim, response = "eim", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

eim9.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = c9dfeim,
                                       response = "eim",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))
pdf("fig/figures4man/Eimeria9_banana.pdf",
        width =5, height = 4)
eim9.plot
dev.off()

fit <- parasiteLoad::analyse(data = c9dfsacc, response = "sacc", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

sacc9.plot <- parasiteLoad::bananaPlot(mod = fit$H0,
                                       data = c9dfsacc,
                                       response = "sacc",
                                       hybridIndex=seq(0,1,0.005),
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#E69F00"))

pdf("fig/figures4man/Saccha9_banana.pdf",
        width =5, height = 4)

sacc9.plot

dev.off()
