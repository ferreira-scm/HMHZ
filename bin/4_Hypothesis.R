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

PS=readRDS(file="tmp/PSpre.R")
PSgen=tax_glom(PS, taxrank="genus")

PSphy=tax_glom(PS, taxrank="phylum")

## prevalence plotting
prevdf = apply(X = otu_table(PS),
               MARGIN = ifelse(taxa_are_rows(PS), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(PS),
                    tax_table(PS))
plyr::ddply(prevdf, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Subset to the remaining phyla
prevdf1 = subset(prevdf, phylum %in% get_taxa_unique(PS, "phylum"))
pdf("fig/prevPhyla.pdf")
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(PS),color=phylum)) +
# Include a guess for parameter
    geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
    geom_point(size = 2, alpha = 0.7) +
    scale_x_log10() +
    xlab("Total Abundance") +
    ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~phylum) + theme(legend.position="none")
dev.off()

# Abundance plotting at phylum level
## Phylum
source("/SAN/Susanas_den/HMHZ/bin/Functions.R")
relphy=summarize_taxa(PSphy, "phylum", "Mouse_ID")
phylnames=unique(relphy$phylum) #36 phyla
phylnames
relphySorted <- relphy %>%
    mutate(phylum = fct_reorder(phylum, meanRA))
relphySorted %>%
    group_by(phylum)%>%
    mutate(phylum_avg=mean(meanRA))%>%
    mutate(Main_phylum=phylum_avg>=0.01) -> relphySorted
relphySorted%>%
    mutate(phylum=fct_reorder(phylum, meanRA, .fun="mean"))%>%
    summarise(no_rows = length(phylum))%>%
    mutate(Prev_phylum= (no_rows/657)*100)%>%
    mutate(High_prev= Prev_phylum>=30)-> relphyPrev
(relphySorted)
ggplot(relphySorted, aes(x = phylum, y = meanRA, color= phylum)) +
    coord_flip() +
        scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005)) +
    labs(x = NULL, y = "Relative abundance")+
          theme(legend.position = "none")+
      #theme_classic()+
      #geom_hline(aes(yintercept = sample_avg), color = "gray70", size = 0.6) +
      geom_jitter(size = 2, alpha = 0.5, width = 0.2)+
    stat_summary(fun = mean, geom = "point", size = 1, color= "black")+
    theme_minimal()+
    theme(legend.position="non")-> phyAb
png(filename = "fig/Abundance_phyla.png",
        width =7, height = 5, units = "in", res= 300)
phyAb
dev.off()

## Testing Eimeriidae correlation with OPG
cor.test(PSTSS@sam_data$OPG_Eimeria, rowSums(otu_table(subset_taxa(PSTSS, family%in%"Eimeriidae"))), method ="spearman")
png(filename = "fig/Oocysts_Eimeriidae.png",
        width =4, height = 5, units = "in", res= 300)
ggplot(PSTSS@sam_data, aes(x=PSTSS@sam_data$OPG_Eimeria)) +
    geom_jitter(aes(y=rowSums(otu_table(subset_taxa(PSTSS, family%in%"Eimeriidae")))))+
    xlab("Eimeria OPG")+ ylab("ASV's Eimeriidae rel. abundance")+
    annotate(geom="text", label="Spearman rho=0.49, p<0.001", x=305500, y=0.6, color="gray9")+
    theme_classic()
dev.off()
cor.test(PSTSS@sam_data$Pinworms, rowSums(otu_table(subset_taxa(PSTSS, family%in%"Oxyuridae"))), method ="spearman")
png(filename = "fig/Oocysts_Oxyuris.png",
        width =4, height = 5, units = "in", res= 300)
ggplot(PSTSS@sam_data, aes(x=PSTSS@sam_data$Pinworms)) +
    geom_point(aes(y=rowSums(otu_table(subset_taxa(PSTSS, family%in%"Oxyuridae")))))+
    theme_classic()+
    xlab("Adult pinworms")+ ylab("ASV's Oxyuridae rel. abundance")+
    annotate(geom="text", label="Spearman rho=0.47, p<0.001",x=250, y=0.15, color="gray9")+
    theme_classic()
dev.off()

get_taxa_unique(PS, "family")
PSTSS1@sam_data$Status

PSnem <- subset_taxa(PS, phylum%in%"Nematoda")
PSapi <- subset_taxa(PS, phylum%in%"Apicomplexa")
PSpla <- subset_taxa(PS, phylum%in%"Platyhelminthes")


# 12 families within phylum Nematoda (+1 NA)
relNem <- summarize_taxa(PSnem, "family", "Mouse_ID")
unique(relNem$family)
# 11 families within phyla Apicomplexa (+1 NA)
relApi <- summarize_taxa(PSapi, "family", "Mouse_ID")
unique(relApi$family)
# 2 families within phyla Apicomplexa (+1 NA)
relpla <- summarize_taxa(PSpla, "family", "Mouse_ID")
unique(relpla$family)

## ploting parasite abundace at family level
relNemSorted <- relNem %>%
      mutate(family = fct_reorder(family, meanRA))
relNemSorted %>%
    group_by(family)%>%
    mutate(family_avg=mean(meanRA))%>%
    mutate(Main_family=family_avg>=0.01) -> relNemSorted
relNemSorted %>%
  mutate(family=fct_reorder(family, meanRA, .fun="mean"))%>%
  ggplot(aes(x = reorder(family, meanRA), y = meanRA, color= family)) +
      coord_flip() +
      scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005)) +
      labs(x = NULL, y = "Relative abundance") +
      theme(legend.position = "none")+
      #theme_classic()+
      #geom_hline(aes(yintercept = sample_avg), color = "gray70", size = 0.6) +
      geom_jitter(size = 2, alpha = 0.5, width = 0.2)+
    stat_summary(fun = mean, geom = "point", size = 1, color= "black")+
    theme_minimal()+
    theme(legend.position="none")-> NemAb
png(filename = "fig/Nematoda_ab.png",
        width =4, height = 5, units = "in", res= 400)
NemAb
dev.off()

relApiSorted <- relApi %>%
      mutate(family = fct_reorder(family, meanRA))
relApiSorted %>%
    group_by(family)%>%
    mutate(family_avg=mean(meanRA))%>%
    mutate(Main_family=family_avg>=0.01) -> relApiSorted
relApiSorted %>%
  mutate(family=fct_reorder(family, meanRA, .fun="mean"))%>%
  ggplot(aes(x = reorder(family, meanRA), y = meanRA, color= family)) +
      coord_flip() +
      scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005)) +
      labs(x = NULL, y = "Relative abundance") +
      theme(legend.position = "none")+
    theme_minimal()+
    theme(legend.position="none")+
      #geom_hline(aes(yintercept = sample_avg), color = "gray70", size = 0.6) +
      geom_jitter(size = 2, alpha = 0.5, width = 0.2)+
      stat_summary(fun = mean, geom = "point", size = 1, color= "black") -> ApiAb
png(filename = "fig/Apicomplexa_ab.png",
        width =4, height = 5, units = "in", res= 400)
ApiAb
dev.off()

relplaSorted <- relpla %>%
      mutate(family = fct_reorder(family, meanRA))
relplaSorted %>%
    group_by(family)%>%
    mutate(family_avg=mean(meanRA))%>%
    mutate(Main_family=family_avg>=0.01) -> relplaSorted
relplaSorted %>%
  mutate(family=fct_reorder(family, meanRA, .fun="mean"))%>%
  ggplot(aes(x = reorder(family, meanRA), y = meanRA, color= family)) +
      coord_flip() +
      scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005)) +
      labs(x = NULL, y = "Relative abundance") +
      theme(legend.position = "none")+
      #theme_classic()+
      #geom_hline(aes(yintercept = sample_avg), color = "gray70", size = 0.6) +
      geom_jitter(size = 2, alpha = 0.5, width = 0.2)+
    stat_summary(fun = mean, geom = "point", size = 1, color= "black")+
    theme_minimal()-> plaAb
png(filename = "fig/Platyhelm_ab.png",
        width =4, height = 2, units = "in", res= 400)
plaAb
dev.off()

PSTSS = transform_sample_counts(PS, function(x){x / sum(x)})
PSTSS1 = transform_sample_counts(PS, function(x){(x / sum(x))*10000})


### beta diversity
#I. Bray-Curtis distances

Bra_MDSRA <- ordinate(PSTSS, "NMDS", "bray")

Bra_MDSRA

BC_itri <-plot_ordination(PSTSS, Bra_MDSRA, color = "itriASV")

BC_itri <- ggplot(BC_itri$data, BC_itri$mapping)+ geom_point(size=1, alpha=0.5)+
    labs(x="NMDS1", y="NMDS2", color="Trichuriidae intensity")+
    theme_minimal()

BC_itri

sample_data(PSTSS)$hexASV
psmds <- metaMDS(otu_table(PSTSS))
plot(psmds, display="sites")
names(sample_data(PS))

# goodness of fit for the nmds plot
stressplot(psmds)
gofp=goodness(psmds)
plot(psmds, display="sites", main="OTU goodness of fit")
points(psmds, display="sites", cex=gofp*50, pch = 21, bg="gray")
head(Bra_MDSRA$vectors)

png(filename = "fig/nMDS_itri.png",
        width =6, height = 4, units = "in", res= 400)
BC_itri
dev.off()

head(sdata)

#### permanovas on beta diversisty
PSTSSBMI  <- PSTSS %>%
    subset_samples(BMI != "NA")

PCdis=distance(PSTSS, method="bray")
sdata=sample_data(PSTSS)
permaPSTSS=adonis2(PCdis~sdata$iEimASV+
            sdata$ioxyASV +
            sdata$itriASV+
            sdata$ihexASV +
            sdata$iascASV +
            sdata$icriASV +
            sdata$ihekASV+
            sdata$ihymASV,
        permutations = 10000, method = "bray")

permaPSTSS

sdata$itriASV

anosimPSTSS=anosim(PCdis,sdata$EimASV, permutations=10000, distance="bray")

### testing AKP
Bra_MDSRA <- ordinate(PSTSS, method = "MDS", distance = "bray", na.rm=TRUE)
evals <- Bra_MDSRA$values[,1]

summary(lm(evals~sdata$EimASV))
sdata$dis <- evals
sdata$EimASV <- as.factor(sdata$EimASV)

png(filename = "fig/dis_eim.png",
        width =4, height = 6, units = "in", res= 400)
ggplot(sdata,aes(y=dis, x=EimASV))+
    geom_violin()+
    annotate("text", x= 1.5,y=31, size=14, label="***")+
    geom_jitter(width=0.3, alpha=0.5)+
    labs(y="Bray-distance", x= "Eimeriidae")+
    theme_minimal()
dev.off()

png(filename="fig/dis_bmi.png",
    width=4, height=6, units="in", res=400)
ggplot(sdata,aes(y=dis, x=BMI))+
    geom_jitter(width=0.3, alpha=0.5)+
    labs(y="Bray-distance")+
    annotate("text", x=0, y=20, label="Spearman rho=-0.01, p=0.7")+
    theme_minimal()
dev.off()

summary(sdata$dis)
wilcox.test(sdata$dis~ sdata$EimASV)
cor.test(sdata$dis, sdata$BMI, method="spearman")

wilcox.test(sdata$dis~ sdata$triASV)
wilcox.test(sdata$dis~ sdata$oxyASV)
wilcox.test(sdata$dis~ sdata$hexASV)
wilcox.test(sdata$dis~ sdata$ascASV)
wilcox.test(sdata$dis~ sdata$hymASV)

library(breakaway)
library(DivNet)
#### richness test
library(vegan)


    ## divnet takes forever and then aborts itself when I run at ASV level

#highP <- phyloseq_filter_prevalence(PSgen, prev.trh = 0.95, abund.trh = NULL,
#                                   threshold_condition = "OR", abund.type = "total")
#taxa_names(highP)
#
#divnet_gen <- PSgen %>%
#    divnet(ncores=20, tuning="careful", base="ACCAAGGAGTCTAACGTCTATGCGAGTGTTTGGGTGTAAAACCCGTACGCGGAATGAAAGTGAACGTAGGTTGGGGCCCCCCGGGGTGCACAATCGACCGATCCTGATGTTTTCGGATGGATTTGAGTAGGAGCATAGCTGTTGGGACCCGAAAGATGGTGAACTATGCCTGAATAGGGTGAAGCCAGAGGAAACTCTGGTGGAGGCTCGTAGCGGTTCTGACGTGCAAATCGATCGTCGAATTTGGGTATAGGGGCGAAAGACTAATCGAACCATCT")
#saveRDS(divnet_gen,  file="tmp/divnet_gen.R")
divnetGen=readRDS(file="tmp/divnet_gen.R")

divsha=(divnetGen$shannon)
df <- data.frame(matrix(unlist(divsha), nrow=length(divsha), byrow=T))
df=df[c("X1")]
rownames(df)=names(divsha)
names(df)=c("shannon")
sample_data(PSgen)$shannon = df$shannon

PSgen@sam_data$shannon=as.numeric(PSgen@sam_data$shannon)

ggplot(PSgen@sam_data, aes(BMI, shannon))+
    geom_point()

sort(phyloseq::sample_sums(PS))

pdf("fig/phshanGen.pdf")
plot_richness(PSgen, x="BMI", measures="Shannon")
dev.off()

cor.test(PSgen@sam_data$BMI, PSgen@sam_data$shannon, method="spearman")

# at ASV level - FAIL
#divnet_asv <- PS %>%
#    divnet(ncores=20, tuning="careful", base="TCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTGGCATGCTAACTAGTTACGCGACCCCCGAGCGGTCGGCGTCCCCCAACTTCTTAGAGGGACAAGTGGCGTTCAGCCACCCGAGATTGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCGGGGCTGCACGCGCGCTACACTGACTGGCTCAGCGTGTGCCTACCCTACGCCGGCAGGCGCGGGTAACCCGTTGAACCCCATTCGTGATGGGGATCGGGGATTGCAATTATTCCCCATGAACGAGGAATTCCCAGTAAGTGCGGGTCATAAGCTTGCGTTGATTAAGTCCCTGCC")
#saveRDS(divnet_asv,  file="tmp/divnet_asv.R")

## plotting DA
# Differencial abundance analysis
#remotes::install <- github("FrederickHuangLin/ANCOMBC")
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

#### Calculate normalized stochasticity ration (NST)
library("NST")
tNST(otups, sample_data(PS)$EimASV)

