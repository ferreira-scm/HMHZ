# plotting community

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

library(metagMisc)

PS=readRDS(file="tmp/PSpre.R")
pPS=readRDS(file="tmp/pPSpre.R")
get_taxa_unique(pPS, "superkingdom")

## prevalence plotting
prevdf = apply(X = otu_table(pPS),
               MARGIN = ifelse(taxa_are_rows(pPS), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame

prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(pPS),
                    tax_table(pPS))
plyr::ddply(prevdf, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Subset to the remaining phyla

prevdf1 = subset(prevdf, phylum %in% get_taxa_unique(pPS, "phylum"))

pdf("fig/prevPhyla.pdf")

ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(pPS),color=phylum)) +
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
relphy=summarize_taxa(pPS, "phylum", "Mouse_ID")
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

png(filename = "fig/BacAbundance_phyla.png",
        width =7, height = 5, units = "in", res= 300)
phyAb
dev.off()

## prevalence plotting
prevdf = apply(X = otu_table(pPS),
               MARGIN = ifelse(taxa_are_rows(pPS), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame

prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(pPS),
                    tax_table(pPS))
plyr::ddply(prevdf, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Subset to the remaining phyla

prevdf1 = subset(prevdf, phylum %in% get_taxa_unique(pPS, "phylum"))

pdf("fig/prevPhyla.pdf")
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(pPS),color=phylum)) +
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
relphy=summarize_taxa(bPS, "phylum", "Mouse_ID")
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

png(filename = "fig/BacAbundance_phyla.png",
        width =7, height = 5, units = "in", res= 300)
phyAb
dev.off()



get_taxa_unique(PS, "family")
PSTSS1@sam_data$Status

PSnem <- subset_taxa(pPS, phylum%in%"Nematoda")
PSapi <- subset_taxa(pPS, phylum%in%"Apicomplexa")
PSpla <- subset_taxa(pPS, phylum%in%"Platyhelminthes")

# 7 families within phylum Nematoda (+1 NA)
relNem <- summarize_taxa(PSnem, "family", "Mouse_ID")
unique(relNem$family)
# 5 families within phyla Apicomplexa (+1 NA)
relApi <- summarize_taxa(PSapi, "family", "Mouse_ID")
unique(relApi$family)
# 1 families within phyla Apicomplexa
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

## Testing Eimeriidae correlation with OPG
cor.test(PSrare@sam_data$OPG_Eimeria, rowSums(otu_table(subset_taxa(PSrare, family%in%"Eimeriidae"))), method ="spearman")

cor.test(PSrare@sam_data$OPG_Eimeria, sample_data(PSrare)$eimASV, method ="spearman")

png(filename = "fig/Oocysts_Eimeriidae-rare.png",
        width =4, height = 5, units = "in", res= 300)

ggplot(PSrare@sam_data, aes(x=PSrare@sam_data$OPG_Eimeria)) +
    geom_jitter(aes(y=sample_data(PSrare)$eimASV))+
    xlab("Eimeria OPG")+ ylab("ASV's Eimeriidae abundance - rarefied")+
    annotate(geom="text", label="Spearman rho=0.48, p<0.001", x=305500, y=0.6, color="gray9")+
    theme_classic()

dev.off()


### Average relative OTU abundances per domain
euPS <- subset_taxa(PS, superkingdom%in%"Eukaryota")
baPS<- subset_taxa(PS, superkingdom%in%"Bacteria")
ePSphy <- tax_glom(euPS, "phylum")

# abundance filtering to 1%? Or keep prevalence filtering?
x = taxa_sums(euPSphy)
keepTaxa = (x / sum(x) > 0.01)

KT <- which(keepTaxa==FALSE)

KT

euPS0 = merge_taxa(ePSphy, KT, 1)



euPS1 = prune_taxa(keepTaxa, euPS)

ePSphy <- tax_glom(euPS, "phylum")

ePSphy1 <- tax_glom(euPS1, "phylum")

get_taxa_unique(ePSphy, "phylum")

get_taxa_unique(ePSphy1, "phylum")

sample_data(ePSphy)$all <- "one"

ePS1 <- merge_samples(ePSphy, "all")

ePS0 <- transform_sample_counts(ePS1, function(x) x / sum(x))

head(otu_table(ePS0))

plot_bar(ePS0, fill="phylum")

##start remove
sample_data(ePSphy1)$all <- "one"

ePS11 <- merge_samples(ePSphy1, "all")

ePS01 <- transform_sample_counts(ePS11, function(x) x / sum(x))

head(otu_table(ePS01))

plot_bar(ePS01, fill="phylum")

##### end
ps2 <- transform <- sample <- counts(ps1, function(x) x / sum(x))
plot <- bar(ps2, fill="Phylum")

ePS <- subset_taxa(pPS, superkingdom%in%"Eukaryota")
bPS<- subset_taxa(pPS, superkingdom%in%"Bacteria")

nb.cols <- length(get_taxa_unique(ePS, "phylum"))
coul <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
ePS1 <- ePS %>%
        aggregate_taxa(level = "phylum") %>%
                      microbiome::transform(transform = "compositional")
ePS1 %>%
    plot_composition(sample.sort ="neatmap", otu.sort = "abundance") +
    scale_fill_manual(values = coul) +
    labs(fill="Phyla")+
    scale_y_continuous(label = scales::percent)+
    theme_classic()+
    theme(axis.text.x=element_blank()) -> eukaryaComp

png(filename = "fig/EukAbundance_phyla.png",
        width =7, height = 5, units = "in", res= 300)
eukaryaComp
dev.off()

eukaryaComp

nb.cols <- length(get_taxa_unique(bPS, "phylum"))
coul <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
bPS1 <- bPS %>%
        aggregate_taxa(level = "phylum") %>%
                      microbiome::transform(transform = "compositional")
bPS1 %>%
    plot_composition(sample.sort ="neatmap", otu.sort = "abundance") +
    scale_fill_manual(values = coul) +
    labs(fill="Phyla")+
    scale_y_continuous(label = scales::percent)+
    theme_classic()+
    theme(axis.text.x=element_blank()) -> bacComp

png(filename = "fig/BacAbundance_phyla.png",
        width =7, height = 5, units = "in", res= 300)
bacComp
dev.off()



eukaryaComp

