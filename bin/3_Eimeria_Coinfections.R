# Eimeria species

source("bin/3_Eim_species_Alignment_Tree_Correlation.R")



######### Preparing phyloseq objects for plotting and so on
# removing empty samples
Eim <- phyloseq::prune_samples(sample_sums(Eim)>0, Eim)
Eim.TSS <- phyloseq::prune_samples(sample_sums(Eim.TSS)>0, Eim.TSS)
Eim.T <- phyloseq::prune_samples(sample_sums(Eim.T)>0, Eim.T)

Eim.T@tax_table[,6] <- "Eimeria"
Eim_sp <- tax_glom(Eim.T, "Species")

amp_names <- gsub("_ASV.*", "", names18S)
Eim@tax_table[,6] <- amp_names
Eim.TSS@tax_table[,6] <- amp_names
Eim.T@tax_table[,6] <- amp_names

# separating Eimeria ASV's by gene
Eim18 <- Eim
Eim.TSS18 <- Eim.TSS
Eim.T18 <- Eim.T
Eim28 <- Eim
Eim.TSS28 <- Eim.TSS
Eim.T28 <- Eim.T

Eim18 <- subset_taxa(Eim18, !Genus=="D3A_5Mod_46_F.D3B_5Mod_46_R")
Eim.TSS18 <- subset_taxa(Eim.TSS18, !Genus=="D3A_5Mod_46_F.D3B_5Mod_46_R")
Eim.T18 <- subset_taxa(Eim.T18, !Genus=="D3A_5Mod_46_F.D3B_5Mod_46_R")
Eim28 <- subset_taxa(Eim28, Genus=="D3A_5Mod_46_F.D3B_5Mod_46_R")
Eim.TSS28 <- subset_taxa(Eim.TSS28, Genus=="D3A_5Mod_46_F.D3B_5Mod_46_R")
Eim.T28 <- subset_taxa(Eim.T28, Genus=="D3A_5Mod_46_F.D3B_5Mod_46_R")

#Eim28@tax_table[10,6]

get_taxa_unique(Eim18, "Species")
get_taxa_unique(Eim28, "Species")

Eim18 <- phyloseq::prune_samples(sample_sums(Eim18)>0, Eim18)
Eim.TSS18 <- phyloseq::prune_samples(sample_sums(Eim.TSS18)>0, Eim.TSS18)
Eim.T18 <- phyloseq::prune_samples(sample_sums(Eim.T18)>0, Eim.T18)
#eim18.TSS <- transform_sample_counts(Eim18, function(x) x / sum(x)) 

#sanity check
colnames(Eim18_sp@otu_table)==rownames(Eim18_sp@tax_table)
#colnames(Eim18_sp@otu_table) <- Eim18_sp@tax_table[,5]

eim_sp <- psmelt(Eim18_sp)
eim <- psmelt(Eim.T18)
eim2 <- psmelt(Eim18)

# specieal for plotting, Eimeria ASV proportion within all Eimeria ASVs
Eim.m1 <- phyloseq::prune_samples(sample_sums(Eim)>0, Eim)
Eim.m1 <- tax_glom(Eim.m1, "Genus")
Eim.m1 <- transform_sample_counts(Eim.m1, function(x) x / sum(x)) 

eim.m1 <- psmelt(Eim.m1)
eim.m0 <- psmelt(Eim_sp)
eim.m <- psmelt(Eim.T)
eim$Genus <- as.factor(eim$Genus)

# relevel
dist_bc18 <- (vegdist(Eim.T18@otu_table, method="bray"))
dist_bc182 <- (vegdist(Eim18_sp@otu_table, method="bray"))
dist_bc <- (vegdist(Eim.T@otu_table, method="bray"))
dist_bc2 <- (vegdist(Eim_sp@otu_table, method="bray"))
res18 <- pcoa(dist_bc18)
res18.2 <- pcoa(dist_bc182)
res <- pcoa(dist_bc)
res2 <- pcoa(dist_bc2)

res

#plot_ordination(Eim18, ordinate(Eim18, "MDS"))
#plot_ordination(Eim.TSS.m, ordinate(Eim.TSS.m, "MDS"))
#plot_ordination(Eim.TSS.m, ordinate(Eim.TSS.m, "MDS"))

EH_sort18 <- names(sort(res18$vectors[,1]))
EH_sort182 <- names(sort(res18.2$vectors[,1]))
EH_sort <- names(sort(res$vectors[,1]))
EH_sort2 <- names(sort(res2$vectors[,1]))

eim$Sample <- factor(eim$Sample, levels= EH_sort18)
eim_sp$Sample <- factor(eim_sp$Sample, levels=EH_sort182)
eim.m$Sample <- factor(eim.m$Sample, levels= EH_sort)
eim.m0$Sample <- factor(eim.m0$Sample, levels= EH_sort) # I will still sort with the distances of all amplicons, so that I can align plots
eim.m1$Sample <- factor(eim.m1$Sample, levels= EH_sort)

EH_sort18==levels(eim$Sample)
EH_sort182==levels(eim_sp$Sample)
EH_sort==levels(eim.m$Sample)

eim$Genus <- as.factor(eim$Genus)
eim.m$Genus <- as.factor(eim.m$Genus)

levels(eim$Genus)
levels(eim.m$Genus)

nb.cols <- 9+1
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
mycolors2 <- mycolors[c(1, 2, 4:10)]


#### plotting 18S amplicons
com_plot_Amp <- ggplot(eim, aes(x=Sample, y=Abundance, fill=Genus))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=mycolors2)+
    labs(fill="Amplicon", x="Sample", y="Proportion within all ASVs/ngDNA")+
    theme_bw(base_size=14)+
    theme(axis.text.y = element_text(colour = 'black', size = 14, face = 'italic'),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text = element_text(colour = 'black', size = 14, face = 'italic'),
      legend.position="top")
#    coord_flip()

com_plot_Amp

####  now same but fill with species
com_plot <- ggplot(eim, aes(x=Sample, y=Abundance, fill=Species))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    labs(fill="Eimeria ASV species", x="Sample", y="Proportion within all ASVs/ngDNA")+
    theme_bw(base_size=14)+
    theme(axis.text.y = element_text(colour = 'black', size = 14, face = 'italic'),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text = element_text(colour = 'black', size = 14, face = 'italic'),
      legend.position="bottom")
#    coord_flip()

com_plot

#### litle cheat for plotting
eim.m$Species2 <- eim.m$Species

eim.m$Species2[eim.m$Genus=="D3A_5Mod_46_F.D3B_5Mod_46_R"] <- paste(eim.m$Species2[eim.m$Genus=="D3A_5Mod_46_F.D3B_5Mod_46_R"], "28S", sep="_")


################# now plotting both 18S and 28S genes
Com.m.all <- ggplot(eim.m, aes(x=Sample, y=Abundance, fill=Species))+
    geom_bar(position="stack", stat="identity")+
#    scale_fill_manual(values=c("forestgreen", "pink", "dodgerblue4",  "darkgray", "darkred"))+
    scale_fill_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    labs(fill="Eimeria", x="Sample", y="Proportion within all ASVs")+
    theme_bw(base_size=14)+
    theme(axis.text.y = element_text(colour = 'black', size = 14, face = 'italic'),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text = element_text(colour = 'black', size = 14, face = 'italic'),
      legend.position="top")
#    coord_flip()

Com.m.all

Com.m.all_amp <- ggplot(eim.m, aes(x=Sample, y=Abundance, fill=Genus))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=mycolors)+
    labs(fill="Amplicon", x="Sample", y="Proportion within all ASVs/ng DNA")+
    theme_bw(base_size=14)+
    theme(axis.text.y = element_text(colour = 'black', size = 14, face = 'italic'),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text = element_text(colour = 'black', size = 14, face = 'italic'),
      legend.position="top")
#    coord_flip()


Com.m.all_amp

library(cowplot)

Comp_amplicon <- plot_grid(com_plot, com_plot_Amp, ncol=1, align="v", rel_heights=c(1, 1))

Comp_amplicon2 <- plot_grid(Com.m.all, Com.m.all_amp, ncol=1, align="v", rel_heights=c(1, 1))

Comp_amplicon2

ggsave("fig/Eimeria_Amplicon_composition.pdf", Comp_amplicon2, height=9, width=14, dpi=400)
ggsave("fig/Eimeria_Amplicon_composition.png", Comp_amplicon2, height=9, width=14, dpi=200)

library(viridis)
library(wesanderson)
library(scales)

pal <- wes_palette("Zissou1", 500, type = "continuous")

Eim_heat_all <- ggplot(eim.m0, aes(Sample, Species, fill=Abundance))+
    geom_tile()+
    labs(y="Eimeria species", x="Sample", fill="Proportion within all ASVs/ng DNA")+
    scale_fill_gradientn(colours = pal, values=rescale(c(0,0.5,1)),)+
    theme_bw(base_size=14)+
      theme(axis.text.y = element_text(colour = 'black', size = 14, face = 'italic'),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text=element_text(size=10),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="bottom")

Eim_heat_all

Amplicon_ASV <- ggplot(eim.m1, aes(Sample, Genus, fill=Abundance))+
    geom_tile()+
      labs(y="Eimeria ASV species", x="Sample", fill="Proportion within Eimeria ASVs")+
    scale_fill_gradientn(colours = pal, values=rescale(c(0,0.5,1)),)+
    theme_bw(base_size=14)+
      theme(axis.text.y = element_text(colour = 'black', size = 14, face = 'italic'),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text=element_text(size=10),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="bottom")

Amplicon_ASV

Sp.m <-ggplot(eim.m[eim.m$Abundance>0,], aes(x=Genus, y=Abundance, fill=Species))+
#    geom_bar(position="stack", stat="identity")+
    geom_point(size=4, colour="gray", shape=21, position=position_jitterdodge(dodge.width=0.8, jitter.width=0.15), alpha=0.4)+
    geom_boxplot(alpha=0.3, colour="black", outlier.shape = NA)+
    scale_fill_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    labs(fill="Amplicon", x="", y="Sum of proportion within all ASVs/ng DNA")+
    theme_bw(base_size=10)+
    guides(fill=guide_legend(ncol=4))+
    theme(axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
       legend.key = element_blank(),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text = element_text(colour = 'black', size = 10, face = 'italic'),
      legend.position="none"
      )+
    coord_flip()

Sp.m

# sanity check
levels(Eim_heat_all$data$Sample)==levels(Com.m.all_amp$data$Sample)

M.dis <- plot_grid(plot_grid(Eim_heat_all, Com.m.all_amp, ncol=1, align="v", rel_heights=c(1, 1)), Sp.m, ncol=1, rel_heights=c(1, 0.6))

M.dis

ggsave("fig/Eimeria_18S_28S_amplicons.pdf", M.dis, height=12, width=14, dpi=400)
ggsave("fig/Eimeria_18S_28S_amplicons.png", M.dis, height=12, width=14, dpi=400)

ggsave("fig/Eimeria_amplicon_sp.pdf", Sp.m, height=4, width=10, dpi=400)
ggsave("fig/Eimeria_amplicon_sp.png", Sp.m, height=4, width=10, dpi=400)


ggsave("fig/Eimeria_amplicon_heatplot.pdf", Amplicon_ASV, height=4, width=12, dpi=400)
ggsave("fig/Eimeria_amplicon_heatplot.png", Amplicon_ASV, height=4, width=12, dpi=400)

### OPG and species abundance
#eim.m0

Eim_sp@sam_data$Eimeira_ASVs <- sample_sums(Eim_sp)
df.opg <- Eim_sp@sam_data
class(df.opg) <- "data.frame"

## first OPG and Eimeria spp
df.opg0 <- df.opg[df.opg$OPG>0,]
df.opg0 <- df.opg0[!is.na(df.opg0$OPG),]
cor.test(df.opg$OPG, df.opg$Eimeira_ASVs)
cor.test(df.opg0$OPG, df.opg0$Eimeira_ASVs)
cor.test(log(df.opg0$OPG), log(df.opg0$Eimeira_ASVs))

ggplot(df.opg0, aes(x=log(OPG), y=log(Eimeira_ASVs)))+
#    geom_bar(position="stack", stat="identity")+
    geom_point()

# now OPG and Eimeria sp (each species separately)
cor.test(eim.m0$Abundance[eim.m0$Species=="ferrisi"],eim.m0$OPG[eim.m0$Species=="ferrisi"])
cor.test(eim.m0$Abundance[eim.m0$Species=="falciformis"],eim.m0$OPG[eim.m0$Species=="falciformis"])
cor.test(eim.m0$Abundance[eim.m0$Species=="vermiformis"],eim.m0$OPG[eim.m0$Species=="vermiformis"])

eim.opg <- eim.m0[eim.m0$OPG>0,]
eim.opg <- eim.opg[eim.opg$Abundance>0,]

eim.opg <- eim.opg[!is.na(eim.opg$OPG),]
eim.opg <- eim.opg[!is.na(eim.opg$Abundance),]

cor.test(eim.opg$Abundance[eim.opg$Species=="ferrisi"],eim.opg$OPG[eim.opg$Species=="ferrisi"])

cor.test(eim.opg$Abundance[eim.opg$Species=="falciformis"],eim.opg$OPG[eim.opg$Species=="falciformis"])

cor.test(eim.opg$Abundance[eim.opg$Species=="vermiformis"],eim.opg$OPG[eim.opg$Species=="vermiformis"])


Eim.OPG <- ggplot(eim.opg, aes(y=log(OPG), x=log(Abundance), fill=Species))+
    geom_point(shape=21, size=4, alpha=0.7)+
    scale_fill_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    scale_colour_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    labs(fill="Eimeria", y="Oocyst/g faeces (log)", x="Proportion within all ASVs/ng DNA (log)")+
    annotate(geom="text", x=-3, y=17, label="Ferrisi: Pearson's rho=0.25, p=0.09\nFalciformis: rho=0.80, p<.001\nVermiformis: rho=-0.06, p=0.93", size=2)+ 
    #geom_smooth(method=lm, aes(colour=Species))+
    theme_bw(base_size=10)+
    guides(fill=guide_legend(nrow=1))+
    theme(axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
#      axis.title.x=element_blank(),
      legend.key = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(colour = 'black', size = 10, face = 'italic'),
      legend.position="top"
      )

Eim.OPG

table(eim.m0$Abundance[eim.m0$Species=="ferrisi"]>0, eim.m0$OPG[eim.m0$Species=="ferrisi"]>0)
table(eim.m0$Abundance[eim.m0$Species=="falciformis"]>0, eim.m0$OPG[eim.m0$Species=="falciformis"]>0)
table(eim.m0$Abundance[eim.m0$Species=="vermiformis"]>0, eim.m0$OPG[eim.m0$Species=="vermiformis"]>0)

ggsave("fig/Eimeria_OPG.pdf", Eim.OPG, height=4, width=4, dpi=400)
ggsave("fig/Eimeria_OPG.png", Eim.OPG, height=4, width=4, dpi=400)


#####################################################################
# Ok, let's try and figure it out what is happening with these co-infections
# removing empty samples



Fer <- subset_taxa(Eim_sp, Species %in% "ferrisi")
sample_data(Eim_sp)$Ferrisi <- sample_sums(Fer)

Fal <- subset_taxa(Eim_sp, Species %in% "falciformis")
sample_data(Eim_sp)$Falciformis <- sample_sums(Fal)

Ver <- subset_taxa(Eim_sp, Species %in% "vermiformis")
sample_data(Eim_sp)$Vermiformis <- sample_sums(Ver)

sample_data(Eim_sp)$Ferrisi

    
#Eimra0 <- subset_samples(Eim.ra0, !BMI=="NA")

Eimdf <- sample_data(Eim_sp)
Eimdf$Locality <- as.factor(Eimdf$Locality)
Eimdf <- as.data.frame(Eimdf)
class(Eimdf) <- "data.frame"

dis=phyloseq::distance(Eim.ra0, method="bray", type="samples")

dis

permaPS=adonis2(dis~
            Eimdf$HI+
            Eimdf$Locality+
            Eimdf$Year,
            permutations = 1000, method = "bray")

permaPS



#chisq.test(table(Eimdf$Ferrisi>0, Eimdf$Falciformis>0))
#chisq.test(table(Eimdf$Ferrisi>0, Eimdf$Vermiformis>0))
#chisq.test(table(Eimdf$Ferrisi>0, Eimdf$Sp>0))
#chisq.test(table(Eimdf$Falciformis>0, Eimdf$Vermiformis>0))

#table(Eimdf$Ferrisi>0, Eimdf$Falciformis>0)
#table(Eimdf$Ferrisi>0, Eimdf$Vermiformis>0)
#table(Eimdf$Falciformis>0, Eimdf$Vermiformis>0)

library(lme4)

Eimdf$fal[Eimdf$Falciformis>0] <- 1
Eimdf$fal[Eimdf$Falciformis==0] <- 0

Eimdf$fer[Eimdf$Ferrisi>0] <- 1
Eimdf$fer[Eimdf$Ferrisi==0] <- 0

Eimdf$ver[Eimdf$Vermiformis>0] <- 1
Eimdf$ver[Eimdf$Vermiformis==0] <- 0

### new variable with amplicon
falModel <- glmer(fal~ver*fer + (1|Locality), family=binomial(), data=Eimdf)

ferModel <- glmer(fer~ver*fal + (1|Locality), family=binomial(), data=Eimdf)

verModel <- glmer(ver~fer*fal + (1|Locality), family=binomial(), data=Eimdf)

summary(falModel)

summary(ferModel)

summary(verModel)

#library("effects")

#F.plot <- plot(effect("ver", ferModel),
#     ylab="Probability of E. falciformis",
#     xlab="E. vermiformis infection",
#     main="")

#png("fig/Vermiformis_effect.png", height=4, width=4, units="in", res=400)
#V.plot
#dev.off()

#png("fig/Falciformis_effect.png", height=4, width=4, units="in", res=400)
#Int
#dev.off()

#Int <- plot(effect("ver:fer", falModel, xlevels=list(fer=0:1)),
#     ylab="Probability E. falciformis +",
#     multiline=TRUE,
#     main="E.falciformis* E.ferrisi")

library(lmerTest)
#ranova(falModel)

#### quantit

Eimdf$EimeriaTotal <- Eimdf$Falciformis+Eimdf$Vermiformis+Eimdf$Ferrisi

Eimdf$EimeriaTotal

#FalQ <- lmer(log(1+Falciformis)~log(1+Vermiformis)*log(1+Ferrisi)*log(1+Sp) + (1|Locality), data=Eimdf)

FalQ <- lmer(Falciformis~Vermiformis*Ferrisi + (1|Locality), data=Eimdf)

summary(FalQ)

FerQ <- lmer(Ferrisi~Vermiformis*Falciformis + (1|Locality), data=Eimdf)

summary(FerQ)

VerQ <- lmer(Vermiformis~Ferrisi*Falciformis + (1|Locality), data=Eimdf)

summary(VerQ)

library(randomForest)

Fal.fit <- randomForest(Falciformis ~ Ferrisi + Vermiformis, data=Eimdf, ntree=500, keep.forest=FALSE, importance=TRUE)
Fal.fit
varImpPlot(Fal.fit)


ImpData <- as.data.frame(importance(Fal.fit))
ImpData$Var.Names <- row.names(ImpData)

ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`))+
    geom_segment(aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
    geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
    theme_light() +
    coord_flip() +
    theme(
        legend.position="bottom",
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()
                      )


Fer.fit <- randomForest(Ferrisi ~ Falciformis + Vermiformis , data=Eimdf, ntree=500, keep.forest=FALSE, importance=TRUE)
Fer.fit
varImpPlot(Fer.fit)

Ver.fit <- randomForest(Vermiformis ~ Falciformis + Ferrisi, data=Eimdf, ntree=500, keep.forest=FALSE, importance=TRUE)
Ver.fit
varImpPlot(Ver.fit)

#FalQ <- glmer.nb(Falciformis~Vermiformis*Ferrisi*Sp + (1|Locality), data=Eimdfq)
#FalQ <- lmer(Falciformis~Vermiformis*Ferrisi*Sp + (1|Locality), data=Eimdfq)

## terrible, even after transforming or using a NB
#plot(FalQ)
#qqnorm(residuals(FalQ))
#qqline(residuals(FalQ))

library(lmerTest)
ranova(FalQ) # not signigicant

#FalModel <- glm(Falciformis~Vermiformis*Ferrisi, data=Eimdf)

#plot(FalModel) # not too terrible

#FalModel <- glm.nb(Falciformis~Vermiformis*Ferrisi*Sp, data=Eimdfq)
summary(FalModel) # bad deviance too
#plot(FalModel) # ugly! More variables?

eim.rao <- psmelt(Eim.ra0)

summary(as.factor(Eimdf$Concatenated))

head(eim_sp)

eim_sp$Concatenated <- as.factor(eim_sp$Concatenated)


eim_c <- eim_sp[eim_sp$Abundance>0,]
eim_c <- eim_c[!is.na(eim_c$Concatenated),]


length(unique(as.factor(eim_c$Mouse_ID)))


eim_c[eim_c$Concatenated=="E_vermiformis", c("Abundance", "Mouse_ID", "Group_18S", "Group_COI_1", "Group_COI_2", "Group_ORF", "Concatenated", "Species")]


TISSUE_MA <- ggplot(eim_c, aes(x=Concatenated, y=Abundance, fill=Species))+
    scale_fill_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    geom_boxplot(alpha=0.3, colour="black", outlier.shape = NA)+
    geom_point(size=4, shape=21, position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), alpha=0.7)+
    labs(x="genotyping PCR from tissue DNA", y="Proportion within all ASVs/ng DNA (log)")+    
#    geom_bar(position="dodge", stat="identity")+
    guides(fill=guide_legend(ncol=4))+
    theme_classic()+
    theme(axis.text.y = element_text(colour = 'black', size = 14, face = 'italic'),
      axis.title.x=element_blank(),
#      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text = element_text(colour = 'black', size = 10, face = 'italic'),
      legend.position="none"
      )

ggsave("fig/Eimeria_qPCR_MA.pdf", TISSUE_MA, height=4, width=5, dpi=400)
ggsave("fig/Eimeria_qPCR_MA.png", TISSUE_MA, height=4, width=5, dpi=400)

Eimdf

library(parasiteLoad)

### Functions
source("bin/TestDistributions.R")

source("bin/bananaplotNoCI.R")

library("fitdistrplus")
library("optimx")
library(FSA)

Eim.T

PS.E <- subset_taxa(fPS, Genus%in%"g__Eimeria")

PS.E@sam_data$EimeriaCounts <- sample_sums(PS.E)

cor.test(PS.E@sam_data$EimeriaCounts, PS.E@sam_data$OPG)

plot(PS.E@sam_data$EimeriaCounts, PS.E@sam_data$OPG)

table(PS.E@sam_data$OPG>0, PS.E@sam_data$EimeriaCounts>0)

summary(PS.E@sam_data$OPG>0)

#devtools::install_github("alicebalard/parasiteLoad@v2.0")

findGoodDist(round(Eimdf$EimeriaTotal*10000),
             distribs = c("normal", "negative binomial", "poisson"),
                          distribs2 = c("norm", "nbinom", "pois"))

Eimdf$EimeriaTotalT <- round(Eimdf$EimeriaTotal*10000)

Eimdf$FerrisiT <- round(Eimdf$Ferrisi*10000)

Eimdf$FalciformisT <- round(Eimdf$Falciformis*10000)

Eimdf$VermiformisT <- round(Eimdf$Vermiformis*10000)

# for eimeria (test)

Eimdf$EimeriaTotalT

#eimraoo <- eimraoo[!is.na(eimraoo$Sex),]

fitp <- parasiteLoad::analyse(data = Eimdf, response = "EimeriaTotalT", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

fitp

Eimp <- parasiteLoad::bananaPlot(mod = fitp$H1,
                                 data = Eimdf,
                                 response = "EimeriaTotalT",
                                 hybridIndex=seq(0,1,0.005),
                                 islog10 = F,
                                 group = "Sex",
                                 cols = c("#E69F00","#006A4E"))


Eimp

Eimdff <-Eimdf[which(Eimdf$fer==1),]


fit.fer <- parasiteLoad::analyse(data = Eimdf, response = "FerrisiT", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Eim.fer <- parasiteLoad::bananaPlot(mod = fit.fer$H2,
                                 data = Eimdf,
                                 response = "FerrisiT",
                                 hybridIndex=seq(0,1,0.005),
                                 islog10 = F,
                                 group = "Sex",
                                 cols = c("#E69F00","#006A4E"))


Eim.fer

Eimdffa <-Eimdf[which(Eimdf$fal==1),]

fit.fal <- parasiteLoad::analyse(data = Eimdf, response = "FalciformisT", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Eim.fal <- parasiteLoad::bananaPlot(mod = fit.fal$H1,
                                 data = Eimdffa,
                                 response = "FalciformisT",
                                 hybridIndex=seq(0,1,0.005),
                                 islog10 = F,
                                 group = "Sex",
                                 cols = c("#006A4E", "#E69F00"))

Eim.fal

Eimdfv <-Eimdf[which(Eimdf$ver==1),]

fit.ver <- parasiteLoad::analyse(data = Eimdfv, response = "VermiformisT", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Eim.ver <- parasiteLoad::bananaPlot(mod = fit.ver$H0,
                                 data = Eimdfv,
                                 response = "VermiformisT",
                                 hybridIndex=seq(0,1,0.005),
                                 islog10 = F,
                                 group = "Sex",
                                 cols = c("#006A4E", "#E69F00"))

Eim.ver



library(SpiecEasi)
library(igraph)

# parallel multicores
pargs <- list(rep.num=1000, seed=10010, ncores=90, thresh=0.05)
## mb
t1 <- Sys.time()
netEim <- spiec.easi(Eim18, method="mb", pulsar.params=pargs)
t2 <- Sys.time()
t2-t1


#t1 <- Sys.time()
#netEimgl <- spiec.easi(Eim18, method="glasso", pulsar.params=pargs)
#t2 <- Sys.time()
#t2-t1

saveRDS(netEim, "tmp/netEim.RDS")
#saveRDS(netEimgl, "tmp/netEimgl.RDS")
                                        #saveRDS(se.gl10net, "tmp/se.gl10net.RDS")

netEim <- readRDS("tmp/netEim.RDS")
#netEimgl <- readRDS("tmp/netEimgl.RDS")


                                        # looking at lambda path
netEim$select$stars$summary

# coding nodes
eim.id=Eim.m18@tax_table[,5]

#### now we improve our visualization
## I want weighted edges

bm=symBeta(getOptBeta(netEim), mode="maxabs")
diag(bm) <- 0
                                        #weights <- Matrix::summary(t(bm))[,3] # includes negative weights
weights <- (1-Matrix::summary(t(bm))[,3])/2 # or

bm

weights

netE <- adj2igraph(Matrix::drop0(getRefit(netEim)),
                    edge.attr=list(weight=weights),
                    vertex.attr = list(name=eim.id))

betaMat=as.matrix(symBeta(getOptBeta(netEim)))

# we want positive edges to be green and negative to be red
edges <- E(netE)
edge.colors=c()

for (e.index in 1:length(edges)){
    adj.nodes=ends(netE, edges[e.index])
    xindex=which(eim.id==adj.nodes[1])
    yindex=which(eim.id==adj.nodes[2])
    beta=betaMat[xindex, yindex]
    if (beta>0){
        edge.colors=append(edge.colors, "#1B7837")
    }else if(beta<0){
        edge.colors=append(edge.colors, "#762A83")
    }
}

E(netE)$color=edge.colors 

V(netE)$species=Eim18@tax_table[,7]
V(netE)$amplicon<- Eim.m18@tax_table[,6]

E(netE)$weight

# we also want the node color to code for family
nb.col <- length(levels(as.factor(V(netE)$amplicon)))
coul <- colorRampPalette(brewer.pal(8, "Accent"))(nb.col)

coul

plot(netE,
#     layout=layout <- with <- fr(net10b),
     vertex.label=V(netE)$species,
#     vertex.size=as.integer(cut(hub.sb, breaks=10))+2,
     vertex.color=adjustcolor(coul,0.8),
     edge.width=as.integer(cut(E(netE)$weight, breaks=6))/3,
     margin=c(0,1,0,0))
legend(x=-2, y=1, legend=levels(as.factor(V(netE)$amplicon)), col=coul, bty="n",x.intersp=0.25,text.width=0.045, pch=20, pt.cex=1.5)

###### with spring

devtools::install_github("stefpeschel/NetCoMi",
            dependencies = c("Depends", "Imports", "LinkingTo"),
            repos = c("https://cloud.r-project.org/",
            BiocManager::repositories()))

library(NetCoMi)
library(mixedCCA)

net.s <- netConstruct(Eim18@otu_table,
                              measure = "spring",
                              measurePar = list(nlambda=10,
                                                rep.num=10),
                              normMethod = "none",
                              zeroMethod = "none",
                              sparsMethod = "none",
                              dissFunc = "signed",
                              verbose = 3,
                              seed = 123456)

net.s

net.is <- igraph::graph_from_adjacency_matrix(net.s$adjaMat1, weighted = TRUE, diag=FALSE)

#net.is <- adj2igraph(net.s$adjaMat1, diag=FALSE)

net.is

E(net.is)$weight

net.s$adjaMat[(net.s$adjaMat1>0&net.s$adjaMat1<1)]

net.s$adjaMat1

E(net.is)$weight==net.s$adjaMat[(net.s$adjaMat1>0&net.s$adjaMat1<1)]

e.col <- NULL
assM <- net.s$assoMat1[(net.s$adjaMat1>0&net.s$adjaMat1<1)]


net.s$assoMat1[(net.s$adjaMat1>0&net.s$adjaMat1<1)]

net.s$assoMat1[which(net.s$adjaMat1%in%E(net.is)$weight)]

(net.s$adjaMat1%in%E(net.is)$weight)

net.s$adjaMat1 %in% E(net.is)$weight


for (i in seq(length(assM))){
    if (net.s$assoMat1[(net.s$adjaMat1>0)][i]>0){
        e.col[i] <-  "#1B7837"}
    else if (net.s$assoMat1[(net.s$adjaMat1>0)][i]<0){
        e.col[i] <-   "#762A83"}
}


V(net.is)

E(net.is)$color=e.col 

e.col


net.is <- as.undirected(net.is)

V(net.is)$species=Eim18@tax_table[,5]
V(net.is)$amplicon<- Eim.m18@tax_table[,6]

# we also want the node color to code for family
nb <- length(levels(as.factor(V(net.is)$amplicon)))
col <- colorRampPalette(brewer.pal(8, "Accent"))(nb.col)

col

E(net.is)$color

plot(net.is,
#     layout=layout <- with <- fr(net10b),
     vertex.label=V(net.is)$species,
#     vertex.size=as.integer(cut(hub.sb, breaks=10))+2,
     vertex.color=adjustcolor(col,0.8),
#     edge.width=as.integer(cut(E(net.is)$weight, breaks=6)),
     margin=c(0,1,0,0))
legend(x=-2, y=1, legend=levels(as.factor(V(net.is)$amplicon)), col=col, bty="n",x.intersp=0.25,text.width=0.045, pch=20, pt.cex=1.5)





props <- netAnalyze(net.s,
                              centrLCC = TRUE,
                              clustMethod = "cluster_fast_greedy",
                              hubPar = "eigenvector",
                                                         weightDeg = FALSE, normDeg = FALSE)



p <- plot(props,
          nodeColor = "cluster",
          nodeSize = "eigenvector",
          title1 = "Network on OTU level with SPRING associations",
          showTitle = TRUE,
          cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"),
              bty = "n", horiz = TRUE)


