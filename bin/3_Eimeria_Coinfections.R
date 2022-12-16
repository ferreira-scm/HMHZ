# Eimeria species
source("bin/3_Eim_species_Alignment_Tree_Correlation.R")

library(cowplot)

######### Preparing phyloseq objects for plotting and so on
# removing empty samples
Eim.Tw@tax_table[,6] <- "Eimeria"
Eim_sp <- tax_glom(Eim.Tw, "Species")

amp_names <- gsub("_ASV.*", "", names18S)
Eim.TSSw@tax_table[,6] <- amp_names
Eim.Tw@tax_table[,6] <- amp_names

# separating Eimeria ASV's by gene
Eim.T18 <- Eim.Tw
Eim.T28 <- Eim.Tw

Eim.T18 <- subset_taxa(Eim.T18, !Genus=="D3A_5Mod_46_F.D3B_5Mod_46_R")
Eim.T28 <- subset_taxa(Eim.T28, Genus=="D3A_5Mod_46_F.D3B_5Mod_46_R")

get_taxa_unique(Eim18, "Species")
get_taxa_unique(Eim28, "Species")

Eim.T18 <- phyloseq::prune_samples(sample_sums(Eim.T18)>0, Eim.T18)
#eim18.TSS <- transform_sample_counts(Eim18, function(x) x / sum(x)) 

eim.m0 <- psmelt(Eim_sp)
eim.m <- psmelt(Eim.Tw)
eim$Genus <- as.factor(eim$Genus)

# relevel
dist_bc <- (vegdist(Eim.Tw@otu_table, method="bray"))
dist_bc2 <- (vegdist(Eim_sp@otu_table, method="bray"))
res <- pcoa(dist_bc)
res2 <- pcoa(dist_bc2)

plot_ordination(Eim.Tw, ordinate(Eim.Tw, "MDS"))

EH_sort <- names(sort(res$vectors[,1]))
EH_sort2 <- names(sort(res2$vectors[,1]))

eim.m$Sample <- factor(eim.m$Sample, levels= EH_sort)
eim.m0$Sample <- factor(eim.m0$Sample, levels= EH_sort) # I will still sort with the distances of all amplicons, so that I can align plots

# sanity check
EH_sort==levels(eim.m$Sample)

eim.m$Genus <- as.factor(eim.m$Genus)

levels(eim.m$Genus)

nb.cols <- 11+1
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

################# now plotting both 18S and 28S genes
Com.m.all <- ggplot(eim.m, aes(x=Sample, y=Abundance, fill=Species))+
    geom_bar(position="stack", stat="identity")+
#    scale_fill_manual(values=c("forestgreen", "pink", "dodgerblue4",  "darkgray", "darkred"))+
    scale_fill_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    labs(fill="Eimeria", x="Sample", y="Eimeria ASV abundance", tag="a")+
    theme_bw(base_size=12)+
    theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
      panel.grid.major = element_blank(),
      legend.position="none")
#    coord_flip()

Com.m.all_amp <- ggplot(eim.m, aes(x=Sample, y=Abundance, fill=Genus))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=mycolors)+
    labs(fill="Amplicon", x="Sample", y="Eimeria ASV abundance", tag="b")+
    theme_bw(base_size=12)+
    theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
      panel.grid.major = element_blank(),
      legend.position="none")

legend <- get_legend(Com.m.all_amp+
                     guides(fill=guide_legend(override.aes=list(size=6),nrow=4))+
                     theme(legend.text = element_text(size = 9),
      legend.position="bottom"))

legend2 <- get_legend(Com.m.all+
                     guides(fill=guide_legend(override.aes=list(size=4), nrow=1))+
                     theme(legend.text = element_text(face = 'italic'),
                           legend.title=element_text(face="italic"),
                           legend.position="top"))



Comp_amplicon2 <- plot_grid(legend2, Com.m.all, Com.m.all_amp, legend,  ncol=1, rel_heights=c(0.04,0.8, 0.8, 0.25))
Comp_amplicon2

ggsave("fig/S3_Eimeria_amplicon_species_distribution.pdf", Comp_amplicon2, height=8, width=10, dpi=400)
ggsave("fig/S3_Eimeria_amplicon_species_distribution.png", Comp_amplicon2, height=9, width=10, dpi=400)

library(viridis)
library(wesanderson)
library(scales)

pal <- wes_palette("Zissou1", 500, type = "continuous")
Eim_heat_all <- ggplot(eim.m0, aes(Sample, Species, fill=Abundance))+
    geom_tile()+
    labs(y="Eimeria", x="Sample", fill="Eimeria ASV abundance")+
    scale_fill_gradientn(colours = pal)+
    theme_bw(base_size=12)+
    theme(axis.text.y = element_text(face = 'italic'),
      axis.title.y=element_text(face="italic"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text=element_text(size=10),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="bottom")

Eim_heat_all

## something seems weird in this plot
Amplicon_ASV <- ggplot(eim.m, aes(Sample, Genus, fill=Abundance))+
    geom_tile()+
      labs(y="Amplicon", x="Sample", fill="Eimeria ASV abundance")+
    scale_fill_gradientn(colours = pal, values=rescale(c(0,0.0005,1)),)+
    theme_bw(base_size=12)+
      theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text=element_text(size=10),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="bottom")

Sp.m <-ggplot(eim.m[eim.m$Abundance>0,], aes(x=Genus, y=Abundance, fill=Species))+
#    geom_bar(position="stack", stat="identity")+
    geom_point(size=4, colour="gray", shape=21, position=position_jitterdodge(dodge.width=0.8, jitter.width=0.15), alpha=0.4)+
    geom_boxplot(alpha=0.3, colour="black", outlier.shape = NA)+
    scale_fill_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    labs(fill="Eimeria", x="", y="Eimeria abundance")+
    theme_bw(base_size=12)+
    guides(fill=guide_legend(ncol=4))+
    theme(axis.text.y = element_text(colour = 'black', size = 10),
       legend.key = element_blank(),
      panel.grid.major = element_blank(),
#      panel.grid.minor = element_blank(),
      legend.text = element_text(face = 'italic'),
      legend.title=element_text(face="italic"),
      legend.position="top")+
    coord_flip()

Sp.m

# sanity check
levels(Eim_heat_all$data$Sample)==levels(Com.m.all_amp$data$Sample)

ggsave("fig/Eimeria_amplicon_sp.pdf", Sp.m, height=4, width=10, dpi=400)
ggsave("fig/Eimeria_amplicon_sp.png", Sp.m, height=4, width=10, dpi=400)

### OPG and species abundance
#eim.m0

Eim_sp@sam_data$Eimeira_ASVs <- sample_sums(Eim_sp)
df.opg <- Eim_sp@sam_data
class(df.opg) <- "data.frame"

## first OPG and Eimeria spp
df.opg0 <- df.opg[df.opg$OPG>0,]
df.opg0 <- df.opg0[!is.na(df.opg0$OPG),]

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

cor.test(log(eim.opg$Abundance[eim.opg$Species=="ferrisi"]),log(eim.opg$OPG[eim.opg$Species=="ferrisi"]))

cor.test(log(eim.opg$Abundance[eim.opg$Species=="falciformis"]),log(eim.opg$OPG[eim.opg$Species=="falciformis"]))

cor.test(log(eim.opg$Abundance[eim.opg$Species=="vermiformis"]),log(eim.opg$OPG[eim.opg$Species=="vermiformis"]))


Eim.OPG <- ggplot(eim.opg, aes(y=log(OPG), x=log(Abundance), fill=Species))+
    geom_point(shape=21, size=4, alpha=0.7)+
    scale_fill_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    scale_colour_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    labs(fill="Eimeria", y="Oocyst/g faeces (log)", x="Eimeria ASV abundance (log)")+
    annotate(geom="text", x=min(log(eim.opg$Abundance)), y=max(log(eim.opg$OPG))-0.2,hjust=0.05, label="Eimeria spp. Pearson's rho=0.47, p<0.001\nFerrisi: Pearson's rho=0.36, p=0.02\nFalciformis: rho=0.40, p=0.04\nVermiformis: rho=0.33, p=0.59", size=2)+ 
    #geom_smooth(method=lm, aes(colour=Species))+
    theme_bw(base_size=10)+
    guides(fill=guide_legend(nrow=1))+
    theme(legend.key = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(face = 'italic'),
      legend.title=element_text(face="italic"),
      legend.position="none")

Eim.OPG

OPG_ab <- readRDS("/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/OPG_Abundance_MA_panel.R")


legend <- get_legend(Eim.OPG+
          guides(fill=guide_legend(nrow=1))+
          theme(legend.text = element_text(face = 'italic'),
          legend.title=element_text(face="italic"),
          legend.position="top"))



Eim.OPG2 <- plot_grid(legend, Eim.OPG, rel_heights=c(0.1,0.8), nrow=2, labels=c("", "c"))
Eim.OPG2

Eim.OPG_p <- plot_grid(OPG_ab, Eim.OPG2, rel_widths=c(1, 0.5))

Eim.OPG_p


ggsave("fig/Eimeria_OPG.pdf", Eim.OPG_p, height=4, width=10, dpi=400)
ggsave("fig/Eimeria_OPG.png", Eim.OPG_p, height=4, width=10, dpi=400)


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

dis2 <- phyloseq::distance(prune_samples(rownames(Eimdf[!is.na(Eimdf$BMI),]), Eim_sp), method="bray", type="samples")

Eimdf1 <- Eimdf[!is.na(Eimdf$BMI),]


permaPS=adonis2(dis2~
#            Eimdf1$HI+
            Eimdf1$Locality+
            Eimdf1$Year+
            Eimdf1$Sex+
           Eimdf1$BMI,
            permutations = 1000, method = "bray")

permaPS

library(merTools)
library(MuMIn)

Eimdf1$logFer <- log(1+Eimdf1$Ferrisi)
Eimdf1$logVer <- log(1+Eimdf1$Vermiformis)
Eimdf1$logFal <- log(1+Eimdf1$Falciformis)

BMIm <- lmer(BMI~logFer * logFal * logVer +(1|Locality), data=Eimdf1)

BMIm0 <- lmer(BMI~1 + (1|Locality), data=Eimdf1)

#BMIm <- lmer(BMI~Ferrisi * Vermiformis * Falciformis + (1|Locality), data=Eimdf1)#, very similar
#summary(BMIm)

summary(BMIm)

plot(BMIm)

qqnorm(resid(BMIm))
qqline(resid(BMIm))

plotREsim(REsim(BMIm))

r.squaredGLMM(BMIm)

anova(BMIm, BMIm0)

ranova(BMIm)

library(sjPlot) #for plotting lmer and glmer mods

library(effects)
F_effect <- as.data.frame(effects::effect(term="logFer", mod=BMIm))
Fal_effect <- as.data.frame(effects::effect(term="logFal", mod=BMIm))
Fer_plot <- ggplot()+
    geom_point(data=Eimdf1, aes(logFer, BMI), shape=21, size=2)+
#    geom_point(data=F_effect, aes(x=logFer, y=fit), fill="dodgerblue4", shape=21, size=4)+
    geom_line(data=F_effect, aes(x=logFer, y=fit), colour="dodgerblue4")+
    geom_ribbon(data=F_effect, aes(x=logFer, ymin=lower, ymax=upper), alpha=0.3, fill="dodgerblue4")+
    theme_bw(base_size=10)+
    labs(x= "Eimeria ferrisi abundance, log(1+)", y="Body mass index")+
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
Fal_plot <- ggplot()+
    geom_point(data=Eimdf1, aes(logFal, BMI), shape=21, size=2)+
#    geom_point(data=Fal_effect, aes(x=logFal, y=fit), fill="forestgreen", shape=21, size=4)+
    geom_line(data=Fal_effect, aes(x=logFal, y=fit), colour="forestgreen")+
    geom_ribbon(data=Fal_effect, aes(x=logFal, ymin=lower, ymax=upper), alpha=0.3, fill="forestgreen")+
        theme_bw(base_size=10)+
    labs(x= "Eimeria falciformis abundance, log(1+)", y="Body mass index")+
    theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

Fig6 <- cowplot::plot_grid(Fer_plot, Fal_plot)

ggsave("fig/Fig6_BMI_Eimeria.pdf", Fig6, height=4, width=9, dpi=400)
ggsave("fig/Fig6_BMI_Eimeria.png", Fig6, height=4, width=9, dpi=400)

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

library(parasiteLoad)

### Functions
source("bin/TestDistributions.R")

source("bin/bananaplotNoCI.R")

library("fitdistrplus")
library("optimx")
library(FSA)

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

Eimp <- parasiteLoad::bananaPlot(mod = fitp$H3,
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
                                 data = Eimdff,
                                 response = "FerrisiT",
                                 hybridIndex=seq(0,1,0.005),
                                 islog10 = F,
                                 group = "Sex",
                                 cols = c("#E69F00","#006A4E"))


Eim.fer

Eimdffa <-Eimdf[which(Eimdf$fal==1),]

fit.fal <- parasiteLoad::analyse(data = Eimdf, response = "FalciformisT", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Eim.fal <- parasiteLoad::bananaPlot(mod = fit.fal$H1,
                                 data = Eimdf,
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



