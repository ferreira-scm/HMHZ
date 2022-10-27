# high level analysis of the "wild microbiome"

source("bin/Pre_Functions.R")
# high level analysis of the "wild microbiome"

eim_sp <- psmelt(Eim18_sp)
eim <- psmelt(Eim.TSS18)
eim2 <- psmelt(Eim18)

Eim.m_sp

Eim.m0 <- phyloseq::prune_samples(sample_sums(Eim.TSS.m)>0, Eim.TSS.m)

Eim.m0.sp <- tax_glom(Eim.m0, "genus")

Eim.m1 <- phyloseq::prune_samples(sample_sums(Eim.m)>0, Eim.m)
Eim.m1 <- tax_glom(Eim.m1, "genus")
Eim.m1 <- transform_sample_counts(Eim.m1, function(x) x / sum(x)) 

eim.m1 <- psmelt(Eim.m1)
eim.m0 <- psmelt(Eim.m0.sp)
eim.m <- psmelt(Eim.m0)
eim$genus <- as.factor(eim$genus)

# relevel
dist_bc18 <- (vegdist(Eim.TSS18@otu_table, method="bray"))
dist_bc182 <- (vegdist(Eim18_sp@otu_table, method="bray"))
dist_bc <- (vegdist(Eim.m0@otu_table, method="bray"))
res18 <- pcoa(dist_bc18)
res18.2 <- pcoa(dist_bc182)
res <- pcoa(dist_bc)

#plot_ordination(Eim18, ordinate(Eim18, "MDS"))
#plot_ordination(Eim.TSS.m, ordinate(Eim.TSS.m, "MDS"))
#plot_ordination(Eim.TSS.m, ordinate(Eim.TSS.m, "MDS"))

EH_sort18 <- names(sort(res18$vectors[,1]))
EH_sort182 <- names(sort(res18.2$vectors[,1]))
EH_sort <- names(sort(res$vectors[,1]))

eim$Sample <- factor(eim$Sample, levels= EH_sort18)
eim_sp$Sample <- factor(eim_sp$Sample, levels=EH_sort182)
eim.m$Sample <- factor(eim.m$Sample, levels= EH_sort)
eim.m0$Sample <- factor(eim.m0$Sample, levels= EH_sort)
eim.m1$Sample <- factor(eim.m1$Sample, levels= EH_sort)

EH_sort18==levels(eim$Sample)
EH_sort182==levels(eim_sp$Sample)
EH_sort==levels(eim.m$Sample)

nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

eim$species <- as.factor(eim$species)

eim.m$species <- as.factor(eim.m$species)

mycolors2 <- mycolors[c(1, 2, 4:9, 11:19)]
mycolors2

com_plot_Amp <- ggplot(eim, aes(x=Sample, y=Abundance, fill=species))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=mycolors2)+
    labs(fill="Amplicon", x="Sample", y="Proportion within all ASVs")+
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
      legend.position="none")
#    coord_flip()

com_plot <- ggplot(eim, aes(x=Sample, y=Abundance, fill=genus))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=c("forestgreen", "dodgerblue4",  "darkgray", "darkred"))+
    labs(fill="Eimeria ASV species", x="Sample", y="Proportion within all ASVs")+
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

Com.m.all <- ggplot(eim.m, aes(x=Sample, y=Abundance, fill=species))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=mycolors)+
    labs(fill="Amplicon", x="Sample", y="Proportion within all ASVs")+
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
      legend.position="none")
#    coord_flip()

Com.m.all

library(cowplot)

Comp_amplicon <- plot_grid(com_plot, com_plot_Amp, ncol=1, align="v", rel_heights=c(1, 1))
ggsave("fig/Eimeria_Amplicon_composition.pdf", Comp_amplicon, height=12, width=10, dpi=400)
ggsave("fig/Eimeria_Amplicon_composition.png", Comp_amplicon, height=12, width=10, dpi=400)



library(viridis)
library(wesanderson)
library(scales)

pal <- wes_palette("Zissou1", 500, type = "continuous")
Eim_heat <-ggplot(eim_sp, aes(Sample, genus, fill=Abundance))+
    geom_tile()+
      labs(y="Eimeria ASV species", x="Sample", fill="Proportion within all ASVs")+
    scale_fill_gradientn(colours = pal, values=rescale(c(0, 0.25, 1)))+
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
#    scale_fill_distiller(palette = "RdPu")
#    scale_fill_gradient2(low="#00AFBB", mid="#E7B800", high="#FC4E07", midpoint=0.5)
#    scale_fill_viridis(discrete=FALSE)

Eim_heat_all <- ggplot(eim.m1, aes(Sample, genus, fill=Abundance))+
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

ggplot(eim.m0, aes(Sample, genus, fill=Abundance))+
    geom_tile()+
      labs(y="Eimeria ASV species", x="Sample", fill="Proportion within Eimeria ASVs")+
    scale_fill_gradientn(colours = pal, values=rescale(c(0,0.05,1)),)+
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

                                        #    scale_fill_distiller(palette = "RdPu")
#    scale_fill_gradient2(low="#00AFBB", mid="#E7B800", high="#FC4E07", midpoint=0.5)
#    scale_fill_viridis(discrete=FALSE)

eim.m$Abundance[eim.m$genus=="28S"]

Eim_heat

Eim_heat_all

summary(eim.m$Abundance)

Sp.m <- ggplot(eim.m, aes(x=genus, y=Abundance, fill=species))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=mycolors)+
    labs(fill="Amplicon", x="Sample", y="Proportion within all ASVs")+
    theme_bw(base_size=14)+
    theme(axis.text.y = element_text(colour = 'black', size = 14, face = 'italic'),
      axis.title.x=element_blank(),
#      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text = element_text(colour = 'black', size = 14, face = 'italic'),
      legend.position="bottom")+
    coord_flip()

M.dis <- plot_grid(plot_grid(Com.m.all, Eim_heat_all, ncol=1, align="v", rel_heights=c(1, 1)), Sp.m, ncol=1)

ggsave("fig/Eimeria_18S_28S_amplicons.pdf", M.dis, height=10, width=12, dpi=400)
ggsave("fig/Eimeria_18S_28S_amplicons.png", M.dis, height=10, width=12, dpi=400)


#####
EimRA <- subset_taxa(fPS.TSS, genus%in%"Eimeria")

EimRA <- subset_samples(EimRA, Mouse_ID%in%EH_sort182)

EimRA <- tax_glom(EimRA, "genus")

# sanity check
levels(as.factor(eimRA$Mouse_ID))==EH_sort182 # ups, we need to reorder
eimRA <- psmelt(EimRA)

eimRA$Mouse_ID <- factor(eimRA$Mouse_ID, levels=EH_sort182)

RA_EIM <- ggplot(eimRA, aes(x=Mouse_ID, y=Abundance))+
    geom_point(shape=21, colour="gray", fill="#c51b8a", size=3, alpha=0.8)+
    geom_bar(position="stack", stat="identity", fill="#dd1c77", alpha=0.2)+
#    scale_fill_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    labs(x="Sample", y="Proportion within total ASVs")+
    scale_y_continuous(limits=c(0,0.3))+
    theme_bw(base_size=14)+
    theme(axis.text.y = element_text(colour = 'black', size = 14, face = 'italic'),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill="white"))

RA_EIM

library(cowplot)

Eim_species_comp <- plot_grid(com_plot, Eim_heat, RA_EIM, ncol=1, align="v", rel_heights=c(0.9,0.6,0.6))

Eim_species_comp

Eim_species_comp23 <- plot_grid(Eim_heat, RA_EIM, ncol=1, align="v", rel_heights=c(1, 1))

Eim_species_comp23

ggsave("fig/Eimeria_heat_composition.pdf", Eim_species_comp, height=9, width=10, dpi=400)
ggsave("fig/Eimeria_heat_composition.png", Eim_species_comp, height=9, width=10, dpi=400)

ggsave("fig/Eimeria_heat_composition23.pdf", Eim_species_comp23, height=6, width=8, dpi=600)
ggsave("fig/Eimeria_heat_composition23.png", Eim_species_comp23, height=6, width=8, dpi=600)

ggsave("fig/Eimeria_heat_composition1.png", com_plot, height=4, width=8, dpi=600)
ggsave("fig/Eimeria_heat_composition2.png", Eim_heat, height=4, width=8, dpi=600)
ggsave("fig/Eimeria_heat_composition3.png", RA_EIM, height=4, width=8, dpi=600)


ggsave("fig/Eimeria_heat_composition1.pdf", com_plot, height=4, width=8, dpi=600)
ggsave("fig/Eimeria_heat_composition2.pdf", Eim_heat, height=4, width=8, dpi=600)
ggsave("fig/Eimeria_heat_composition3.pdf", RA_EIM, height=4, width=8, dpi=600)

ggsave("fig/Eimeria_heat_composition_species.pdf", Eim_heat_sp, height=6, width=10, dpi=600)

phyloseq::plot_bar(eim.TSS, fill="genus")+
    geom_bar(stat="identity", position="stack")+
    scale_fill_manual(values=c("white", "darkgray", "forestgreen", "dodgerblue4", "deepskyblue", "darksalmon", "darkred"))


#sample_data(Eim.TSS18)$Mouse_ID <- factor(sample_data(Eim.TSS18)$Mouse_ID, levels=EH_sort182)
#levels(as.factor(sample_data(Eim.TSS18)$Mouse_ID))==EH_sort182
#####################################################################
# Ok, let's try and figure it out what is happening with these co-infections
# removing empty samples

Eim.TSS18

Eim.ra0 <- phyloseq::prune_samples(sample_sums(Eim.TSS18)>0, Eim.TSS18)

Eim.ra0 <- tax_glom(Eim.TSS18, "genus")


Fer <- subset_taxa(Eim.ra0, genus %in% "ferrisi")
sample_data(Eim.ra0)$Ferrisi <- sample_sums(Fer)

Fal <- subset_taxa(Eim.ra0, genus %in% "falciformis")
sample_data(Eim.ra0)$Falciformis <- sample_sums(Fal)

Ver <- subset_taxa(Eim.ra0, genus %in% "vermiformis")
sample_data(Eim.ra0)$Vermiformis <- sample_sums(Ver)

Sp <- subset_taxa(Eim.ra0, genus %in% "sp.")
sample_data(Eim.ra0)$Sp <- sample_sums(Sp)

sample_data(Eim.ra0)$Sp

#Eimra0 <- subset_samples(Eim.ra0, !BMI=="NA")

sdata=sample_data(Eim.ra0)
dis=phyloseq::distance((Eim.ra0), method="bray", type="samples")

permaPS=adonis2(dis~
            sdata$HI+
            sdata$Locality+
            sdata$Year,
        permutations = 1000, method = "bray")

permaPS


Eimdf <- (sample_data(Eim.ra0))

chisq.test(table(Eimdf$Ferrisi>0, Eimdf$Falciformis>0))

chisq.test(table(Eimdf$Ferrisi>0, Eimdf$Vermiformis>0))

chisq.test(table(Eimdf$Ferrisi>0, Eimdf$Sp>0))

chisq.test(table(Eimdf$Falciformis>0, Eimdf$Vermiformis>0))

table(Eimdf$Ferrisi>0, Eimdf$Falciformis>0)

table(Eimdf$Ferrisi>0, Eimdf$Vermiformis>0)

table(Eimdf$Falciformis>0, Eimdf$Vermiformis>0)

Eimdf$Locality <- as.factor(Eimdf$Locality)

Eimdf <- as.data.frame(Eimdf)

class(Eimdf) <- "data.frame"

library(lme4)

Eimdf$fal[Eimdf$Falciformis>0] <- 1
Eimdf$fal[Eimdf$Falciformis==0] <- 0

Eimdf$fer[Eimdf$Ferrisi>0] <- 1
Eimdf$fer[Eimdf$Ferrisi==0] <- 0

Eimdf$ver[Eimdf$Vermiformis>0] <- 1
Eimdf$ver[Eimdf$Vermiformis==0] <- 0

### new variable with amplicon

Eimdf$qPCRsummary

FalModel <- glmer(fal~ver*fer + (1|Locality), family=binomial(), data=Eimdf)

summary(FalModel)

summary(FerModel)

summary(VerModel)

FerModel <- glmer(fer~ver*fal*Sp + (1|Locality), family=binomial(), data=Eimdf)

Eimdf$Sp

SpModel <- glmer(Sp~fer*fal*ver + (1|Locality), family=binomial(), data=Eimdf)

summary(SpModel)

VerModel <- glmer(ver~fer*fal + (1|Locality), family=binomial(), data=Eimdf)

sink("fig/Falciformis_infection.txt")
summary(FalModel)
sink()

library("effects")

V.plot <- plot(effect("ver", FalModel),
     ylab="Probability of E. falciformis",
     xlab="E. vermiformis infection",
     main="")

F.plot <- plot(effect("ver", FerModel),
     ylab="Probability of E. falciformis",
     xlab="E. vermiformis infection",
     main="")


F.plot <- plot(effect("ver", FerModel),
     ylab="Probability of E. falciformis",
     xlab="E. vermiformis infection",
     main="")


V.plot

F.plot

png("fig/Vermiformis_effect.png", height=4, width=4, units="in", res=400)
V.plot
dev.off()

png("fig/Falciformis_effect.png", height=4, width=4, units="in", res=400)
Int
dev.off()


Int <- plot(effect("ver:fer", FalModel, xlevels=list(fer=0:1)),
     ylab="Probability E. falciformis +",
     multiline=TRUE,
     main="E.falciformis* E.ferrisi")


summary(VerModel)

anova(FalModel)

library(lmerTest)

ranova(FalModel)

#### quantit

Eimdf$EimeriaTotal <- Eimdf$Falciformis+Eimdf$Vermiformis+Eimdf$Ferrisi + Eimdf$Sp

Eimdfq <- Eimdf[Eimdf$EimeriaTotal>0,]

nrow(Eimdfq)



FalModel <- glmer(log(1+Falciformis)~log(1+Vermiformis)*log(1+Ferrisi) + (1|Locality), data=Eimdfq)

ranova(FalModel) # not signigicant

FalModel <- glm(Falciformis~Vermiformis*Ferrisi, data=Eimdfq)

FalModelI <- glm(log(1+Falciformis)~log(1+Vermiformis)+log(1+Ferrisi), data=Eimdfq)
FalModelV <- glm(log(1+Falciformis)~log(1+Ferrisi), data=Eimdfq)
FalModelF <- glm(log(1+Falciformis)~log(1+Vermiformis), data=Eimdfq)

summary(FalModel)

plot(FalModel) # baaaad, let's try something else

FalModel <- glm(log(1+Falciformis)~log(1+Vermiformis)*log(1+Ferrisi), data=Eimdfq)
plot(FalModel) # also not good. different distribution?

FalModel <- glm.nb(Falciformis~Vermiformis*Ferrisi, data=Eimdfq)
summary(FalModel) # bad deviance too
plot(FalModel) # ugly! More variables?

FalModel <- glm(log(1+Falciformis)~Vermiformis*Ferrisi, data=Eimdfq)

summary(FalModel) # bad deviance too

plot(FalModel) # ugly! More variables?

names(Eimdf)

table(Eimdfq$qPCRsummary)

Eimdf$qPCRsummary

anova(FalModel)

anova(FalModel, FalModelI)
anova(FalModel, FalModelF)


FerModel <- glmer(I(Ferrisi>0)~I(Vermiformis>0)*I(Falciformis>0) + (1|Locality), family=binomial(), data=Eimdf)

VerModel <- glmer(I(Vermiformis>0)~I(Ferrisi>0)*I(Falciformis>0) + (1|Locality), family=binomial(), data=Eimdf)










Eim.ra0

eim.rao <- psmelt(Eim.ra0)

head(eim.rao)

eim.rao$fer[eim.rao$Ferrisi>0] <- 1
eim.rao$fer[eim.rao$Ferrisi==0] <- 0
eim.rao$fer

eim.rao$fal[eim.rao$Falciformis>0] <- 1
eim.rao$fal[eim.rao$Falciformis==0] <- 0
eim.rao$fal

eim.rao$ver[eim.rao$Vermiformis>0] <- 1
eim.rao$ver[eim.rao$Vermiformis==0] <- 0
eim.rao$ver

eim.rao$EimRich <- eim.rao$fer+eim.rao$ver + eim.rao$fal

eim.rao$EimRich

eim.rao$Falciformis

richModel <- glm(Falciformis~Year + Vermiformis*Ferrisi, data=eim.rao)

summary(richModel)

library(parasiteLoad)

### Functions
source("bin/TestDistributions.R")

source("bin/bananaplotNoCI.R")

library("fitdistrplus")
library("optimx")
library(FSA)

#devtools::install_github("alicebalard/parasiteLoad@v2.0")

library(parasiteLoad)

findGoodDist(eimraoo$Abundance,
             distribs = c("normal", "negative binomial", "poisson"),
                          distribs2 = c("norm", "nbinom", "pois"))


# for eimeria (test)
eimraoo <-eim.rao[eim.rao$Abundance>0,]

eimraoo <- eimraoo[!is.na(eimraoo$Sex),]

# gotta round this
eimraoo$Abundance <- round(eimraoo$Abundance*100000)

fitp <- parasiteLoad::analyse(data = eimraoo, response = "Abundance", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

fitp$H3


Eimp <- parasiteLoad::bananaPlot(mod = fitp$H1,
                                 data = eimraoo,
                                 response = "Abundance",
                                 hybridIndex=seq(0,1,0.005),
                                 islog10 = F,
                                 group = "Sex",
                                 cols = c("#006A4E", "#E69F00"))




Reads.plot <- bananaplotNoCI(mod = fitp$H3,
                                       data = eimraoo,
                                       response = "Abundance",
                                       islog10 = F,
                                       group = "Sex",
                                       cols = c("#006A4E", "#006A4E"))+
    labs(tag = "A)")+
    xlab("Mouse genotype (Hybrid Index)")+
    ylab("Total read counts")+
      theme(legend.position ="none")
