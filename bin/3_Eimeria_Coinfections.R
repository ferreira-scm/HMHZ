# Eimeria species

source("bin/Pre_Functions.R")

eim_sp <- psmelt(Eim18_sp)
eim <- psmelt(Eim.TSS18)
eim2 <- psmelt(Eim18)

Eim.m0 <- phyloseq::prune_samples(sample_sums(Eim.TSS.m)>0, Eim.TSS.m)
Eim.m0@tax_table[,6] <- Eim.m0@tax_table[,7]
Eim.m0.sp <- tax_glom(Eim.m0, "species")
Eim.m0 <- phyloseq::prune_samples(sample_sums(Eim.TSS.m)>0, Eim.TSS.m)

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

nb.cols <- 9
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

eim$genus <- as.factor(eim$genus)
eim.m$genus <- as.factor(eim.m$genus)

mycolors2 <- mycolors[c(1, 2, 4:9)]

mycolors2

com_plot_Amp <- ggplot(eim, aes(x=Sample, y=Abundance, fill=genus))+
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
      legend.position="top")
#    coord_flip()

com_plot_Amp

com_plot <- ggplot(eim, aes(x=Sample, y=Abundance, fill=species))+
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

levels(as.factor(eim$species))

Com.m.all <- ggplot(eim.m, aes(x=Sample, y=Abundance, fill=species))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=c("pink", "forestgreen", "dodgerblue4",  "darkgray", "darkred"))+
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

Com.m.all_amp <- ggplot(eim.m, aes(x=Sample, y=Abundance, fill=genus))+
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
      legend.position="top")
#    coord_flip()


Com.m.all_amp

library(cowplot)

Comp_amplicon <- plot_grid(com_plot, com_plot_Amp, ncol=1, align="v", rel_heights=c(1, 1))

ggsave("fig/Eimeria_Amplicon_composition.pdf", Comp_amplicon, height=9, width=14, dpi=400)

ggsave("fig/Eimeria_Amplicon_composition.png", Comp_amplicon, height=9, width=14, dpi=200)

library(viridis)
library(wesanderson)
library(scales)

pal <- wes_palette("Zissou1", 500, type = "continuous")
Eim_heat <-ggplot(eim_sp, aes(Sample, species, fill=Abundance))+
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

Eim_heat

eim.m0

Eim_heat_all <- ggplot(eim.m0, aes(Sample, species, fill=Abundance))+
    geom_tile()+
      labs(y="Eimeria ASV species", x="Sample", fill="Proportion within Eimeria ASVs")+
    scale_fill_gradientn(colours = pal, values=rescale(c(0,0.25,1)),)+
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

ggplot(eim.m1, aes(Sample, genus, fill=Abundance))+
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

levels(as.factor(eim.m$species))

Sp.m <- ggplot(eim.m, aes(x=genus, y=Abundance, fill=species))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=c("pink","forestgreen", "dodgerblue4", "darkgray", "darkred"))+
    labs(fill="Amplicon", x="Sample", y="Proportion within all ASVs")+
    theme_bw(base_size=14)+
    guides(fill=guide_legend(ncol=4))+
    theme(axis.text.y = element_text(colour = 'black', size = 14, face = 'italic'),
      axis.title.x=element_blank(),
#      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text = element_text(colour = 'black', size = 10, face = 'italic'),
      legend.position="none")+
    coord_flip()

Sp.m

M.dis <- plot_grid(plot_grid(Com.m.all, Eim_heat_all, ncol=1, align="v", rel_heights=c(1, 1)), Sp.m, ncol=1, rel_heights=c(1, 0.6))

M.dis

ggsave("fig/Eimeria_18S_28S_amplicons.pdf", M.dis, height=9, width=14, dpi=400)
ggsave("fig/Eimeria_18S_28S_amplicons.png", M.dis, height=9, width=14, dpi=400)


#####################################################################
# Ok, let's try and figure it out what is happening with these co-infections
# removing empty samples

Eim.ra0 <- phyloseq::prune_samples(sample_sums(Eim18_sp)>0, Eim18_sp)

Fer <- subset_taxa(Eim.ra0, species %in% "ferrisi")
sample_data(Eim.ra0)$Ferrisi <- sample_sums(Fer)

Fal <- subset_taxa(Eim.ra0, species %in% "falciformis")
sample_data(Eim.ra0)$Falciformis <- sample_sums(Fal)

Ver <- subset_taxa(Eim.ra0, species %in% "vermiformis")
sample_data(Eim.ra0)$Vermiformis <- sample_sums(Ver)

Sp <- subset_taxa(Eim.ra0, species %in% "sp")
sample_data(Eim.ra0)$Sp <- sample_sums(Sp)

sample_data(Eim.ra0)$Sp

#Eimra0 <- subset_samples(Eim.ra0, !BMI=="NA")

Eimdf <- sample_data(Eim.ra0)
Eimdf$Locality <- as.factor(Eimdf$Locality)
Eimdf <- as.data.frame(Eimdf)
class(Eimdf) <- "data.frame"

dis=phyloseq::distance(Eim.ra0, method="bray", type="samples")

dis

#permaPS=adonis2(dis~
#            Eimdf$HI+
#            Eimdf$Locality+
#            Eimdf$Year,
#            permutations = 1000, method = "bray")

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

Eimdf$sp[Eimdf$Sp>0] <- 1
Eimdf$sp[Eimdf$Sp==0] <- 0


### new variable with amplicon
falModel <- glmer(fal~ver*fer*sp + (1|Locality), family=binomial(), data=Eimdf)

ferModel <- glmer(fer~ver*fal*sp + (1|Locality), family=binomial(), data=Eimdf)

spModel <- glmer(sp~fer*fal*ver + (1|Locality), family=binomial(), data=Eimdf)

verModel <- glmer(ver~fer*fal*sp + (1|Locality), family=binomial(), data=Eimdf)

summary(falModel)
summary(ferModel)
summary(verModel)
summary(spModel)

library("effects")

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
ranova(falModel)

#### quantit

Eimdf$EimeriaTotal <- Eimdf$Falciformis+Eimdf$Vermiformis+Eimdf$Ferrisi + Eimdf$Sp

Eimdf$EimeriaTotal

#FalQ <- lmer(log(1+Falciformis)~log(1+Vermiformis)*log(1+Ferrisi)*log(1+Sp) + (1|Locality), data=Eimdf)

FalQ <- lmer(Falciformis~Vermiformis*Ferrisi*Sp + (1|Locality), data=Eimdf)

summary(FalQ)

FerQ <- lmer(Ferrisi~Vermiformis*Falciformis*Sp + (1|Locality), data=Eimdf)

summary(FerQ)

VerQ <- lmer(Vermiformis~Ferrisi*Falciformis*Sp + (1|Locality), data=Eimdf)

SpQ <- lmer(Sp~Ferrisi*Falciformis*Vermiformis + (1|Locality), data=Eimdf)

summary(VerQ)

summary(SpQ)

library(randomForest)

Fal.fit <- randomForest(Falciformis ~ Ferrisi + Vermiformis + Sp, data=Eimdf, ntree=500, keep.forest=FALSE, importance=TRUE)
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


Fer.fit <- randomForest(Ferrisi ~ Falciformis + Vermiformis + Sp, data=Eimdf, ntree=500, keep.forest=FALSE, importance=TRUE)
Fer.fit
varImpPlot(Fer.fit)

Ver.fit <- randomForest(Vermiformis ~ Falciformis + Ferrisi + Sp, data=Eimdf, ntree=500, keep.forest=FALSE, importance=TRUE)
Ver.fit
varImpPlot(Ver.fit)

Sp.fit <- randomForest(Sp ~ Vermiformis+Falciformis + Ferrisi, data=Eimdf, ntree=500, keep.forest=FALSE, importance=TRUE)
Sp.fit
varImpPlot(Sp.fit)



#FalQ <- glmer.nb(Falciformis~Vermiformis*Ferrisi*Sp + (1|Locality), data=Eimdfq)
#FalQ <- lmer(Falciformis~Vermiformis*Ferrisi*Sp + (1|Locality), data=Eimdfq)

summary(FalQ)

## terrible, even after transforming or using a NB
plot(FalQ)
qqnorm(residuals(FalQ))
qqline(residuals(FalQ))

library(lmerTest)
ranova(FalQ) # not signigicant

FalModel <- glm(Falciformis~Vermiformis*Ferrisi*Sp, data=Eimdf)

plot(FalModel) # not too terrible

FalModel <- glm.nb(Falciformis~Vermiformis*Ferrisi*Sp, data=Eimdfq)
summary(FalModel) # bad deviance too
#plot(FalModel) # ugly! More variables?

eim.rao <- psmelt(Eim.ra0)

summary(as.factor(Eimdf$Concatenated))

head(eim_sp)

eim_sp$Concatenated <- as.factor(eim_sp$Concatenated)

ggplot(eim_sp, aes(x=Concatenated, y=Abundance, fill=species))+
    geom_bar(position="dodge", stat="identity")


library(parasiteLoad)

### Functions
source("bin/TestDistributions.R")

source("bin/bananaplotNoCI.R")

library("fitdistrplus")
library("optimx")
library(FSA)

PS.E <- subset_taxa(fPS, genus%in%"Eimeria")

PS.E@sam_data$EimeriaCounts <- sample_sums(PS.E)

cor.test(PS.E@sam_data$EimeriaCounts, PS.E@sam_data$OPG)

plot(PS.E@sam_data$EimeriaCounts, PS.E@sam_data$OPG)

table(PS.E@sam_data$OPG>0, PS.E@sam_data$EimeriaCounts>0)

summary(PS.E@sam_data$OPG>0)

head(PS.E@sam_data)

PS.E@sam_data$eimeriaSpecies

fPS

#devtools::install_github("alicebalard/parasiteLoad@v2.0")

library(parasiteLoad)

findGoodDist(round(Eimdf$EimeriaTotal*10000),
             distribs = c("normal", "negative binomial", "poisson"),
                          distribs2 = c("norm", "nbinom", "pois"))

Eimdf$EimeriaTotalT <- round(Eimdf$EimeriaTotal*10000)

Eimdf$FerrisiT <- round(Eimdf$Ferrisi*10000)

Eimdf$FalciformisT <- round(Eimdf$Falciformis*10000)

Eimdf$VermiformisT <- round(Eimdf$Vermiformis*10000)

Eimdf$SpT <- round(Eimdf$Sp*10000)

# for eimeria (test)

Eimdf$EimeriaTotalT

#eimraoo <- eimraoo[!is.na(eimraoo$Sex),]

Eimdf[is.na(Eimdf$Sex),]

fitp <- parasiteLoad::analyse(data = Eimdf, response = "EimeriaTotalT", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

fitp

Eimp <- parasiteLoad::bananaPlot(mod = fitp$H0,
                                 data = Eimdf,
                                 response = "EimeriaTotalT",
                                 hybridIndex=seq(0,1,0.005),
                                 islog10 = F,
                                 group = "Sex",
                                 cols = c("#E69F00","#006A4E"))


Eimp

Eimdff <-Eimdf.s[which(Eimdf.s$fer==1),]


fit.fer <- parasiteLoad::analyse(data = Eimdff, response = "FerrisiT", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Eim.fer <- parasiteLoad::bananaPlot(mod = fit.fer$H0,
                                 data = Eimdff,
                                 response = "FerrisiT",
                                 hybridIndex=seq(0,1,0.005),
                                 islog10 = F,
                                 group = "Sex",
                                 cols = c("#E69F00","#006A4E"))


Eim.fer

Eimdffa <-Eimdf.s[which(Eimdf.s$fal==1),]
fit.fal <- parasiteLoad::analyse(data = Eimdffa, response = "FalciformisT", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Eim.fal <- parasiteLoad::bananaPlot(mod = fit.fal$H1,
                                 data = Eimdffa,
                                 response = "FalciformisT",
                                 hybridIndex=seq(0,1,0.005),
                                 islog10 = F,
                                 group = "Sex",
                                 cols = c("#006A4E", "#E69F00"))

Eim.fal

Eimdfv <-Eimdf[which(Eimdf.s$ver==1),]

Eimdfv

fit.ver <- parasiteLoad::analyse(data = Eimdfv, response = "VermiformisT", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Eim.ver <- parasiteLoad::bananaPlot(mod = fit.ver$H1,
                                 data = Eimdfv,
                                 response = "VermiformisT",
                                 hybridIndex=seq(0,1,0.005),
                                 islog10 = F,
                                 group = "Sex",
                                 cols = c("#006A4E", "#E69F00"))

Eim.ver

Eimdfs <-Eimdf[which(Eimdf.s$sp==1),]

fit.sp <- parasiteLoad::analyse(data = Eimdfs, response = "SpT", model = "negbin", group = "Sex", hybridIndex = "HI", myparamBounds="default")

Eim.sp <- parasiteLoad::bananaPlot(mod = fit.sp$H0,
                                 data = Eimdfs,
                                 response = "SpT",
                                 hybridIndex=seq(0,1,0.005),
                                 islog10 = F,
                                 group = "Sex",
                                 cols = c("#006A4E", "#E69F00"))

Eim.sp



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


