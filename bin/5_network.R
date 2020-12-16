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
library(SpiecEasi)
library(igraph)

library(taxonomizr)

PS=readRDS(file="tmp/PSpre.R")
PSeim <- subset_taxa(PS, family%in%"Eimeriidae")

PSpara <- subset_taxa(PS, phylum%in%c("Apicomplexa", "Nematoda", "Platyhelminthes"))

# network for Eimeria
#otu.sub <- t(otu_table(PSeim)@.Data)
#tax.sub <- as.data.frame(tax_table(PSeim)@.Data)
# parallel multicores
pargs <- list(rep.num=50, seed=10010, ncores=20)

t1 <- Sys.time()
se.mb.PSeim <- spiec.easi(PSeim, method="mb", pulsar.params=pargs)
t2 <- Sys.time()
t2-t1

ig.mb <- adj2igraph(getRefit(se.mb.PSeim), vertex.attr=list(name=taxa_names(PSeim)))

png(filename = "fig/Netwoork_eim.png",
        width =15, height = 15, units = "in", res= 300)
plot_network(ig.mb, PSeim, type="taxa", label=NULL, color="species")
dev.off()

# network for parasites
t1 <- Sys.time()
se.mb.PSpara <- spiec.easi(PSpara, method="mb", pulsar.params=pargs)
t2 <- Sys.time()
t2-t1

para.ig.mb <- adj2igraph(getRefit(se.mb.PSpara), vertex.attr=list(name=taxa_names(PSpara)))

png(filename = "fig/Netwoork_para.png",
        width =15, height = 15, units = "in", res= 300)
plot_network(para.ig.mb, PSpara, type="taxa", label=NULL, color="family")
dev.off()

# network for overall
t1 <- Sys.time()
se.mb.PS <- spiec.easi(PS, method="mb", pulsar.params=pargs)
t2 <- Sys.time()
t2-t1

saveRDS(se.mb.PS, file="tmp/se.mb.PS.rds")

ps.ig.mb <- adj2igraph(getRefit(se.mb.PS), vertex.attr=list(name=taxa_names(PS)))
png(filename = "fig/Netwoork_all.png",
        width =15, height = 15, units = "in", res= 300)
plot_network(ps.ig.mb, PS, type="taxa", label=NULL)
dev.off()

(se.mb.PS$refit)

otu.ids=taxa_names(PS)

edges <- E(ps.ig.mb)

length(edges)

a=ends(ps.ig.mb, edges[1])

which(otu.ids==a[1])

which(otu.ids==a)

edge.colors=c()
for (e.index in 1:length(edges)){
    adj.nodes=ends(ps.ig.mb, edges[e.index])
    xindex=which(otu.ids==adj.nodes[1])
    yindex=which(otu.ids==adj.nodes[2])
    beta=betaMat[xindex, yindex]
    if (beta>0){
        edge.colors=append(edge.colors, "forestgreen")
    }else if(beta<0){
              edge.colors=append(edge.colors, "red")
    }
}


E(ps.ig.mb)$color=edge.colors

ps.ig.mb.b=ps.ig.mb

nodenames=V(ps.ig.mb.b)$name


#V(spiec.graph.b)$name=getTaxonomy(nodenames, taxa.f, useRownames=TRUE)
library("igraph")



E(ps.ig.mb.b)$arrow.size=5
V(ps.ig.mb.b)$color="white"
V(ps.ig.mb.b)$frame.color="black"

png(filename = "fig/Netwoork_all_test.png",
        width =15, height = 15, units = "in", res= 300)

plot(ps.ig.mb.b)

dev.off()


betaMat=as.matrix(symBeta(getOptBeta(se.mb.PS)))

clusters=cluster_fast_greedy(ps.ig.mb)

clusterOneIndices=which(clusters$membership==1)
clusterOneOtus=clusters$names[clusterOneIndices]
clusterTwoIndices=which(clusters$membership==2)
clusterTwoOtus=clusters$names[clusterTwoIndices]

st1.mb <- degree.distribution(ps.ig.mb)

class(ps.ig.mb)

plot(0:(length(st1.mb)-1), st1.mb)

plot_network(ps.ig.mb, PS, type="taxa", label=NULL, color="family")




# We divide by two since an edge is represented by two entries in the matrix.
positive=length(betaMat[betaMat>0])/2

negative=length(betaMat[betaMat<0])/2

total=length(betaMat[betaMat!=0])/2

se.mb.PS
