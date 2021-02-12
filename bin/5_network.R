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
library(RColorBrewer)
library(ggpubr)

PS=readRDS(file="tmp/PSpre.R")
pPS=readRDS(file="tmp/pPSpre.R")

pPSeim <- subset_taxa(pPS, family%in%"Eimeriidae")
pPSpara <- subset_taxa(pPS, phylum%in%c("Apicomplexa", "Nematoda", "Platyhelminthes"))
pPSparagen <- tax_glom(pPSpara, taxrank="genus")
pPSgen <- tax_glom(pPS, taxrank="genus")

####################### network for Eimeria
#otu.sub <- t(otu_table(PSeim)@.Data)
#tax.sub <- as.data.frame(tax_table(PSeim)@.Data)
# parallel multicores
pargs <- list(rep.num=1000, seed=10010, ncores=20)

t1 <- Sys.time()
se.mb.PSeim <- spiec.easi(pPSeim, method="mb", pulsar.params=pargs)
t2 <- Sys.time()
t2-t1

ig.mb <- adj2igraph(getRefit(se.mb.PSeim), vertex.attr=list(name=taxa_names(pPSeim)))

png(filename = "fig/Network_eim.png",
        width =15, height = 15, units = "in", res= 300)
plot_network(ig.mb, pPSeim, type="taxa", label=NULL, color="species")
dev.off()

################### network for parasites
t1 <- Sys.time()
se.mb.PSpara <- spiec.easi(pPSpara, method="mb", pulsar.params=pargs)
t2 <- Sys.time()
t2-t1

para.ig.mb <- adj2igraph(getRefit(se.mb.PSpara), vertex.attr=list(name=taxa_names(pPSpara)))

betaMat=as.matrix(symBeta(getOptBeta(se.mb.PSpara)))

otu.ids=taxa_names(pPSpara)
edges <- E(para.ig.mb)
edge.colors=c()
for (e.index in 1:length(edges)){
    adj.nodes=ends(para.ig.mb, edges[e.index])
    xindex=which(otu.ids==adj.nodes[1])
    yindex=which(otu.ids==adj.nodes[2])
    beta=betaMat[xindex, yindex]
    if (beta>0){
        edge.colors=append(edge.colors, "forestgreen")
    }else if(beta<0){
              edge.colors=append(edge.colors, "red")
    }
}

E(para.ig.mb)$color=edge.colors

para.ig.mb.b=para.ig.mb
nodenames=V(para.ig.mb.b)$name

### generate colors based on family
# we got 15 families
nb.col=nrow(table(pPSpara@tax_table[,c(5)]))+1
coul <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.col)
V(para.ig.mb.b)$family=pPSpara@tax_table[,5]
mycolor <- coul[as.numeric(as.factor(V(para.ig.mb.b)$family))]
E(para.ig.mb.b)$arrow.size=1

#V(para.ig.mb.b)$frame.color="black"
library("network")
library("visNetwork")
library("sna")
library("ndtv")

jpeg(filename = "fig/Network_para.jpg",
    width =15, height = 15, units = "in", res= 300)
plot(para.ig.mb.b, vertex.label=NA, vertex.size=5, vertex.color=mycolor, margin=c(2,2,0,0))
legend(x=-1.6, y=-0.2, legend=levels(as.factor(V(para.ig.mb.b)$family)), col=coul, bty="n", pch=20, pt.cex=2)
dev.off()

######################################## network for overall
#t1 <- Sys.time()
#se.mb.PS <- spiec.easi(pPS, method="mb", pulsar.params=pargs)
#t2 <- Sys.time()
#t2-t1
#saveRDS(se.mb.PS, file="tmp/se.mb.PS.rds")

se.mb.PS <- readRDS(file="tmp/se.mb.PS.rds")
betaMat=symBeta(getOptBeta(se.mb.PS), mode="maxabs")
betaMat

#weights <- Matrix::summary(t(betaMat))[,3]
#weights that are not negative
weights <- (1-Matrix::summary(t(betaMat))[,3])/2
ps.ig.mb <- adj2igraph(getRefit(se.mb.PS),
#                       diag=TRUE,
                       edge.attr=list(weight=weights),
                       vertex.attr=list(name=taxa_names(pPS)))
# plot network eith edges coloured by negative or positive correlation
otu.ids=taxa_names(pPS)
edges <- E(ps.ig.mb)
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
### generate colors based on phyla
get_taxa_unique(pPS, "phylum")
# we got 21 phyla
nb.col=nrow(table(pPS@tax_table[,c(2)]))+1
coul <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.col)
V(ps.ig.mb.b)$phyla=pPS@tax_table[,2]
mycolor <- coul[as.numeric(as.factor(V(ps.ig.mb.b)$phyla))]
E(ps.ig.mb.b)$arrow.size=1

#V(spiec.graph.b)$name=getTaxonomy(nodenames, taxa.f, useRownames=TRUE)
library("igraph")
#E(ps.ig.mb.b)$arrow.size=weight
V(ps.ig.mb.b)$color="white"
V(ps.ig.mb.b)$frame.color="black"

deg <- igraph::degree(ps.ig.mb.b, mode="all")

png(filename = "fig/Network_all_test.png",
        width =15, height = 15, units = "in", res= 300)
plot(ps.ig.mb.b, vertex.label=NA, vertex.size=deg/2, vertex.color=mycolor, margin=c(2,2,0,0), edge.width=E(ps.ig.mb.b)$weight)
legend(x=-1.6, y=-0.2, legend=levels(as.factor(V(ps.ig.mb.b)$phyla)), col=coul, bty="n", pch=20, pt.cex=2)
dev.off()

centr_degree(ps.ig.mb.b, mode="in", normalized=T)
#igraph::closeness(ps.ig.mb.b, mode="all")
#igraph::centr_clo(ps.ig.mb.b, mode="all", normalized=T) 

######## other stuff
edge_density(ps.ig.mb.b)
transitivity(ps.ig.mb.b, type="global")
diameter(ps.ig.mb.b, directed=F)
#cluster analysis
betaMat=as.matrix(symBeta(getOptBeta(se.mb.PS)))
clusters=cluster_fast_greedy(ps.ig.mb.b)
clusters$modularity
clusters$membership
sizes(clusters)
d <- as.dendrogram(clusters)

length(clusters)

plot(d, xlab=NULL)

png(filename = "fig/Network_all_clusters.png",
        width =15, height = 15, units = "in", res= 300)
plot(clusters,ps.ig.mb.b, vertex.label=NA, vertex.size=deg/2, vertex.color=mycolor, margin=c(2,2,0,0), edge.width=E(ps.ig.mb.b)$weight)
legend(x=-1.6, y=-0.2, legend=levels(as.factor(V(ps.ig.mb.b)$phyla)), col=coul, bty="n", pch=20, pt.cex=2)
dev.off()

### which taxa are in the clusters?
cluster1Indices=which(clusters$membership==1)
clusterOneOtus=clusters$names[clusterOneIndices]
cluster2Indices=which(clusters$membership==2)
clusterTwoOtus=clusters$names[clusterTwoIndices]
cluster3Indices=which(clusters$membership==3)
cluster4Indices=which(clusters$membership==4)
cluster5Indices=which(clusters$membership==5)
cluster6Indices=which(clusters$membership==6)
cluster7Indices=which(clusters$membership==7)
cluster8Indices=which(clusters$membership==8)
cluster9Indices=which(clusters$membership==9)
cluster10Indices=which(clusters$membership==10)
cluster11Indices=which(clusters$membership==11)
cluster12Indices=which(clusters$membership==12)
cluster13Indices=which(clusters$membership==13)
cluster14Indices=which(clusters$membership==14)
cluster15Indices=which(clusters$membership==15)
cluster16Indices=which(clusters$membership==16)

sort(table(tax_table(pPS)[cluster1Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster2Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster3Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster4Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster5Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster6Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster7Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster8Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster9Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster10Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster11Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster12Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster13Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster14Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster15Indices,5]),decreasing = TRUE)
sort(table(tax_table(pPS)[cluster16Indices,5]),decreasing = TRUE)

st1.mb <- degree.distribution(ps.ig.mb)
plot(0:(length(st1.mb)-1), st1.mb)
E(ps.ig.mb.b)$weight

nrow(pPS@tax_table[,2])

ceb <- cluster_edge_betweenness(ps.ig.mb.b) # based on edge betweenness (Newman-Girvan)
length(ceb)
dendPlot(ceb, mode="hclust")

# We divide by two since an edge is represented by two entries in the matrix.
positive=length(betaMat[betaMat>0])/2
negative=length(betaMat[betaMat<0])/2
total=length(betaMat[betaMat!=0])/2
negative
positive
total

## Cluster positive edges only
#otu.ids=taxa_names(pPS)
#edges <- E(ps.ig.mb)
filtered.edges=c()
for (e.index in 1:length(edges)){
    adj.nodes=ends(ps.ig.mb, edges[e.index])
    xindex=which(otu.ids==adj.nodes[1])
    yindex=which(otu.ids==adj.nodes[2])
    beta=betaMat[xindex, yindex]
    if (beta>0){
        filtered.edges=c(filtered.edges, edges[e.index])
    }
}

ps.ig.mb.neg=delete_edges(ps.ig.mb.b, filtered.edges)
ps.ig.mb.neg

nfiltered.edges=c()
for (e.index in 1:length(edges)){
    adj.nodes=ends(ps.ig.mb, edges[e.index])
    xindex=which(otu.ids==adj.nodes[1])
    yindex=which(otu.ids==adj.nodes[2])
    beta=betaMat[xindex, yindex]
    if (beta<0){
        nfiltered.edges=c(nfiltered.edges, edges[e.index])
    }
}

#ps.ig.mb.pos=delete_edges(ps.ig.mb.b, nfiltered.edges)

E(ps.ig.mb.neg)$width <- E(ps.ig.mb.neg)$weight/6

plot(ps.ig.mb.pos, vertex.label=NA, vertex.size=5, vertex.color=mycolor, margin=c(2,2,0,0), edge.width=E(ps.ig.mb.pos)$weight,layout=layout_in_circle(ps.ig.mb.pos))
legend(x=-1.6, y=-0.2, legend=levels(as.factor(V(ps.ig.mb.b)$phyla)), col=coul, bty="n", pch=20, pt.cex=2)

plot(ps.ig.mb.neg, vertex.label=NA, vertex.size=5, vertex.color=mycolor, margin=c(2,2,0,0), edge.width=E(ps.ig.mb.neg)$weight, layout=layout_in_circle(ps.ig.mb.neg))
legend(x=-1.6, y=-0.2, legend=levels(as.factor(V(ps.ig.mb.neg)$phyla)), col=coul, bty="n", pch=20, pt.cex=2)


########## genus level

# network
#pargs <- list(rep.num=1000, seed=10010, ncores=20)
#se.mb.PSgen <- spiec.easi(pPSgen, method="mb", pulsar.params=pargs)
#saveRDS(se.mb.PSgen, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/se.mb.PSgen")
se.mb.PSgen <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/se.mb.PSgen")

gen.ig.mb <- adj2igraph(getRefit(se.mb.PSgen), vertex.attr=list(name=taxa_names(pPSgen)))

betaMat=as.matrix(symBeta(getOptBeta(se.mb.PSgen)))

otu.ids=taxa_names(pPSgen)
edges <- E(gen.ig.mb)
edge.colors=c()
for (e.index in 1:length(edges)){
    adj.nodes=ends(gen.ig.mb, edges[e.index])
    xindex=which(otu.ids==adj.nodes[1])
    yindex=which(otu.ids==adj.nodes[2])
    beta=betaMat[xindex, yindex]
    if (beta>0){
        edge.colors=append(edge.colors, "forestgreen")
    }else if(beta<0){
              edge.colors=append(edge.colors, "red")
    }
}

E(gen.ig.mb)$color=edge.colors
gen.ig.mb.b=gen.ig.mb
nodenames=V(gen.ig.mb.b)$name

### generate colors based on family
# we got 20 phyla
nb.col=nrow(table(pPSgen@tax_table[,c(2)]))+1

coul <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.col)

V(gen.ig.mb.b)$family=pPSgen@tax_table[,2]

mycolor <- coul[as.numeric(as.factor(V(gen.ig.mb.b)$family))]

E(gen.ig.mb.b)$arrow.size=1

#V(para.ig.mb.b)$frame.color="black"

jpeg(filename = "fig/Network_gen.jpg",
    width =15, height = 15, units = "in", res= 300)
plot(gen.ig.mb.b, vertex.label=NA, vertex.size=5, vertex.color=mycolor, margin=c(0,0,0,0))
legend(x=-1.6, y=-0.2, legend=levels(as.factor(V(gen.ig.mb.b)$family)), col=coul, bty="n", pch=20, pt.cex=2)
dev.off()

############################### network without parasite taxa
pargs <- list(rep.num=1000, seed=10010, ncores=20)

t1 <- Sys.time()
se.pr.PS <- spiec.easi(pPS, method="mb", pulsar.params=pargs)
t2 <- Sys.time()
t2-t1
saveRDS(se.pr.PS, file="tmp/se.pr.PS.rds")
