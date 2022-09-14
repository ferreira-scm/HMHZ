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

library(metagMisc)

PS=readRDS(file="tmp/PSpre.R")
pPS=readRDS(file="tmp/pPSpre.R")

p10PPS=phyloseq_filter_prevalence(pPS, prev.trh=0.1)

pPSnoP10 <- subset_taxa(p10PPS, !family%in%c("Eimeriidae", "Cryptosporidiidae", "Oxyuridae", "Trichuridae", "Heteroxynematidae", "Ascaridiidae", "Spirocercidae", "Heterakidae", "Hymenolepididae", "Sarcocystidae", "Strongyloididae"))

pPSnoP10

se.pr.PS10 = readRDS(file="tmp/se.pr.PS10.rds")
se.PS10 = readRDS(file="tmp/se.PS10.rds")

se.pr.PS10

se.PS10

boot.PS10 <- readRDS(file="tmp/netPS")
randPS <- readRDS(file="tmp/randPS")


net.boot <- vector(mode="list", length=10)
for(i in 1:10) {
    net.boot[[i]] <- boot.PS10[[i]]
}

### adding weight for no parasite network
betaMat=symBeta(getOptBeta(se.pr.PS10), mode="maxabs")
#weights <- Matrix::summary(t(betaMat))[,3]
#weights that are not negative
weights <- (1-Matrix::summary(t(betaMat))[,3])/2
                                        # igraph format
weights

se.pr.PS10

ps.pr.mb <- adj2igraph(getRefit(se.pr.PS10),
                       diag=TRUE,
                       edge.attr=list(weight=weights),
                       vertex.attr=list(name=taxa_names(pPSnoP10)))

### same for complete network
betaMat=symBeta(getOptBeta(se.PS10), mode="maxabs")
weights <- Matrix::summary(t(betaMat))[,3]
#weights that are not negative
weights <- (1-Matrix::summary(t(betaMat))[,3])/2
ps.mb <- adj2igraph(getRefit(se.PS10),
                       diag=TRUE,
                       edge.attr=list(weight=weights),
                    vertex.attr=list(name=taxa_names(p10PPS)))

### and now for bootstrap networks
inet.boot <- vector(mode="list", length=10)

for(i in 1:10) {
    betaMat=symBeta(getOptBeta(net.boot[[i]]), mode="maxabs")
    weights <- (1-Matrix::summary(t(betaMat))[,3])/2
    inet.boot[[i]] <- adj2igraph(getRefit(net.boot[[i]]),
                       diag=TRUE,
                       edge.attr=list(weight=weights),
                       vertex.attr=list(name=taxa_names(randPS[[i]])))
}

# plot complete network eith edges coloured by negative or positive correlation
otu.ids=taxa_names(p10PPS)

edges <- E(ps.mb)
edge.colors=c()

for (e.index in 1:length(edges)){
    adj.nodes=ends(ps.mb, edges[e.index])
    xindex=which(otu.ids==adj.nodes[1])
    yindex=which(otu.ids==adj.nodes[2])
    beta=betaMat[xindex, yindex]
    if (beta>0){
        edge.colors=append(edge.colors, "forestgreen")
    }else if(beta<0){
              edge.colors=append(edge.colors, "red")
    }
}

edge.colors

E(ps.mb)$color=edge.colors
ps.mb.b=ps.mb
nodenames=V(ps.mb.b)$name

### generate colors based on phyla
get_taxa_unique(p10PPS, "phylum")

# we got 12 phyla
nb.col=nrow(table(p10PPS@tax_table[,c(2)]))+1

coul <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.col)
V(ps.mb.b)$phyla=p10PPS@tax_table[,2]
mycolor <- coul[as.numeric(as.factor(V(ps.mb.b)$phyla))]
E(ps.mb.b)$arrow.size=1


#V(spiec.graph.b)$name=getTaxonomy(nodenames, taxa.f, useRownames=TRUE)
library("igraph")
#E(ps.ig.mb.b)$arrow.size=weight
V(ps.mb.b)$color="white"
V(ps.mb.b)$frame.color="black"

deg <- igraph::degree(ps.mb.b, mode="all")

png(filename = "fig/figures4man/Network_all.png",
        width =15, height = 15, units = "in", res= 300)
plot(ps.mb.b, vertex.label=NA, vertex.size=5, vertex.color=mycolor, margin=c(2,2,0,0), edge.width=E(ps.mb.b)$weight)
legend(x=-1.6, y=-0.2, legend=levels(as.factor(V(ps.mb.b)$phyla)), col=coul, bty="n", pch=20, pt.cex=2)
dev.off()

centr_degree(ps.mb.b, mode="in", normalized=T)
#igraph::closeness(ps.ig.mb.b, mode="all")
#igraph::centr_clo(ps.ig.mb.b, mode="all", normalized=T) 

######## other stuff
edge_density(ps.mb.b)
transitivity(ps.mb.b, type="global")
diameter(ps.mb.b, directed=F)

#cluster analysis
betaMat=as.matrix(symBeta(getOptBeta(se.PS10)))
clusters=cluster_fast_greedy(ps.mb.b)

clusters$modularity
clusters$membership
sizes(clusters)
d <- as.dendrogram(clusters)

length(clusters)

plot(d, mode="hclust", labels=FALSE, xlab=NULL)

mycolor

png(filename = "fig/figures4man/Network_all_clusters.png",
        width =15, height = 15, units = "in", res= 300)
plot(clusters,ps.mb.b, vertex.label=NA, vertex.size=5, vertex.color=mycolor, margin=c(2,2,0,0))
legend(x=-1.6, y=-0.2, legend=levels(as.factor(V(ps.mb.b)$phyla)), col=coul, bty="n", pch=20, pt.cex=2)
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

write.csv(sort(table(tax_table(p10PPS)[cluster1Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster1.csv")
write.csv(sort(table(tax_table(p10PPS)[cluster2Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster2.csv")
write.csv(sort(table(tax_table(p10PPS)[cluster3Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster3.csv")
write.csv(sort(table(tax_table(p10PPS)[cluster4Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster4.csv")
write.csv(sort(table(tax_table(p10PPS)[cluster5Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster5.csv")
write.csv(sort(table(tax_table(p10PPS)[cluster6Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster6.csv")
write.csv(sort(table(tax_table(p10PPS)[cluster7Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster7.csv")
write.csv(sort(table(tax_table(p10PPS)[cluster8Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster8.csv")
write.csv(sort(table(tax_table(p10PPS)[cluster9Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster9.csv")
write.csv(sort(table(tax_table(p10PPS)[cluster10Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster10.csv")
write.csv(sort(table(tax_table(p10PPS)[cluster11Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster11.csv")
write.csv(sort(table(tax_table(p10PPS)[cluster12Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster12.csv")
write.csv(sort(table(tax_table(p10PPS)[cluster13Indices,5]),decreasing = TRUE), file="fig/figures4man/cluster13.csv")


st1.mb <- degree.distribution(ps.mb.b)
plot(0:(length(st1.mb)-1), st1.mb)
E(ps.mb.b)$weight

nrow(pPS@tax_table[,2])

ceb <- cluster_edge_betweenness(ps.mb.b, weights=E(ps.mb.b)$weight) # based on edge betweenness (Newman-Girvan)
length(ceb)

ceb$membership

dendPlot(ceb, mode="hclust")

# We divide by two since an edge is represented by two entries in the matrix.
positive=length(betaMat[betaMat>0])/2
negative=length(betaMat[betaMat<0])/2
total=length(betaMat[betaMat!=0])/2
negative
positive
total

