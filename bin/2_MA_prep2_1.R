#!/usr/bin/Rscript

#require(devtools)
#devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)
library(ggplot2)
library(dada2)
#library(MultiAmplicon, lib.loc="/usr/local/lib/R/site-library/")
library(reshape)
library(phyloseq)
library(data.table, lib.loc="/usr/local/lib/R/site-library/")
library(taxonomizr)
library(taxize)
library(parallel)

## using the devel
devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- TRUE
doMultiAmp <- TRUE    
doTax <- FALSE

###################Full run Microbiome#######################
#Preparation of files
##These are the same steps that are followed by the DADA2 pipeline
path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_2_1/"
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads
samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("s\\d+-", "\\1", basename(samples)) ##For Pool 1
samples<- gsub("-", "_", basename(samples))

                                        #Quality plots of the reads
#pdf("fig/quality/qualityProfileF1_2_1.pdf", height = 7, width = 7)
plotQualityProfile(fastqF[[1]])
#dev.off()
#pdf("fig/quality/qualityProfileR1_2_1.pdf", height = 7, width = 7)
#plotQualityProfile(fastqR[[1]])
#dev.off()

                                        #Creation of a folder for filtrated reads
filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/interData/filtered2_1"
#Pipeline filtration of pair-end reads
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples
## some files will be filtered out completely
if(doFilter){
    lapply(seq_along(fastqF),  function (i) {
        filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
#                      truncLen=c(260,230),
                      minLen=c(200,200),
                      maxN=0, maxEE=2, truncQ=2,
                      compress=TRUE, verbose=TRUE)
    })
}
names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)
#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel
#Primers used in the arrays
ptable <- read.csv(file = "/SAN/Susanas_den/HMHZ/data/primer_list.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "Seq_F"]
primerR <- ptable[, "Seq_R"]
names(primerF) <- as.character(ptable[, "Name_F"])
names(primerR) <- as.character(ptable[, "Name_R"])
primer <- PrimerPairsSet(primerF, primerR)

##Multi amplicon pipeline
#We start by sorting our amplicons by primer sequences cutting off the latter from sequencing reads. The directory for sorted amplicons must be empty before that.
if(doMultiAmp){
    MA <- MultiAmplicon(primer, files)
    filedir <- "tmp/interData/stratified_files_2_1"
    if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
    MA <- sortAmplicons(MA, n=1e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
    errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                         verbose=0, multithread = 90)
    errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                        verbose=0, multithread = 90)
#    MA <- derepMulti(MA, mc.cores=90) deprecated, no longer needed
    MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                    verbose=0, mc.cores=90)
    MA <- mergeMulti(MA, mc.cores=90)
    propMerged <- MultiAmplicon::calcPropMerged(MA)
    summary(propMerged)
    table(propMerged<0.8)
    MA <- mergeMulti(MA, mc.cores=90) 
    MA <- makeSequenceTableMulti(MA, mc.cores=90)
    MA <- removeChimeraMulti(MA, mc.cores=90)
    saveRDS(MA, "tmp/interData/MA2_1.RDS")
} else{
    MA <- readRDS("tmp/interData/MA2_1.RDS")
}


                                        #trackingF <- getPipelineSummary(MA)
#PipSum <- plotPipelineSummary(trackingF)+scale_y_log10() 
#ggsave("Sequencing_summary_HMHZ_2_1.pdf", PipSum, path = "fig/quality/", height = 15, width = 15)
# save fasta file with all sequences for taxonomic analyses
#MA1 <- getSequenceTableNoChime(MA)
#all.dada.seq <- DNAStringSet(unlist(lapply(MA1, colnames)))
#head(all.dada.seq)
#writeFasta(all.dada.seq, "/SAN/Susanas_den/HMHZ/results/2020May/HMHZ1_1.fasta")
#err_F <- plotErrors(errF, nominalQ=TRUE)
#pdf("fig/quality/Estimeted_error_ratesF_2_1.pdf",
#    height = 7, width = 7)
#err_F
#dev.off()
#err_R <- plotErrors(errR, nominalQ=TRUE)
#pdf("fig/quality/Estimeted_error_ratesR_2_1.pdf",
#    height = 7, width = 7)
#err_R
#dev.off()
#Heatmap <- plotAmpliconNumbers(MA)
#pdf("fig/quality/heat_Sequencing_summary_HMHZ_2_1.pdf",
#    height = 15, width = 15)
#Heatmap
#dev.off()
###New taxonomic assignment
#MA <- blastTaxAnnot(MA,
#                    db = "/SAN/db/blastdb/nt/nt",
#                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
#                    infasta = "tmp/interData/HMHZ2_1.fasta",
#                    outblast = "tmp/interData/blast2_1_out.fasta",
#                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql",
#                    num_threads = 90)

source("bin/Primer_target.R")
p.df <- p.df[match(primer@names, p.df$Primer_name),]
primer@names==p.df$Primer_name

taxT1 <- list()
seqs <- getSequencesFromTable(MA)
seqs <- lapply(seqs, DNAStringSet)
for (i in 1:48){
    if (p.df$Gen[i]=="16S"){
        try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
          "/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/Slv138.dada2.fa",
          multithread=90,
                                    tryRC = TRUE,
                                   verbose=TRUE))
    }
    else if (p.df$Gen[i]=="18S"){
        try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
                 "/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/Slv138.dada2.fa",
                                     multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE))
    }
    else if (p.df$Gen[i]=="28S"){
        try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
                      "/SAN/Susanas_den/AmpMarkers/RESCRIPt/LSURef_NR99/Fastas/Slv138LSU.dada2.fa",
                                     multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE))
    }   
    else if (p.df$Gen[i]=="ITS"){
     try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
                                      "/SAN/Susanas_den/AmpMarkers/UNITE/sh_general_release_s_all_10.05.2021/sh_general_release_dynamic_s_all_10.05.2021.fasta",
                                     multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE))
    }
    else {
     try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
           "/SAN/Susanas_den/AmpMarkers/RESCRIPt/other/Fastas/other.dada2.fa",
                                     multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE))
    }   
}

MA@taxonTable <- taxT1

saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/interData/MA2_1Tax.Rds") #





##Start from here after the taxonomic annotation
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/interData/MA2_1Tax.Rds") ###Test run
##To phyloseq
PS <- TMPtoPhyloseq(MA, colnames(MA)) ##Now it work
############# add metadata
###Load sample information
source("bin/1_LoadingSOTA.R")
### eh we need to adjust rownames
sample_names(PS) <- gsub("S\\d*_", "", rownames(PS@otu_table))   
# sanity check
rownames(PS@otu_table)
rownames(PS@sam_data)
meta <- sota[match(rownames(PS@sam_data), sota$Mouse_ID),] 
#sanity checks
rownames(PS@otu_table)[!rownames(PS@otu_table)==rownames(meta)]
PS@sam_data <- sample_data(meta)
#
#sanity check
#PS@sam_data[which(!rownames(PS@sam_data)==rownames(PS@otu_table))]
#PS@sam_data[which(rownames(PS@sam_data)==rownames(PS@otu_table))]
#
rownames(PS@sam_data) <- rownames(PS@otu_table)
#
# another sanity check
sample_names(PS)==rownames(PS@sam_data)
#
PS_neg <- subset_samples(PS, grepl("NE",rownames(PS@otu_table)))   
#
PS@sam_data$Control <- FALSE
PS@sam_data$Control[which(sample_names(PS)%in%sample_names(PS_neg))] <- TRUE
# sanity check
PS@sam_data$Mouse_ID[PS@sam_data$Control==FALSE]
rownames(PS@sam_data)[PS@sam_data$Control==TRUE]
#
library("decontam")
###### removing contaminants
## assuming that negative controls have 0 DNA
#PS@sam_data$Concentration[PS@sam_data$Control==TRUE] <- 0.0001
#
## ----see-depths---------------------------------------------------------------
#df <- as.data.frame(sample_data(PS)) # Put sample_data into a ggplot-friendly data.frame
#df$LibrarySize <- sample_sums(PS)
#df <- df[order(df$LibrarySize),]
#df$Index <- seq(nrow(df))
#ggplot(data=df, aes(x=Index, y=LibrarySize, color=Control)) + geom_point()
#
ps <- phyloseq::prune_samples(sample_sums(PS)>0, PS)
#
contamdf.freq <- isContaminant(ps, method="either", conc="Concentration", neg="Control", threshold=c(0.1,0.5), normalize=TRUE)
#
table(contamdf.freq$contaminant)
#
### taxa to remove
ps@tax_table[rownames(contamdf.freq[contamdf.freq$contaminant==TRUE,]),5]
#
## quick plotting for contaminant prevalence
#ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
#ps.pa.neg <- prune_samples(sample_data(ps.pa)$Control == TRUE, ps.pa)
#ps.pa.pos <- prune_samples(sample_data(ps.pa)$Control == FALSE, ps.pa)
# Make data.frame of prevalence in positive and negative samples
#df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
#                    contaminant=contamdf.prev$contaminant)
#ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
#      xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
#
## let's remove them now and negative controls
Keep <- rownames(contamdf.freq[contamdf.freq$contaminant==FALSE,])
PS <- prune_samples(sample_data(PS)$Control == FALSE, PS)
PS <- prune_taxa(Keep, PS)
#
saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/interData/PhyloSeqCombi_HMHZ_2_1.Rds") ###Results from preliminary analysis (Sample data)
#
sum(otu_table(PS)) ##Total denoised reads
#
##Primer data
PS.l <- TMPtoPhyloseq(MA, colnames(MA),  multi2Single=FALSE) ##It work
#
sample_names(PS.l[[1]])
#
### eh we need to adjust rownames
for (i in 1:48){
try(sample_names(PS.l[[i]]) <- gsub("S\\d*_", "", rownames(PS.l[[i]]@otu_table))   , silent = TRUE)
}
#
## adding metadata, removing contaminants and controls
neg <- sample_names(subset_samples(PS.l[[1]], !grepl("NE",rownames(PS.l[[1]]@otu_table))))
#
for (i in 1:48) {
    try(PS.l[[i]] <- prune_taxa(Keep, PS.l[[i]]), silent=TRUE)
    try(PS.l[[i]] <- prune_samples(neg, PS.l[[i]]), silent=TRUE)
}
for (i in 1:48) {
    try(PS.l[[i]]@sam_data <- PS@sam_data, silent=TRUE)
}
###For primer analysis (Victor)
saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/interData/PhyloSeqList_HMHZ_2_1.Rds") ###Full run Pool 3

neg

sample_names(PS.l[[1]])


