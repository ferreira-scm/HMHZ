# MultiAmplicon (dada2) preprocessing of wild samples
setwd('/SAN/Susanas_den/gitProj/HMHZ/')

library(rlang)
library(ggplot2)
library(MultiAmplicon, lib.loc="/usr/local/lib/R/site-library/")
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)
library(parallel)
library(dada2)
library(phyloseq)
library(utils)
library(Biostrings)
library("tidyverse")
library(dplyr)



## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- FALSE
doMultiAmp <- FALSE
doTax <- FALSE

###################Full run Microbiome#######################
#Preparation of files
##These are the same steps that are followed by the DADA2 pipeline
path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_1_1/"
#path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_1_2/"
#path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_2_1/"
#path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_2_2/"
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads
samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("s\\d+-", "\\1", basename(samples)) ##For Pool 1
#samples<- gsub("S\\d+-", "\\1", basename(samples)) ##For Pool 2
samples<- gsub("-", "_", basename(samples))
#plotQualityProfile(fastqF[[20]])
#plotQualityProfile(fastqR[[200]])
#Creation of a folder for filtrated reads
filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/filtered1_1"
#filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/filtered1_2"
#filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/filtered2_1"
#filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/filtered2_2"
#Pipeline filtration of pair-end reads
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples
## some files will be filtered out completely
# run 1_1 truncLen=c(200,180)
# run 1_2 truncLen=c(200,200)
# run 2_1 truncLen=c(200,200)
# run 2_2 truncLen=c(200,180)
if(doFilter){
    lapply(seq_along(fastqF),  function (i) {
        filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                      truncLen=c(200,200),
                      maxN=0, maxEE=2, truncQ=2,
                      compress=TRUE, verbose=TRUE)
    })
}
names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)
#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel
#Primers used in the arrays
ptable <- read.csv(file = "/SAN/Susanas_den/gitProj/HMHZ/data/primer_list.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "Seq_F"]
primerR <- ptable[, "Seq_R"]
names(primerF) <- as.character(ptable[, "Name_F"])
names(primerR) <- as.character(ptable[, "Name_R"])
primer <- PrimerPairsSet(primerF, primerR)
##Multi amplicon pipeline
#We start by sorting our amplicons by primer sequences cutting off the latter from sequencing reads. The directory for sorted amplicons must be empty before that.
if(doMultiAmp){
    MA <- MultiAmplicon(primer, files)
    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/stratified_files_1_1"
#    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/stratified_files_1_2"
#    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/stratified_files_2_1"
#    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/stratified_files_2_2"
    if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
    MA <- sortAmplicons(MA, n=2e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
    errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                         verbose=0, multithread = 20)
    errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                        verbose=0, multithread = 20)
    MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                    verbose=0, mc.cores=20)
    MA <- mergeMulti(MA, mc.cores=20)
    propMerged <- MultiAmplicon::calcPropMerged(MA)
    MA <- mergeMulti(MA, mc.cores=20, justConcatenate=propMerged<0.8)
    MA <- makeSequenceTableMulti(MA, mc.cores=20)
    MA <- removeChimeraMulti(MA, mc.cores=20)
    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1.RDS")
#    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2.RDS")
#    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1.RDS")
#    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2.RDS")
} else{
    MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1.RDS")
#   MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2.RDS")
#   MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1.RDS")
#   MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2.RDS")
}
###New taxonomic assignment

if(doTax){MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/1_1.fasta",
                    outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/blast1_1_out.fasta",
                    #infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/1_2.fasta",
                    #outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/blast1_2_out.fasta",
                    #infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/2_1.fasta",
                    #outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/blast2_1_out.fasta",
                    #infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/2_2.fasta",
                    #outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/blast2_2_out.fasta",
                    taxonSQL = "/SAN/db/taxonomy/accessionTaxa.sql",
                    num_threads = 20)
saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1Tax.Rds")
#saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2Tax.Rds")
#saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1Tax.Rds")
#saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2Tax.Rds")
##Start from here after the taxonomic annotation
} else{
    MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1Tax.Rds")
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2Tax.Rds")
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1Tax.Rds")
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2Tax.Rds")
}

colnames(MA)

if(!exists("sample.data")){
    source("/SAN/Susanas_den/gitProj/HMHZ/bin/1_data_prep.R")
}

##Adding sample data
#MAsample <- addSampleData(MA, sample.data)

##To phyloseq
##Sample data
PS <- toPhyloseq(MA, colnames(MA)) ##Now it work
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi11.Rds")

PS@sam_data <- sample_data(sample.data)

head(PS@sam_data)


saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi11.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi12.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi21.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi22.Rds")
sum(otu_table(PS)) ##Total denoised reads
### this is giving errors
##Primer data

PS.l <- toPhyloseq(MAsample, colnames(MAsample),  multi2Single=FALSE)
###For primer analysis (Victor)
saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist11.Rds")
#saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist12.Rds")
#saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist21.Rds")
#saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist22.Rds")
#### to erase or organize

###################Full run Microbiome 1_2#######################
#Preparation of files
##These are the same steps that are followed by the DADA2 pipeline
#path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_1_1/"
path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_1_2/"
#path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_2_1/"
#path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_2_2/"
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads
samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("s\\d+-", "\\1", basename(samples)) ##For Pool 1
#samples<- gsub("S\\d+-", "\\1", basename(samples)) ##For Pool 2
samples<- gsub("-", "_", basename(samples))
#plotQualityProfile(fastqF[[200]])
#plotQualityProfile(fastqR[[30]])
#Creation of a folder for filtrated reads
#filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/filtered1_1"
filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/filtered1_2"
#filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/filtered2_1"
#filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/filtered2_2"
#Pipeline filtration of pair-end reads
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples
## some files will be filtered out completely
# run 1_1 truncLen=c(200,180)
# run 1_2 truncLen=c(200,200)
# run 2_1 truncLen=c()
# run 2_2 truncLen=c()
if(doFilter){
    lapply(seq_along(fastqF),  function (i) {
        filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                      truncLen=c(200,200),
                      maxN=0, maxEE=2, truncQ=2,
                      compress=TRUE, verbose=TRUE)
    })
}
names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)
#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel
#Primers used in the arrays
ptable <- read.csv(file = "/SAN/Susanas_den/gitProj/HMHZ/data/primer_list.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "Seq_F"]
primerR <- ptable[, "Seq_R"]
names(primerF) <- as.character(ptable[, "Name_F"])
names(primerR) <- as.character(ptable[, "Name_R"])
primer <- PrimerPairsSet(primerF, primerR)
##Multi amplicon pipeline
#We start by sorting our amplicons by primer sequences cutting off the latter from sequencing reads. The directory for sorted amplicons must be empty before that.
if(doMultiAmp){
    MA <- MultiAmplicon(primer, files)
#    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/stratified_files_1_1"
    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/stratified_files_1_2"
#    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/stratified_files_2_1"
#    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/stratified_files_2_2"
    if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
    MA <- sortAmplicons(MA, n=2e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
    errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                         verbose=0, multithread = 20)
    errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                        verbose=0, multithread = 20)
    MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                    verbose=0, mc.cores=20)
    MA <- mergeMulti(MA, mc.cores=20)
    propMerged <- MultiAmplicon::calcPropMerged(MA)
    MA <- mergeMulti(MA, mc.cores=20, justConcatenate=propMerged<0.8)
    MA <- makeSequenceTableMulti(MA, mc.cores=20)
    MA <- removeChimeraMulti(MA, mc.cores=20)
#    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1.RDS")
    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2.RDS")
#    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1.RDS")
#    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2.RDS")
} else{
#    MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1.RDS")
   MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2.RDS")
#   MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1.RDS")
#   MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2.RDS")
}
###New taxonomic assignment
MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    #infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/1_1.fasta",
                    #outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/blast1_1_out.fasta",
                    infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/1_2.fasta",
                    outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/blast1_2_out.fasta",
                    #infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/2_1.fasta",
                    #outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/blast2_1_out.fasta",
                    #infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/2_2.fasta",
                    #outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/blast2_2_out.fasta",
                    taxonSQL = "/SAN/db/taxonomy/accessionTaxa.sql",
                    num_threads = 20)
#saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1Tax.Rds")
saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2Tax.Rds")
#saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1Tax.Rds")
#saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2Tax.Rds")
##Start from here after the taxonomic annotation
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1Tax.Rds")
MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2Tax.Rds")
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1Tax.Rds")
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2Tax.Rds")
if(!exists("sample.data")){
    source("/SAN/Susanas_den/gitProj/HMHZ/bin/1_data_prep.R")
}
##Adding sample data
#MAsample <- addSampleData(MA, sample.data)
##To phyloseq
##Sample data
PS <- toPhyloseq(MA, colnames(MA)) ##Now it work
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi11.Rds")
PS@sam_data <- sample_data(sample.data)
saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi12.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi21.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi22.Rds")
sum(otu_table(PS)) ##Total denoised reads
### this is giving errors
##Primer data
PS.l <- toPhyloseq(MA, colnames(MA),  multi2Single=FALSE)
###For primer analysis (Victor)
#saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist11.Rds")
saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist12.Rds")
#saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist21.Rds")
#saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist22.Rds")
### and again
###################Full run Microbiome 2_1#######################
#Preparation of files
##These are the same steps that are followed by the DADA2 pipeline
#path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_1_1/"
#path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_1_2/"
path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_2_1/"
#path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_2_2/"
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads
samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
#samples<- gsub("s\\d+-", "\\1", basename(samples)) ##For Pool 1
samples<- gsub("S\\d+-", "\\1", basename(samples)) ##For Pool 2
samples<- gsub("-", "_", basename(samples))
#plotQualityProfile(fastqF[[10]])
#plotQualityProfile(fastqR[[30]])
#Creation of a folder for filtrated reads
#filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/filtered1_1"
#filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/filtered1_2"
filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/filtered2_1"
#filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/filtered2_2"
#Pipeline filtration of pair-end reads
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples
## some files will be filtered out completely
# run 1_1 truncLen=c(200,180)
# run 1_2 truncLen=c(200,200)
# run 2_1 truncLen=c(200,200)
# run 2_2 truncLen=c(200,180)
if(doFilter){
    lapply(seq_along(fastqF),  function (i) {
        filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                      truncLen=c(200,200),
                      maxN=0, maxEE=2, truncQ=2,
                      compress=TRUE, verbose=TRUE)
    })
}
names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)
#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel
#Primers used in the arrays
ptable <- read.csv(file = "/SAN/Susanas_den/gitProj/HMHZ/data/primer_list.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "Seq_F"]
primerR <- ptable[, "Seq_R"]
names(primerF) <- as.character(ptable[, "Name_F"])
names(primerR) <- as.character(ptable[, "Name_R"])
primer <- PrimerPairsSet(primerF, primerR)
##Multi amplicon pipeline
#We start by sorting our amplicons by primer sequences cutting off the latter from sequencing reads. The directory for sorted amplicons must be empty before that.
if(doMultiAmp){
    MA <- MultiAmplicon(primer, files)
#    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/stratified_files_1_1"
#    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/stratified_files_1_2"
    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/stratified_files_2_1"
#    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/stratified_files_2_2"
    if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
    MA <- sortAmplicons(MA, n=2e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
    errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                         verbose=0, multithread = 20)
    errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                        verbose=0, multithread = 20)
    MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                    verbose=0, mc.cores=20)
    MA <- mergeMulti(MA, mc.cores=20)
    propMerged <- MultiAmplicon::calcPropMerged(MA)
    MA <- mergeMulti(MA, mc.cores=20, justConcatenate=propMerged<0.8)
    MA <- makeSequenceTableMulti(MA, mc.cores=20)
    MA <- removeChimeraMulti(MA, mc.cores=20)
#    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1.RDS")
#    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2.RDS")
    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1.RDS")
#    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2.RDS")
} else{
#    MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1.RDS")
#   MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2.RDS")
   MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1.RDS")
#   MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2.RDS")
}
###New taxonomic assignment
MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    #infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/1_1.fasta",
                    #outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/blast1_1_out.fasta",
                    #infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/1_2.fasta",
                    #outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/blast1_2_out.fasta",
                    infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/2_1.fasta",
                    outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/blast2_1_out.fasta",
                    #infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/2_2.fasta",
                    #outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/blast2_2_out.fasta",
                    taxonSQL = "/SAN/db/taxonomy/accessionTaxa.sql",
                    num_threads = 20)
#saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1Tax.Rds")
#saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2Tax.Rds")
saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1Tax.Rds")
#saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2Tax.Rds")
##Start from here after the taxonomic annotation
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1Tax.Rds")
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2Tax.Rds")
MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1Tax.Rds")
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2Tax.Rds")
if(!exists("sample.data")){
    source("/SAN/Susanas_den/gitProj/HMHZ/bin/1_data_prep.R")
}
##Adding sample data
#MAsample <- addSampleData(MA, sample.data)
##To phyloseq
##Sample data
PS <- toPhyloseq(MA, colnames(MA)) ##Now it work
PS@sam_data <- sample_data(sample.data)
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi11.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi12.Rds")
saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi21.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi22.Rds")
sum(otu_table(PS)) ##Total denoised reads
##Primer data
PS.l <- toPhyloseq(MA, colnames(MA),  multi2Single=FALSE)
###For primer analysis (Victor)
#saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist11.Rds")
#saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist12.Rds")
saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist21.Rds")
#saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist22.Rds")
### and again###################Full run Microbiome 2_2#######################
#Preparation of files
##These are the same steps that are followed by the DADA2 pipeline
#path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_1_1/"
#path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_1_2/"
#path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_2_1/"
path <- "/SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_2_2/"
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads
samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
#samples<- gsub("s\\d+-", "\\1", basename(samples)) ##For Pool 1
samples<- gsub("S\\d+-", "\\1", basename(samples)) ##For Pool 2
samples<- gsub("-", "_", basename(samples))
#plotQualityProfile(fastqF[[10]])
#plotQualityProfile(fastqR[[10]])
#Creation of a folder for filtrated reads
#filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/filtered1_1"
#filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/filtered1_2"
#filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/filtered2_1"
filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/filtered2_2"
#Pipeline filtration of pair-end reads
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples
## some files will be filtered out completely
# run 1_1 truncLen=c(200,180)
# run 1_2 truncLen=c(200,200)
# run 2_1 truncLen=c(200,200)
# run 2_2 truncLen=c(200,180)
if(doFilter){
    lapply(seq_along(fastqF),  function (i) {
        filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                      truncLen=c(200,200),
                      maxN=0, maxEE=2, truncQ=2,
                      compress=TRUE, verbose=TRUE)
    })
}
names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)
#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel
#Primers used in the arrays
ptable <- read.csv(file = "/SAN/Susanas_den/gitProj/HMHZ/data/primer_list.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "Seq_F"]
primerR <- ptable[, "Seq_R"]
names(primerF) <- as.character(ptable[, "Name_F"])
names(primerR) <- as.character(ptable[, "Name_R"])
primer <- PrimerPairsSet(primerF, primerR)
##Multi amplicon pipeline
#We start by sorting our amplicons by primer sequences cutting off the latter from sequencing reads. The directory for sorted amplicons must be empty before that.
if(doMultiAmp){
    MA <- MultiAmplicon(primer, files)
#    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/stratified_files_1_1"
#    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/stratified_files_1_2"
#    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/stratified_files_2_1"
    filedir <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/stratified_files_2_2"
    if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
    MA <- sortAmplicons(MA, n=2e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
    errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                         verbose=0, multithread = 20)
    errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                        verbose=0, multithread = 20)
    MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                    verbose=0, mc.cores=20)
    MA <- mergeMulti(MA, mc.cores=20)
    propMerged <- MultiAmplicon::calcPropMerged(MA)
    MA <- mergeMulti(MA, mc.cores=20, justConcatenate=propMerged<0.8)
    MA <- makeSequenceTableMulti(MA, mc.cores=20)
    MA <- removeChimeraMulti(MA, mc.cores=20)
#    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1.RDS")
#    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2.RDS")
#    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1.RDS")
    saveRDS(MA, "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2.RDS")
} else{
#    MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1.RDS")
#   MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2.RDS")
#   MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1.RDS")
   MA <- readRDS("/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2.RDS")
}
###New taxonomic assignment
MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    #infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/1_1.fasta",
                    #outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/blast1_1_out.fasta",
                    #infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/1_2.fasta",
                    #outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/blast1_2_out.fasta",
                    #infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/2_1.fasta",
                    #outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/blast2_1_out.fasta",
                    infasta = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/2_2.fasta",
                    outblast = "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/blast2_2_out.fasta",
                    taxonSQL = "/SAN/db/taxonomy/accessionTaxa.sql",
                    num_threads = 20)
#saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1Tax.Rds")
#saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2Tax.Rds")
#saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1Tax.Rds")
saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2Tax.Rds")
##Start from here after the taxonomic annotation
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_1Tax.Rds")
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/MA1_2Tax.Rds")
#MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_1Tax.Rds")
MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/MA2_2Tax.Rds")
if(!exists("sample.data")){
    source("/SAN/Susanas_den/gitProj/HMHZ/bin/1_data_prep.R")
}
##Adding sample data
#MAsample <- addSampleData(MA, sample.data)
##To phyloseq
##Sample data
PS <- toPhyloseq(MA, colnames(MA)) ##Now it work
PS@sam_data <- sample_data(sample.data)
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi11.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi12.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi21.Rds")
saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi22.Rds")
sum(otu_table(PS)) ##Total denoised reads
### this is giving errors
##Primer data
PS.l <- toPhyloseq(MA, colnames(MA),  multi2Single=FALSE)
###For primer analysis (Victor)
#saveRDS(PS.l11, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist11.Rds")
#saveRDS(PS.l12, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist12.Rds")
#saveRDS(PS.l21, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist21.Rds")
saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist22.Rds")


