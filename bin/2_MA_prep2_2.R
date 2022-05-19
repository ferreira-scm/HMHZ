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
doTax <- TRUE

###################Full run Microbiome#######################
#Preparation of files
##These are the same steps that are followed by the DADA2 pipeline
path <- "/SAN/Susanas_den/HMHZ/data/2018_22_HMHZ_1_1/"

fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads

samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("s\\d+-", "\\1", basename(samples)) ##For Pool 1
samples<- gsub("-", "_", basename(samples))

#Quality plots of the reads
pdf("fig/quality/qualityProfileF1_2_2.pdf", height = 7, width = 7)
plotQualityProfile(fastqF[[1]])
dev.off()

pdf("fig/quality/qualityProfileR1_2_2.pdf", height = 7, width = 7)
plotQualityProfile(fastqR[[1]])
dev.off()

#Creation of a folder for filtrated reads
filt_path <- "/SAN/Susanas_den/gitProj/HMHZ/tmp/interData/filtered2_2"

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
                      truncLen=c(200,200),
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
    filedir <- "tmp/interData/stratified_files_2_2"
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
    MA <- mergeMulti(MA, justConcatenate=propMerged<0.8, mc.cores=90) 
    MA <- makeSequenceTableMulti(MA, mc.cores=90)
    MA <- removeChimeraMulti(MA, mc.cores=90)
    saveRDS(MA, "tmp/interData/MA2_2.RDS")
} else{
    MA <- readRDS("tmp/interData/MA2_2.RDS")
}

trackingF <- getPipelineSummary(MA)
PipSum <- plotPipelineSummary(trackingF)+scale_y_log10() 
ggsave("Sequencing_summary_HMHZ_2_2.pdf", PipSum, path = "fig/quality/", height = 15, width = 15)


# save fasta file with all sequences for taxonomic analyses
#MA1 <- getSequenceTableNoChime(MA)
#all.dada.seq <- DNAStringSet(unlist(lapply(MA1, colnames)))
#head(all.dada.seq)
#writeFasta(all.dada.seq, "/SAN/Susanas_den/HMHZ/results/2020May/HMHZ1_1.fasta")

err_F <- plotErrors(errF, nominalQ=TRUE)
pdf("fig/quality/Estimeted_error_ratesF_2_2.pdf",
    height = 7, width = 7)
err_F
dev.off()

err_R <- plotErrors(errR, nominalQ=TRUE)
pdf("fig/quality/Estimeted_error_ratesR_2_2.pdf",
    height = 7, width = 7)
err_R
dev.off()

Heatmap <- plotAmpliconNumbers(MA)
pdf("fig/quality/heat_Sequencing_summary_HMHZ_2_2.pdf",
    height = 15, width = 15)
Heatmap
dev.off()

###New taxonomic assignment
#Sys.setenv("BLASTDB" = "/SAN/db/blastdb/") #To make the annotation work, boss will fix this in the package
#library("vctrs", lib.loc="/usr/local/lib/R/site-library")
#MA <- blastTaxAnnot(MA,  dataBaseDir = Sys.getenv("BLASTDB"), negative_gilist = "/SAN/db/blastdb/uncultured.gi", num_threads = 20)

MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    infasta = "tmp/interData/HMHZ2_2.fasta",
                    outblast = "tmp/interData/blast2_2_out.fasta",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql",
                    num_threads = 90)

saveRDS(MA, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/interData/MA2_2Tax.Rds") ##Just Test run HMHZ 1

##Start from here after the taxonomic annotation
MA<- readRDS(file= "/SAN/Susanas_den/gitProj/HMHZ/tmp/interData/MA2_2Tax.Rds") ###Test run

###Load sample information
if(!exists("sample.data")){
    source("bin/1_data_prep.R")
}

##To phyloseq
PS <- toPhyloseq(MA, colnames(MA)) ##Now it work
##Sample data
PS@sam_data <- sample_data(sample.data)

saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/interData/PhyloSeqCombi_HMHZ_2_2.Rds") ###Results from preliminary analysis (Sample data)

sum(otu_table(PS)) ##Total denoised reads

##Primer data
#PS.l <- toPhyloseq(MA, colnames(MAsample),  multi2Single=FALSE) ##It work
###For primer analysis (Victor)
#saveRDS(PS.l, file="/SAN/Susanas_den/HMHZ/results/2020Aug/PhyloSeqList_HMHZ_1_1.Rds") ###Full run Pool 1

###
#lapply(getTaxonTable(MAsample), function (x) table(as.vector(x[, "phylum"])))
