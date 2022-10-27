## Please uncomment the first time you run this and re-install packages

#require(devtools)
#devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)

library(ggplot2)
library(MultiAmplicon)
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)
library(parallel)

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- FALSE
doMultiAmp <- FALSE    
doTax <- TRUE

###################Full run Microbiome#######################
#Preparation of files
##These are the same steps that are followed by the DADA2 pipeline
#path <- "/SAN//Susanas_den/HMHZ/data/2018_22_HMHZ_1_2/"
path <- "/SAN/Susanas_den/HMHZ/data/2018_22_HMHZ_2_1/"
#path <- "/SAN//Susanas_den/HMHZ/data/2018_22_HMHZ_2_2/"


fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads

samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("S\\d+-", "\\1", basename(samples)) ##For Pool 2
samples<- gsub("-", "_", basename(samples))

#Quality plots of the reads
pdf("/SAN/Susanas_den/HMHZ/results/2020Aug/qualityProfileF1_2_1.pdf",
    height = 7, width = 7)
plotQualityProfile(fastqF[[1]])
dev.off()

pdf("/SAN/Susanas_den/HMHZ/results/2020Aug/qualityProfileR1_2_1.pdf",
    height = 7, width = 7)
plotQualityProfile(fastqR[[1]])
dev.off()

#Creation of a folder for filtrated reads
filt_path <- "/SAN/Susanas_den/HMHZ/results/2020Aug/filtered2_1"

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
                      truncLen=c(200,150),
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
    filedir <- "/SAN/Susanas_den/HMHZ/results/2020Aug/stratified_files_2_1"
    if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
    MA <- sortAmplicons(MA, n=1e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
    errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                         verbose=0, multithread = 12)
    errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                        verbose=0, multithread = 12)
    MA <- derepMulti(MA, mc.cores=12)
    MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                    verbose=0, mc.cores=12)
    MA <- mergeMulti(MA, mc.cores=12)
    propMerged <- MultiAmplicon::calcPropMerged(MA)
    MA <- mergeMulti(MA, mc.cores=12)
    MA <- makeSequenceTableMulti(MA, mc.cores=12)
    MA <- removeChimeraMulti(MA, mc.cores=12)
    saveRDS(MA, "/SAN/Susanas_den/HMHZ/results/2020Aug/MA2_1.RDS")
} else{
    MA <- readRDS("/SAN/Susanas_den/HMHZ/results/2020Aug/MA2_1.RDS")
}

trackingF=getPipelineSummary(MA)
PipSum=plotPipelineSummary(trackingF)+scale_y_log10() 
ggsave("Sequencing_summary_HMHZ_2_1.pdf", PipSum, path = "/SAN/Susanas_den/HMHZ/results/2020Aug/", height = 15, width = 15)


# save fasta file with all sequences for taxonomic analyses
#MA1 <- getSequenceTableNoChime(MA)
#all.dada.seq <- DNAStringSet(unlist(lapply(MA1, colnames)))
#head(all.dada.seq)
#writeFasta(all.dada.seq, "/SAN/Susanas_den/HMHZ/results/2020May/HMHZ1_1.fasta")

err_F <- plotErrors(errF, nominalQ=TRUE)
pdf("/SAN/Susanas_den/HMHZ/results/2020Aug/Estimeted_error_ratesF_2_1.pdf",
    height = 7, width = 7)
err_F
dev.off()
err_R <- plotErrors(errR, nominalQ=TRUE)
pdf("/SAN/Susanas_den/HMHZ/results/2020Aug/Estimeted_error_ratesR_2_1.pdf",
    height = 7, width = 7)
err_R
dev.off()

trackingF <- getPipelineSummary(MA)
PipSum <- plotPipelineSummary(trackingF)
PipSumlog <- plotPipelineSummary(trackingF) + scale_y_log10()
pdf("/SAN/Susanas_den/HMHZ/results/2020May/plots/Sequencing_summary_HMHZ_2_1.pdf",
    height = 6, width = 8)
PipSum
dev.off()

Heatmap <- plotAmpliconNumbers(MA)
pdf("/SAN/Susanas_den/HMHZ/results/2020Aug/heat_Sequencing_summary_HMHZ_2_1.pdf",
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
                    infasta = "/SAN/Susanas_den/HMHZ/tmp/HMHZ2_1.fasta",
                    outblast = "/SAN/Susanas_den/HMHZ/results/2020Aug/blast2_1_out.fasta",
                    taxonSQL = "/SAN/db/taxonomy/update20200827/accessionTaxa.sql",
                    num_threads = 20)
saveRDS(MA, file="/SAN/Susanas_den/HMHZ/results/2020Aug/MA2_1Tax.Rds") ##Just Test run HMHZ 1

##Start from here after the taxonomic annotation
MA<- readRDS(file= "/SAN/Susanas_den/HMHZ/results/2020Aug/MA2_1Tax.Rds") ###Test run

###Load sample information

sample.data <- read.csv("/SAN/Susanas_den/HMHZ/data/Sample_selection_Metabarcoding_Complete.csv", dec=",", stringsAsFactors=FALSE)

summary(sample.data)

##Tiny adjustments to metadata
require(dplyr)

sample.data$X.1<- NULL
sample.data$X<- NULL
sample.data$Run_I_ID <- NULL
sample.data$Run_II_ID <- NULL
sample.data$Year <- as.factor(sample.data$Year)
sample.data$HI <- as.numeric(sample.data$HI)
sample.data$HI<-round(sample.data$HI, 2)
sample.data$OPG_Eimeria <- as.numeric(sample.data$OPG_Eimeria)
sample.data <- sample.data[!is.na(sample.data$HI),] ##Eliminate 3 samples without HI

sample.data$Concentration <- as.numeric(sample.data$Concentration)
sample.data$Chip_number <- as.factor(sample.data$Chip_number)
sample.data$Body_weight <- as.numeric(sample.data$Body_weight)
sample.data$Body_length <- as.numeric(sample.data$Body_length)
sample.data$Sex <- as.factor(sample.data$Sex)

sample.data%>%
    mutate(BMI = Body_weight/((Body_length)^2)) -> sample.data

sample.data%>%
    mutate(Seq_Run = case_when(Chip_number%in% c(1:8) ~ "Pool_1",
                                     Chip_number%in% c(9:15) ~ "Pool_2")) -> sample.data

sample.data$Seq_Run <- as.factor(sample.data$Seq_Run)

sample.data$Longitude<- as.numeric(sample.data$Longitude)
sample.data$Latitude<- as.numeric(sample.data$Latitude)
sample.data$Longitude<-round(sample.data$Longitude, 4)
sample.data$Latitude<-round(sample.data$Latitude, 4)
sample.data$Locality <- paste(sample.data$Latitude, sample.data$Longitude)

sample.data$OPG_Eimeria <- as.numeric(sample.data$OPG_Eimeria)
sample.data$delta_ct_MminusE <- as.numeric(sample.data$delta_ct_MminusE)

sample.data$Transect <- gsub(" ", "", sample.data$Transect)
sample.data$Transect <- as.factor(sample.data$Transect)

sample.data %>%
    mutate(Genotype = case_when(HI >= 0.95 ~ "Mmm",
                                   HI <= 0.05 ~ "Mmd", HI < 0.95 | HI > 0.05 ~ "Hybrid")) -> sample.data

sample.data$Genotype <- as.factor(sample.data$Genotype)

sample.data %>%
    mutate(Pinworms = sample.data$Aspiculuris_tetraptera+sample.data$Syphacia_obvelata) -> sample.data ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

sample.data %>%
    mutate(Eimeria = case_when(Species_tissue== "E_ferrisi"| Species_tissue== "E_falciformis"| Species_tissue== "E_vermiformis"| Species_tissue== "Other" ~ "Positive",
                                  Species_tissue== "Negative" ~ "Negative",
                                  Flot== "FALSE" & Ap5== "FALSE" ~ "Negative")) -> sample.data

sample.data$Eimeria <- as.factor(sample.data$Eimeria)

rownames(sample.data) <- make.unique(sample.data$Mouse_ID) ##Works when MA contains single run data

##Adding sample data

MAsample <- addSampleData(MA, sample.data)

#plotAmpliconNumbers(MAsample[, which(colnames(MAsample)%in%
#sample.data$Mouse_ID)])


##To phyloseq
##Sample data
PS <- toPhyloseq(MAsample, colnames(MAsample)) ##Now it work
saveRDS(PS, file="/SAN/Susanas_den/HMHZ/results/2020Aug/PhyloSeqCombi_HMHZ_2_1.Rds") ###Results from preliminary analysis (Sample data)

sum(otu_table(PS)) ##Total denoised reads

##Primer data
PS.l <- toPhyloseq(MAsample, colnames(MAsample),  multi2Single=FALSE) ##It work
###For primer analysis (Victor)
##saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqList_HMHZ_1_1.Rds") ###Test run pool 1
saveRDS(PS.l, file="/SAN/Susanas_den/HMHZ/results/2020Aug/PhyloSeqList_HMHZ_2_1.Rds") ###Full run Pool 1

###
lapply(getTaxonTable(MAsample), function (x) table(as.vector(x[, "phylum"])))
lapply(getTaxonTable(MA), function (x) table(as.vector(x[, "genus"])))
lapply(getTaxonTable(MAsample), function (x) table(as.vector(x[, "species"])))
