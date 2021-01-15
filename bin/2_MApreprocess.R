# MultiAmplicon (dada2) preprocessing of wild samples

#require(devtools)
#devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)
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

blastTaxAnnot <- function (MA, db="nt/nt",
			       num_threads= getOption("mc.cores", 1L),
			       negative_gilist = system.file("extdata", "uncultured.gi",
							     package = "MultiAmplicon"),
			       infasta=paste0(tempfile(), ".fasta"),
			       outblast=paste0(tempfile(), ".blt"),
			       taxonSQL, ...
			       ) {

## Prevent R CMD check from complaining about the use of pipe expressions
## standard data.table variables
#if (getRversion() >= "2.15.1")
#utils::globalVariables(c(".", ".I", ".N", ".SD"), utils::packageName(), add=F)
#
## Make sure data.table knows we know we're using it
#.datatable.aware = TRUE

SEQ <- getSequencesFromTable(MA)
if(!file.exists(infasta)){
Ssplit <- lapply(seq_along(SEQ), function (i) {
			  if(!is.null(SEQ[[i]])){
			  ampSEQ <- SEQ[[i]]
			  SEQnames <- paste("A", i, "S", 1:length(ampSEQ), "R_", sep="_")
			  names(ampSEQ) <- SEQnames
			  Ssplit <- strsplit(ampSEQ, "NNNNNNNNNN")
			  unlist(Ssplit)
			  } else (NULL)
			  })
Biostrings::writeXStringSet(DNAStringSet(unlist(Ssplit)),
					infasta)
message("wrote file ", infasta,
	" storing input sequences for blast\n")
} else {
message("file ", infasta, " exists, using existing file!",
	" To extract input sequences for blast  again delete the file\n")
}

if(!file.exists(outblast)){
DBmessage <-
paste0("using database ", db, " resulting in as blast database full path\n")
message(DBmessage)

## We blast this file against NR with a gi-list excluding all
## uncultured sequences
## create the gi-list as a download from an NCBI Entrez Nucleotide
## search '"environmental samples"[organism] OR metagenomes[orgn]'
        ## or via command line: esearch -db nucleotide -query '"environmental
        ## samples"[organism] OR metagenomes[orgn]' | efetch -format gi -mode
## text > /SAN/db/blastdb/uncultured_gilist.txt
command <- paste("blastn",
		 "-negative_gilist",  negative_gilist,
		 "-query", infasta,
		 "-db", db,
		 "-evalue",  1e-15,
		 "-num_threads", num_threads,
		 "-max_hsps", 1,
		 "-max_target_seqs 25",
		 "-out", outblast,
		 "-outfmt", "'10 qaccver saccver pident length mismatch gapopen",
		 "qstart qend sstart send evalue bitscore staxid'",
		 ...)
message("STARTED running blast with command:\n",
	command)
system(command)
cat("\nFINISHED running blast\n")
} else {
message("file", outblast, " exists, using existing file!",
	" To run blast again delete file")
}
## we read the blast file
blast <- read.csv(outblast, header=FALSE)
names(blast) <- c("query", "subject", "pident", "length", "mismatch",
		  "gapopen", "qstart", "qend", "sstart", "send", "evalue",
		  "bitscore", "staxid")
blastT <- as.data.table(blast)
blastT$staxid <- as.character(blastT$staxid)
blastT$ampProd <- gsub("_R_\\d+", "", blastT$query)


## retain only the maximal bitscore for each query and staxid
blastT$qtaxid <- as.factor(paste(blastT$query, blastT$staxid, sep="|"))
blastT <- blastT[blastT[, .I[bitscore == max(bitscore)], by=qtaxid][, V1]]

## unique in case of multiple best hits in one taxID
blastT <- blastT[!duplicated(blastT$qtaxid)]

blastT <- blastT[blastT[, .I[bitscore == max(bitscore)], by=qtaxid][, V1]]

## for each amplification product get the sum of the bitscores for
## forward and reverse query (for non-merged sequences that is)
blt <- blastT[, list(bitsum = sum(bitscore)), by=c("ampProd", "staxid")]

## ## now we would hav to generated a taxonomizr sql database
## read.nodes.sql("/SAN/db/taxonomy/nodes.dmp",
		  ##                "/SAN/db/taxonomy/taxonomizr.sql")
## read.names.sql("/SAN/db/taxonomy/names.dmp",
		  ##                "/SAN/db/taxonomy/taxonomizr.sql")
blast.tax <- getTaxonomy(unique(blt$staxid),
			       "/SAN/db/taxonomy/taxonomizr.sql")
blast.tax <- as.data.table(blast.tax, keep.rownames="staxid")
blast.tax$staxid <- gsub("\\s*", "", blast.tax$staxid)
blt <- merge(blt, blast.tax, by="staxid", all=TRUE)

## get the best bitscore for each amplification product
blt <- blt[blt[, .I[bitsum == max(bitsum)], by=ampProd][, V1]]

get.unique.tax <- function(x, Nsupport=1) {
## unique taxa at that level excluding potential NA's
agnostic <- as.character(x)
taxa <- agnostic[!is.na(agnostic)]
ux <- unique(taxa)
## but return NA if they are not unique
if(length(taxa)>=Nsupport && ## number of supporting annotations
	 length(ux)==1){ ## has to be a unique answer
return(ux)
} else {as.character(NA)}
}

## now if multiple taxa with bitscores are given, set NA at any
## level
B <- blt[,
	 {list(superkingdom = get.unique.tax(superkingdom),
			    phylum = get.unique.tax(phylum),
			    class = get.unique.tax(class),
			    order = get.unique.tax(order),
			    family = get.unique.tax(family),
			    genus = get.unique.tax(genus),
			    species = get.unique.tax(species))
	      },
	      by=ampProd]

B$amplicon <- gsub("A_(\\d+)_S_(\\d+)(_R_)?", "\\1", B$ampProd)
B$sequence <- gsub("A_(\\d+)_S_(\\d+)(_R_)?", "\\2", B$ampProd)

annot.l <- by(B, B$amplicon, function (x){
		 taxDF <- as.data.frame(x)
		 n.amp <- as.numeric(unique(x$amplicon))
		 n.seq <- as.numeric(x$sequence)
		 rownames(taxDF) <- SEQ[[n.amp]][n.seq]
		 taxDF <- taxDF[SEQ[[n.amp]],
				   c("superkingdom", "phylum", "class", "order", "family",
				     "genus", "species")]
		 rownames(taxDF) <- SEQ[[n.amp]]
		 as.matrix(taxDF)
		 })
names(annot.l) <- names(SEQ)[as.numeric(names(annot.l))]
annot.l <- lapply(annot.l[names(SEQ)], function (x) x) ## drop array
taxTab.l <- lapply(annot.l, function (x) {
			    if (!is.null(x)) {
			    new("taxonomyTable", x)
			    } else {NULL}
			    })
initialize(MA, taxonTable = taxTab.l)
}



#Sys.setenv("BLASTDB" = "/SAN/db/blastdb/")

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- FALSE
doMultiAmp <- FALSE
doTax <- TRUE

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
                      truncLen=c(200,180),
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
if(!exists("sample.data")){
    source("/SAN/Susanas_den/gitProj/HMHZ/bin/1_data_prep.R")
}
##Adding sample data
MAsample <- addSampleData(MA, sample.data)
##To phyloseq
##Sample data
PS <- toPhyloseq(MAsample, colnames(MAsample)) ##Now it work
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
MAsample <- addSampleData(MA, sample.data)
##To phyloseq
##Sample data
PS <- toPhyloseq(MAsample, colnames(MAsample)) ##Now it work
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi11.Rds")
saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi12.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi21.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi22.Rds")
sum(otu_table(PS)) ##Total denoised reads
### this is giving errors
##Primer data
PS.l <- toPhyloseq(MAsample, colnames(MAsample),  multi2Single=FALSE)
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
MAsample <- addSampleData(MA, sample.data)
##To phyloseq
##Sample data
PS <- toPhyloseq(MAsample, colnames(MAsample)) ##Now it work
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi11.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi12.Rds")
saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi21.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi22.Rds")
sum(otu_table(PS)) ##Total denoised reads
### this is giving errors
##Primer data
PS.l <- toPhyloseq(MAsample, colnames(MAsample),  multi2Single=FALSE)
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
                      truncLen=c(200,180),
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
MAsample <- addSampleData(MA, sample.data)
##To phyloseq
##Sample data
PS <- toPhyloseq(MAsample, colnames(MAsample)) ##Now it work
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi11.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi12.Rds")
#saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi21.Rds")
saveRDS(PS, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi22.Rds")
sum(otu_table(PS)) ##Total denoised reads
### this is giving errors
##Primer data
PS.l <- toPhyloseq(MAsample, colnames(MAsample),  multi2Single=FALSE)
###For primer analysis (Victor)
#saveRDS(PS.l11, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist11.Rds")
#saveRDS(PS.l12, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist12.Rds")
#saveRDS(PS.l21, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist21.Rds")
saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist22.Rds")

