


TMPtoPhyloseq <- function(MA, samples, multi2Single=TRUE, ...){
    STL <- getSequenceTableNoChime(MA)
    TTL <- getTaxonTable(MA, simplify=FALSE)
    if(length(TTL) == nrow(MA)){
        TAX <- TRUE
    } else if(length(TTL) == 0){
        message("No taxon table provided, so your phyloseq object will lack",
                " taxonomical annotations")
        TAX <- FALSE
    } else {
        stop("\nTaxon tables in provided in Multiamplicon object for are",
             " incongruent with the number of amplicons")
    }
    if(multi2Single){
        ## get sample tables filled with zeros for non-assessed
        ## samples for particular amplicons
        filledST <- lapply(STL, .fillSampleTables, samples=samples)
        allST <- as.matrix(Reduce(cbind, filledST))
        if (!all(dim(allST)>0)) {
            stop(paste("\nempty OTU table\n",
                       "rownames",  unlist(lapply(STL, base::rownames)),
                       "don't match sample names:",  samples, "\n"))
        }
        ## The same for taxon annotations
        if(TAX){
            all.tax <- as.matrix(Reduce(rbind, TTL))
            ## to avoid problems with duplicated rownames (same sequences
            ## recovered for different amplicons), this can happen after trimming
            base::rownames(all.tax) <- make.unique(base::rownames(all.tax))
        }
        base::colnames(allST) <- make.unique(base::colnames(allST))
        ## wrap it up into one Phyloseq object
        phyloseq(otu_table(allST, taxa_are_rows=FALSE),
                 sample_data(MA@sampleData),
                 if (TAX) tax_table(all.tax),
                 ...)
    } else {
        PS.l <- lapply(seq_along(STL), function (i) {
            ## currently taxa tables are NULL if empty and
            ## sequence Tables have zero dimensions
            seqExists <- all(dim(STL[[i]])>0)
            if(TAX) {
                taxExists <- !is.null(TTL[[i]])
            } else {
                taxExists <- FALSE
            }
            if(seqExists) {
                allSampleTable <- .fillSampleTables(STL[[i]], samples=samples)
                phyloseq(otu_table(allSampleTable, taxa_are_rows=FALSE),
                         if(TAX && taxExists) tax_table(TTL[[i]]),
                         sample_data(MA@sampleData[base::rownames(allSampleTable),]),
                         ...)
            } else if(!isTRUE(seqExists) && !isTRUE(taxExists)){
                NULL
            } else  {
                stop(paste("inconsistent taxa data provided for",
                           "sequences in amplicon", i,
                           rownames(MA)[[i]], "\n")
                     )
            }
        })
        ## somehow can't use rownames(MA)
        names(PS.l) <- rownames(MA)
        PS.l
    }
}
         
