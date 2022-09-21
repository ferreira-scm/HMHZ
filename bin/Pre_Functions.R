fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # abundance filtering at 0.001%
    keepTaxa = (x / sum(x) > 0.00001)
#    keepTaxa = (x / sum(x) > 0.0005)
    summary(keepTaxa)
    ps = phyloseq::prune_taxa(keepTaxa, ps)
# plus prevalnce filter at 1%
    KeepTaxap <- microbiome::prevalence(ps)>0.01
    ps <- phyloseq::prune_taxa(KeepTaxap, ps)
# subset samples based on total read count (500 reads)
#ps <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 500)
    ps <- phyloseq::prune_samples(sample_sums(ps)>500, ps)
    ps
}


############## old preproc

##Eliminate Unassigned to superkingdom level
#PS <- subset_taxa(PS, !is.na(superkingdom) & !superkingdom %in% c("", "uncharacterized"))

# subset samples based on total read count (1000 reads)
median(phyloseq::sample_sums(PS))

#hist(phyloseq::sample_sums(PS))

#PS <- phyloseq::subset_samples(PS, phyloseq::sample_sums(PS) > 1000)

# prevalence filtering at 5%
#pPPS=phyloseq_filter_prevalence(pPS, prev.trh=0.05)

###A lot of Mus :(
## Host read numbers
#sum(otu_table(subset_taxa(PS, genus%in%"Mus")))/sum(otu_table(PS))
###Eliminate reads assigned as "Mus"
#PS <- subset_taxa(PS, !genus %in% "Mus") ##Eliminate reads :S

# Eliminate samples with no reads
#PS <- prune_samples(sample_sums(PS)>0, PS)

# abundance filtering to 0.01%? Or keep prevalence filtering?
#x = taxa_sums(PS)
#keepTaxa = (x / sum(x) > 0.0001)
#summary(keepTaxa)
#pPS = prune_taxa(keepTaxa, PS)

# Eliminate samples with no reads
#pPS <- prune_samples(sample_sums(pPS)>0, pPS)
