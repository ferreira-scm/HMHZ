# merge all the runs
PS11 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi11.Rds")
PS12 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSCombi12.Rds")
PS21 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi21.Rds")
PS22 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSCombi22.Rds")
PS <- merge_phyloseq(PS11, PS12, PS21, PS22)

PSl11 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist11.Rds")
PSl12 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run1/PSlist12.Rds")
PSl21 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist21.Rds")
PSl22 <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/run2/PSlist22.Rds")

###Eliminate the empty primers first and then merge the phyloseq lists... empty list make the next function bug
PSl11[sapply(PSl11, is.null)]<- NULL
PSl12[sapply(PSl12, is.null)]<- NULL
PSl21[sapply(PSl21, is.null)]<- NULL
PSl22[sapply(PSl22, is.null)]<- NULL
 ##Merge all the information from both experiments
along<- names(PSl11)
PSl <- lapply(along, function(i) merge_phyloseq(PSl11[[i]], PSl12[[i]],PSl21[[i]], PSl22[[i]]))
names(PSl) <- names(PSl11)

# adjusting amplicon names for some reason
x<- names(PSl)
x[6]<-"BGf_132_F.BGr_132_R"
x[22]<-"LSU_Fwd_2_3Mod_55_F.LSU_Rev_4_54_R"
x[28]<-"NLF184cw_74_F.NL818cw_74_R"
names(PSl)<- x

saveRDS(PSl, file="tmp/PSl.Rds")
saveRDS(PS, file="tmp/PSCombi.Rds")

PSl <- readRDS(file="tmp/PSl.Rds")

# load qiime-derived PS

PSq <- readRDS(file="/SAN/Susanas_den/gitProj/HMHZ/tmp/PSqiime200.RDS")
PSq
PS

rm(PS11)
rm(PS12)
rm(PS21)
rm(PS22)
rm(PSl11)
rm(PSl12)
rm(PSl21)
rm(PSl22)
rm(along)

##Eliminate Unassigned to superkingdom level
PS <- subset_taxa(PS, !is.na(superkingdom) & !superkingdom %in% c("", "uncharacterized"))

PSl <- lapply(PSl, function(x){
    subset_taxa(x, !is.na(superkingdom) & !superkingdom %in% c("", "uncharacterized"))})

    
# subset samples based on total read count (1000 reads)
median(phyloseq::sample_sums(PS))

PS <- phyloseq::subset_samples(PS, phyloseq::sample_sums(PS) > 1000)

# can't do it for the PSl
#lapply(PSl, function(x){
#    subset_samples(x, phyloseq::sample_sums(x)>50
#})

