ptable <- read.csv(file = "/SAN/Susanas_den/HMHZ/data/primer_list.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)


# amplicon targets
primerL <- read.table("data/primerInputUnique.csv", head=T, sep=",")
ptable$Primer_name <- paste(ptable$Name_F, ptable$Name_R, sep=".")
which(!ptable$Primer_name %in% primerL$Primer_name)
# manual correction
primerL$Primer_name[6] <- "18S_0067a_deg_3Mod_53_F.NSR399_3Mod_53_R"
primerL$Primer_name[7] <-  "18S_0067a_deg_5Mod_52_F.NSR399_5Mod_52_R" 
primerL$Primer_name[120] <- "Bgf_132_F.Bgr_132_R"
primerL$Primer_name[86] <- "NLF184cw_74_F.NL818cw_74_R"
p.df <- primerL[which(primerL$Primer_name%in%ptable$Primer_name),]

p.df$Primer_name[p.df$Gen=="28S"]
