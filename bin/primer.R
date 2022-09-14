primer1 <- read.table("data/primer_list.csv", head=T, sep=",")

primer2 <- read.table("data/primerInputUnique.csv", head=T, sep=",")

primer1$Primer_name <- paste(primer1$Name_F, primer1$Name_R, sep=".")

primer2$seq <- paste(primer2$Seq_F, primer2$Seq_R)

primer2$seq

primer <- merge(primer1, primer2, by="seq")

primer <- unique(primer)

nrow(primer)

primer

