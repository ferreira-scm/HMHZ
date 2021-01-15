

primerL <- read.csv("data/primer_list.csv")
primerl <- read.csv("data/primerInputUnique.csv")

primerL$seq <- paste(primerL$Seq_F,primerL$Seq_R, sep=".")

primerl$seq <- paste(primerl$Seq_F,primerl$Seq_R, sep=".")

head(primerL)

pri <- merge(primerL, primerl, by="seq")

pri$Primer_comb_ID
