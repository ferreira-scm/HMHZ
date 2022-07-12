library(phyloseq, lib="/usr/local/lib/R/site-library/")
library(mvabund)
library("iNEXT")

### multivariate analyses
## this takes WEEKS to run at the ASV level. I am starting at genus level
PS=readRDS(file="tmp/PSpre.R")
pPS=readRDS(file="tmp/pPSpre.R")

set.seed(1234)
PSrare <- rarefy_even_depth(pPS, sample.size=5000, replace=F)

rPSgen=tax_glom(PSrare, taxrank="genus")

rPSgen

rOTUmv <- mvabund(otu_table(rPSgen))

boxplot(otu_table(rPSgen), horizontal=TRUE, las=2)

#mean variance relationship
meanvar.plot(rOTUmv)

PSgen <- tax_glom(pPS, "genus")

# multivariate model
mod1 <- manyglm(rOTUmv ~
                        sample_data(rPSgen)$BMI +
                        sample_data(rPSgen)$Sex +
                        sample_data(rPSgen)$Transect +
                        sample_data(rPSgen)$Year +
                        sample_data(rPSgen)$eimASV +
                        sample_data(rPSgen)$criASV +
                        sample_data(rPSgen)$oxyASV +
                        sample_data(rPSgen)$hexASV +
                        sample_data(rPSgen)$hymASV +
                        sample_data(rPSgen)$triASV +
                        sample_data(rPSgen)$ascASV +
                        sample_data(rPSgen)$spiASV +
                        sample_data(rPSgen)$hekASV, family="negative_binomial")

plot(mod1)

#table1 <- summary.manyglm(mod1)
table2 <- anova.manyglm(mod1, p.uni = "adjusted", show.time="all")
#table1 <- summary.manyglm(mod1)
# this takes > 500h
sink("tmp/mlm_prel_rPSgen.txt")
table2
sink()
saveRDS(table2, file="tmp/manyanovatable2.RDS")
