###Load sample information

sample.data <- read.csv("/SAN/Susanas_den/gitProj/HMHZ/data/Sample_selection_Metabarcoding_Complete.csv", dec=",", stringsAsFactors=FALSE)

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
