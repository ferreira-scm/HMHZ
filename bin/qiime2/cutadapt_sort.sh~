#!/bin/bash

for p in /SAN/Susanas_den/EimeriaMicrobiome/primer/*;
    do
	PRIMER=$(cat ${p})
	AMPLICON=$(echo ${p} | sed "s/\/SAN\/Susanas_den\/EimeriaMicrobiome\/primer\///")
        mkdir /SAN/Susanas_den/EimeriaMicrobiome/sortedReads/${AMPLICON}/
        for i in /SAN/Susanas_den/EimeriaMicrobiome/data/2018_22_Eie_FullRun_1/*_R1_001.fastq.gz;
do
    SAMPLE=$(echo ${i} | sed "s/_R1_\001\.fastq\.gz//")
    OUT=$(echo ${SAMPLE} | sed "s/\/SAN\/Susanas_den\/EimeriaMicrobiome\/data\/2018_22_Eie_FullRun_1\///")
#    echo ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_0.001.fastq.gz
    #    echo ${PRIMER}
    cutadapt ${PRIMER}\
	     --discard-untrimmed\
    -o /SAN/Susanas_den/EimeriaMicrobiome/sortedReads/${AMPLICON}/${OUT}_R1_001.fastq.gz -p /SAN/Susanas_den/EimeriaMicrobiome/sortedReads/${AMPLICON}/${OUT}_R2_001.fastq.gz ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz 
done    
done
