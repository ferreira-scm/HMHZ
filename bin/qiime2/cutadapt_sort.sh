#!/bin/bash

for p in /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/primer/*;
    do
	PRIMER=$(cat ${p})
	AMPLICON=$(echo ${p} | sed "s/\/SAN\/Susanas_den\/gitProj\/HMHZ\/tmp\/qiime2\/primer\///")
        mkdir /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_1_2/${AMPLICON}/
        for i in /SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_1_2/*_R1_001.fastq.gz;
do
    SAMPLE=$(echo ${i} | sed "s/_R1_\001\.fastq\.gz//")
    OUT=$(echo ${SAMPLE} | sed "s/\/SAN\/Susanas_den\/gitProj\/HMHZ\/data\/2018_22_HMHZ_1_2\///")
    echo ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_0.001.fastq.gz
    echo ${PRIMER}
    cutadapt ${PRIMER}\
	     --discard-untrimmed\
    -o /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_1_2/${AMPLICON}/${OUT}_R1_001.fastq.gz -p /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_1_2/${AMPLICON}/${OUT}_R2_001.fastq.gz ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz 
done    
done

for p in /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/primer/*;
    do
	PRIMER=$(cat ${p})
	AMPLICON=$(echo ${p} | sed "s/\/SAN\/Susanas_den\/gitProj\/HMHZ\/tmp\/qiime2\/primer\///")
        mkdir /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_2_1/${AMPLICON}/
        for i in /SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_2_1/*_R1_001.fastq.gz;
do
    SAMPLE=$(echo ${i} | sed "s/_R1_\001\.fastq\.gz//")
    OUT=$(echo ${SAMPLE} | sed "s/\/SAN\/Susanas_den\/gitProj\/HMHZ\/data\/2018_22_HMHZ_2_1\///")
    echo ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_0.001.fastq.gz
    echo ${PRIMER}
    cutadapt ${PRIMER}\
	     --discard-untrimmed\
    -o /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_2_1/${AMPLICON}/${OUT}_R1_001.fastq.gz -p /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_2_1/${AMPLICON}/${OUT}_R2_001.fastq.gz ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz 
done    
done

for p in /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/primer/*;
    do
	PRIMER=$(cat ${p})
	AMPLICON=$(echo ${p} | sed "s/\/SAN\/Susanas_den\/gitProj\/HMHZ\/tmp\/qiime2\/primer\///")
        mkdir /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_2_2/${AMPLICON}/
        for i in /SAN/Susanas_den/gitProj/HMHZ/data/2018_22_HMHZ_2_2/*_R1_001.fastq.gz;
do
    SAMPLE=$(echo ${i} | sed "s/_R1_\001\.fastq\.gz//")
    OUT=$(echo ${SAMPLE} | sed "s/\/SAN\/Susanas_den\/gitProj\/HMHZ\/data\/2018_22_HMHZ_2_2\///")
    echo ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_0.001.fastq.gz
    echo ${PRIMER}
    cutadapt ${PRIMER}\
	     --discard-untrimmed\
    -o /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_2_2/${AMPLICON}/${OUT}_R1_001.fastq.gz -p /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_2_2/${AMPLICON}/${OUT}_R2_001.fastq.gz ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz 
done    
done
