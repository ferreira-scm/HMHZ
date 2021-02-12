#!/bin/bash

# first we use data sorted and primers removed with cutadapt

cd /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2

### need to change the names directly on the sorted reads directories

#for primer in /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_2_2/*; do
#    cd $primer
#    for file in *; do
#    newfile="$(echo ${file} | sed "s/^./S/")";
#    mv "${file}" "${newfile}";
#    done
#done

# being lazy here. I want to create an fmpty file for each sample
# I know the samples in run 2 start with capital letter, and that is a  problem

cp sortedReads_1_1/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/* merged11/.
cp sortedReads_1_2/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/* merged12/.
cp sortedReads_2_1/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/* merged21/.
cp sortedReads_2_2/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/* merged22/.

for file in merged11/*; do
    > $file
done
for file in merged12/*; do
    > $file
done

for file in merged21/*; do
    > $file
done
for file in merged22/*; do
    > $file
done



# now we concatenate all the samples from all the runs
for primer in /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_1_1/*; do
    cd $primer
    for sample in *; do
	cat "$sample" >> "/SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged11/$sample"
    done
done

for primer in /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_1_2/*; do
    cd $primer
    for sample in *; do
	cat "$sample" >> "/SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged12/$sample"
    done
done

for primer in /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_2_1/*; do
    cd $primer
    for sample in *; do
	cat "$sample" >> "/SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged21/$sample"
    done
done

for primer in /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_2_2/*; do
    cd $primer
    for sample in *; do
	cat "$sample" >> "/SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged22/$sample"
    done
done

cd /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged11/

#sample seem to be empty  too smaller reads and that makes errors when runnung vsearch
mkdir ../emptySamplesMerged
cd /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged11/
find . -name "*.gz" -type 'f' -size -10k -exec mv "{}" ../emptySamplesMerged/ \;
cd /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged12/
find . -name "*.gz" -type 'f' -size -10k -exec mv "{}" ../emptySamplesMerged/ \;
cd /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged21/
find . -name "*.gz" -type 'f' -size -10k -exec mv "{}" ../emptySamplesMerged/ \;
cd /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged22/
find . -name "*.gz" -type 'f' -size -10k -exec mv "{}" ../emptySamplesMerged/ \;
# some  samples have no pair
cd /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/

mv merged12/s310-SK-3135_S310_L001_R2_001.fastq.gz emptySamplesMerged/
mv merged12/s336-AA-0137_S336_L001_R2_001.fastq.gz  emptySamplesMerged/
mv merged12/s334-NE-H021_S334_L001_R2_001.fastq.gz emptySamplesMerged/

mv merged21/S561-NE-H040_S264_L001_R1_001.fastq.gz  emptySamplesMerged/
mv merged21/S318-SK-3094_S21_L001_R2_001.fastq.gz emptySamplesMerged/
mv merged21/S369-AA-0416_S72_L001_R2_001.fastq.gz emptySamplesMerged/
mv merged21/S377-AA-0151_S80_L001_R2_001.fastq.gz emptySamplesMerged/
mv merged21/S372-AA-0100_S75_L001_R2_001.fastq.gz emptySamplesMerged/

#import artifact into qiime2 format
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged11/ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/demux-pooled11.qza

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged12/ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/demux-pooled12.qza

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged21/ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/demux-pooled21.qza

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged22/ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/demux-pooled22.qza

cd /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/

#join read pairs
qiime vsearch join-pairs --i-demultiplexed-seqs demux-pooled11.qza --o-joined-sequences demux-pooled-joined11.qza
qiime demux summarize --i-data demux-pooled-joined11.qza --o-visualization demux-pooled-joined11.qzv

qiime vsearch join-pairs --i-demultiplexed-seqs demux-pooled12.qza --o-joined-sequences demux-pooled-joined12.qza
qiime demux summarize --i-data demux-pooled-joined12.qza --o-visualization demux-pooled-joined12.qzv

qiime vsearch join-pairs --i-demultiplexed-seqs demux-pooled21.qza --o-joined-sequences demux-pooled-joined21.qza
qiime demux summarize --i-data demux-pooled-joined21.qza --o-visualization demux-pooled-joined21.qzv

qiime vsearch join-pairs --i-demultiplexed-seqs demux-pooled22.qza --o-joined-sequences demux-pooled-joined22.qza
qiime demux summarize --i-data demux-pooled-joined22.qza --o-visualization demux-pooled-joined22.qzv


# denoise with deblur
qiime quality-filter q-score --i-demux demux-pooled-joined11.qza --o-filtered-sequences  demux-pooled-FILTERED11.qza --o-filter-stats  demux-pooled-filter-STATS11.qza

qiime quality-filter q-score --i-demux demux-pooled-joined12.qza --o-filtered-sequences  demux-pooled-FILTERED12.qza --o-filter-stats  demux-pooled-filter-STATS12.qza

qiime quality-filter q-score --i-demux demux-pooled-joined21.qza --o-filtered-sequences  demux-pooled-FILTERED21.qza --o-filter-stats  demux-pooled-filter-STATS21.qza

qiime quality-filter q-score --i-demux demux-pooled-joined22.qza --o-filtered-sequences  demux-pooled-FILTERED22.qza --o-filter-stats  demux-pooled-filter-STATS22.qza


qiime deblur denoise-other --i-reference-seqs /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURef_NR99_tax_silva.qza --i-demultiplexed-seqs demux-pooled-FILTERED11.qza --p-trim-length 200 --p-sample-stats --p-jobs-to-start 20 --o-representative-sequences rep-seqs-deblur-pooled11.qza --o-table table-deblur-pooled11.qza --o-stats deblur-stats-pooled11.qza

qiime deblur denoise-other --i-reference-seqs /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURef_NR99_tax_silva.qza --i-demultiplexed-seqs demux-pooled-FILTERED12.qza --p-trim-length 200 --p-sample-stats --p-jobs-to-start 20 --o-representative-sequences rep-seqs-deblur-pooled12.qza --o-table table-deblur-pooled12.qza --o-stats deblur-stats-pooled12.qza

qiime deblur denoise-other --i-reference-seqs /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURef_NR99_tax_silva.qza --i-demultiplexed-seqs demux-pooled-FILTERED21.qza --p-trim-length 200 --p-sample-stats --p-jobs-to-start 20 --o-representative-sequences rep-seqs-deblur-pooled21.qza --o-table table-deblur-pooled21.qza --o-stats deblur-stats-pooled21.qza

qiime deblur denoise-other --i-reference-seqs /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURef_NR99_tax_silva.qza --i-demultiplexed-seqs demux-pooled-FILTERED22.qza --p-trim-length 200 --p-sample-stats --p-jobs-to-start 20 --o-representative-sequences rep-seqs-deblur-pooled22.qza --o-table table-deblur-pooled22.qza --o-stats deblur-stats-pooled22.qza



#qiime vsearch cluster-features-de-novo \
#      --i-table tmp/table-deblur-18S.qza \
#      --i-sequences tmp/req-seqs-deblur-18S.qza \
#      --p-perc-identity 0.99 \
#      --o-clustered-table tmp/table-dn-99-deblur-18S.qza \
#      --o-clustered-sequences tmp/rep-seqs-dn-99-deblur-18S.qza

# remove chimeras denovo

qiime vsearch uchime-denovo --i-table table-deblur-pooled11.qza --i-sequences rep-seqs-deblur-pooled11.qza --output-dir uchime-dn-out-deblur-pooled11
qiime vsearch uchime-denovo --i-table table-deblur-pooled12.qza --i-sequences rep-seqs-deblur-pooled12.qza --output-dir uchime-dn-out-deblur-pooled12
qiime vsearch uchime-denovo --i-table table-deblur-pooled21.qza --i-sequences rep-seqs-deblur-pooled21.qza --output-dir uchime-dn-out-deblur-pooled21
qiime vsearch uchime-denovo --i-table table-deblur-pooled22.qza --i-sequences rep-seqs-deblur-pooled22.qza --output-dir uchime-dn-out-deblur-pooled22


qiime feature-table filter-features --i-table table-deblur-pooled11.qza --m-metadata-file uchime-dn-out-deblur-pooled11/nonchimeras.qza --o-filtered-table uchime-dn-out-deblur-pooled11/table-nonchimeric-wo-borderline.qza
qiime feature-table filter-seqs --i-data rep-seqs-deblur-pooled11.qza --m-metadata-file uchime-dn-out-deblur-pooled11/nonchimeras.qza --o-filtered-data uchime-dn-out-deblur-pooled11/rep-seqs-nonchimeric-wo-borderline.qza
qiime feature-table filter-features --i-table table-deblur-pooled12.qza --m-metadata-file uchime-dn-out-deblur-pooled12/nonchimeras.qza --o-filtered-table uchime-dn-out-deblur-pooled12/table-nonchimeric-wo-borderline.qza
qiime feature-table filter-seqs --i-data rep-seqs-deblur-pooled12.qza --m-metadata-file uchime-dn-out-deblur-pooled12/nonchimeras.qza --o-filtered-data uchime-dn-out-deblur-pooled12/rep-seqs-nonchimeric-wo-borderline.qza
qiime feature-table filter-features --i-table table-deblur-pooled21.qza --m-metadata-file uchime-dn-out-deblur-pooled21/nonchimeras.qza --o-filtered-table uchime-dn-out-deblur-pooled21/table-nonchimeric-wo-borderline.qza
qiime feature-table filter-seqs --i-data rep-seqs-deblur-pooled21.qza --m-metadata-file uchime-dn-out-deblur-pooled21/nonchimeras.qza --o-filtered-data uchime-dn-out-deblur-pooled21/rep-seqs-nonchimeric-wo-borderline.qza
qiime feature-table filter-features --i-table table-deblur-pooled22.qza --m-metadata-file uchime-dn-out-deblur-pooled22/nonchimeras.qza --o-filtered-table uchime-dn-out-deblur-pooled22/table-nonchimeric-wo-borderline.qza
qiime feature-table filter-seqs --i-data rep-seqs-deblur-pooled22.qza --m-metadata-file uchime-dn-out-deblur-pooled22/nonchimeras.qza --o-filtered-data uchime-dn-out-deblur-pooled22/rep-seqs-nonchimeric-wo-borderline.qza


# taxonomic annotation
qiime feature-classifier classify-sklearn --i-classifier /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/classifier-SILVA-138.1.qza --i-reads uchime-dn-out-deblur-pooled11/rep-seqs-nonchimeric-wo-borderline.qza --p-n-jobs -2 --o-classification taxonomy-deblur-pooled11.qza
qiime feature-classifier classify-sklearn --i-classifier /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/classifier-SILVA-138.1.qza --i-reads uchime-dn-out-deblur-pooled12/rep-seqs-nonchimeric-wo-borderline.qza --p-n-jobs -2 --o-classification taxonomy-deblur-pooled12.qza
qiime feature-classifier classify-sklearn --i-classifier /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/classifier-SILVA-138.1.qza --i-reads uchime-dn-out-deblur-pooled21/rep-seqs-nonchimeric-wo-borderline.qza --p-n-jobs -2 --o-classification taxonomy-deblur-pooled21.qza
qiime feature-classifier classify-sklearn --i-classifier /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/classifier-SILVA-138.1.qza --i-reads uchime-dn-out-deblur-pooled22/rep-seqs-nonchimeric-wo-borderline.qza --p-n-jobs -2 --o-classification taxonomy-deblur-pooled22.qza


#export qiime to phyloseq
qiime tools export --input-path uchime-dn-out-deblur-pooled11/table-nonchimeric-wo-borderline.qza --output-path phyloseq_pooled11
biom convert -i phyloseq_pooled11/feature-table.biom -o phyloseq_pooled11/otu_table.txt --to-tsv
qiime tools export --input-path taxonomy-deblur-pooled11.qza --output-path phyloseq_pooled11/

qiime tools export --input-path uchime-dn-out-deblur-pooled12/table-nonchimeric-wo-borderline.qza --output-path phyloseq_pooled12
biom convert -i phyloseq_pooled12/feature-table.biom -o phyloseq_pooled12/otu_table.txt --to-tsv
qiime tools export --input-path taxonomy-deblur-pooled12.qza --output-path phyloseq_pooled12/

qiime tools export --input-path uchime-dn-out-deblur-pooled21/table-nonchimeric-wo-borderline.qza --output-path phyloseq_pooled21
biom convert -i phyloseq_pooled21/feature-table.biom -o phyloseq_pooled21/otu_table.txt --to-tsv
qiime tools export --input-path taxonomy-deblur-pooled21.qza --output-path phyloseq_pooled21/

qiime tools export --input-path uchime-dn-out-deblur-pooled22/table-nonchimeric-wo-borderline.qza --output-path phyloseq_pooled22
biom convert -i phyloseq_pooled22/feature-table.biom -o phyloseq_pooled22/otu_table.txt --to-tsv
qiime tools export --input-path taxonomy-deblur-pooled22.qza --output-path phyloseq_pooled22/


# Manually change #OTUID to OTUID fin otu_table.txt
# Manually change "feature ID" to "OTUID" in taxonomy-pooled.tsv

