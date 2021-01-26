#!/bin/bash

# first we use data sorted and primers removed with cutadapt

cd /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2

cp sortedReads_1_1/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/* merged/.

for primer in /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/sortedReads_1_1/*; do
    cd $primer
    for sample in *; do
	cat "$sample" >> "/SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged/$sample"
    done
done

cd /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged/

#sample seem to be empty  too smaller reads and that makes errors when runnung vsearch
mkdir ../emptySamplesMerged
find . -name "*.gz" -type 'f' -size -10k -exec mv "{}" ../emptySamplesMerged/ \;

#import artifact into qiime2 format
qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/merged/ \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt \
      --output-path /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/demux-pooled.qza

cd /SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/

#join read pairs
qiime vsearch join-pairs \
     --i-demultiplexed-seqs demux-pooled.qza \
     --o-joined-sequences demux-pooled-joined.qza

qiime demux summarize   \
      --i-data demux-pooled-joined.qza \
      --o-visualization demux-pooled-joined.qzv

# denoise with deblur
qiime quality-filter q-score \
      --i-demux demux-pooled-joined.qza \
      --o-filtered-sequences  demux-pooled-FILTERED.qza \
      --o-filter-stats  demux-pooled-filter-STATS.qza

qiime deblur denoise-other --i-reference-seqs /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURef_NR99_tax_silva.qza --i-demultiplexed-seqs demux-pooled-FILTERED.qza --p-trim-length 300 --p-sample-stats --p-jobs-to-start 20 --o-representative-sequences rep-seqs-deblur-pooled.qza --o-table table-deblur-pooled.qza --o-stats deblur-stats-pooled.qza


qiime deblur denoise-other --i-reference-seqs /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURef_NR99_tax_silva.qza --i-demultiplexed-seqs demux-pooled-FILTERED.qza --p-trim-length 200 --p-sample-stats --p-jobs-to-start 20 --o-representative-sequences rep-seqs200-deblur-pooled.qza --o-table table200-deblur-pooled.qza --o-stats deblur200-stats-pooled.qza


#qiime vsearch cluster-features-de-novo \
#      --i-table tmp/table-deblur-18S.qza \
#      --i-sequences tmp/req-seqs-deblur-18S.qza \
#      --p-perc-identity 0.99 \
#      --o-clustered-table tmp/table-dn-99-deblur-18S.qza \
#      --o-clustered-sequences tmp/rep-seqs-dn-99-deblur-18S.qza

# remove chimeras denovo

qiime vsearch uchime-denovo --i-table table-deblur-pooled.qza --i-sequences rep-seqs-deblur-pooled.qza --output-dir uchime-dn-out-deblur-pooled

qiime vsearch uchime-denovo --i-table table200-deblur-pooled.qza --i-sequences rep-seqs200-deblur-pooled.qza --output-dir uchime-dn-out-deblur200-pooled


qiime feature-table filter-features --i-table table-deblur-pooled.qza --m-metadata-file uchime-dn-out-deblur-pooled/nonchimeras.qza --o-filtered-table uchime-dn-out-deblur-pooled/table-nonchimeric-wo-borderline.qza
qiime feature-table filter-seqs --i-data rep-seqs-deblur-pooled.qza --m-metadata-file uchime-dn-out-deblur-pooled/nonchimeras.qza --o-filtered-data uchime-dn-out-deblur-pooled/rep-seqs-nonchimeric-wo-borderline.qza

qiime feature-table filter-features --i-table table200-deblur-pooled.qza --m-metadata-file uchime-dn-out-deblur200-pooled/nonchimeras.qza --o-filtered-table uchime-dn-out-deblur200-pooled/table-nonchimeric-wo-borderline.qza
qiime feature-table filter-seqs --i-data rep-seqs200-deblur-pooled.qza --m-metadata-file uchime-dn-out-deblur200-pooled/nonchimeras.qza --o-filtered-data uchime-dn-out-deblur200-pooled/rep-seqs-nonchimeric-wo-borderline.qza


# taxonomic annotation
qiime feature-classifier classify-sklearn --i-classifier /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/classifier-SILVA-138.1.qza --i-reads uchime-dn-out-deblur-pooled/rep-seqs-nonchimeric-wo-borderline.qza --p-n-jobs -2 --o-classification taxonomy-deblur-pooled.qza

qiime feature-classifier classify-sklearn --i-classifier /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/classifier-SILVA-138.1.qza --i-reads uchime-dn-out-deblur200-pooled/rep-seqs-nonchimeric-wo-borderline.qza --p-n-jobs -2 --o-classification taxonomy-deblur200-pooled.qza

#export qiime to phyloseq
qiime tools export --input-path uchime-dn-out-deblur-pooled/table-nonchimeric-wo-borderline.qza --output-path phyloseq_pooled
biom convert -i phyloseq_pooled/feature-table.biom -o phyloseq_pooled/otu_table.txt --to-tsv
qiime tools export --input-path taxonomy-deblur-pooled.qza --output-path phyloseq_pooled/

qiime tools export --input-path uchime-dn-out-deblur200-pooled/table-nonchimeric-wo-borderline.qza --output-path phyloseq_pooled200
biom convert -i phyloseq_pooled200/feature-table.biom -o phyloseq_pooled200/otu_table.txt --to-tsv
qiime tools export --input-path taxonomy-deblur200-pooled.qza --output-path phyloseq_pooled200/


# Manually change #OTUID to OTUID fin otu_table.txt
# Manually change "feature ID" to "OTUID" in taxonomy-pooled.tsv

