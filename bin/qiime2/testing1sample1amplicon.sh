# import qiime2 artifact for one sample and one primer pair

qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/reads \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt \
      --output-path /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/demux-paired-end-S018-OUW.qza

# join reads
qiime vsearch join-pairs \
      --i-demultiplexed-seqs /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/demux-paired-end-S018-OUW.qza \
      --o-joined-sequences /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/demux-paired-end-joined-S018-OUW.qza

#make summary of reads
qiime demux summarize   \
      --i-data /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/demux-paired-end-joined-S018-OUW.qza \
      --o-visualization /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/demux-paired-end-joined-S018-OUW.qzv

#qualily filter
qiime quality-filter q-score \
      --i-demux  /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/demux-paired-end-joined-S018-OUW.qza \
      --o-filtered-sequences  /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/demux-paired-end-joined-S018-OUW-Filtered.qza \
      --o-filter-stats  /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/demux-paired-end-joined-S018-OUW-Filtered-stats.qza

# visualize first quality filter
qiime metadata tabulate \
      --m-input-file /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/demux-paired-end-joined-S018-OUW-Filtered-stats.qza \
      --o-visualization /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/demux-paired-end-joined-S018-OUW-filtered-stats.qzv

qiime deblur denoise-other \
      --i-reference-seqs /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURef_NR99_tax_silva.qza \
      --i-demultiplexed-seqs /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/demux-paired-end-joined-S018-OUW-Filtered.qza \
      --p-trim-length 420 \
      --p-sample-stats \
      --p-jobs-to-start 10 \
      --o-representative-sequences /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/req-seqs.qza \
      --o-table /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/table.qza \
      --o-stats /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/qiime2_sorted/18S_0067a_deg_3Mod_53_FNSR399_3Mod_53_R/deblur-stats.qza
