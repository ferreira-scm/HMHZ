#select only the single 18S amplicons from the stratified samples
cp tmp/stratified_ma_Fullrun/*27M_F_98_F.Klin0341_CR_18_R* tmp/ma_sorted_3amplicon18S/27M_F_98_F.Klin0341_CR_18_R/.
cp tmp/stratified_ma_Fullrun/*wang1141_13_F.Nem_0425_6_3_R* tmp/ma_sorted_3amplicon18S/wang1141_13_F.Nem_0425_6_3_R/.
cp tmp/stratified_ma_Fullrun/*515F_Y_118_F.806R_118_R* tmp/ma_sorted_3amplicon18S/515F_Y_118_F.806R_118_R/.

#now we remove the primer names
for file in tmp/ma_sorted_3amplicon18S/wang1141_13_F.Nem_0425_6_3_R/*; do
    [ -f "$file" ] || continue
    mv "$file" "${file//wang1141_13_F.Nem_0425_6_3_R.fastq.gz/}";
done

for file in tmp/ma_sorted_3amplicon18S/27M_F_98_F.Klin0341_CR_18_R/*; do
    [ -f "$file" ] || continue
    mv "$file" "${file//27M_F_98_F.Klin0341_CR_18_R.fastq.gz/}";
done

for file in tmp/ma_sorted_3amplicon18S/515F_Y_118_F.806R_118_R/*; do
    [ -f "$file" ] || continue
    mv "$file" "${file//515F_Y_118_F.806R_118_R.fastq.gz/}";
done

# concatenate into one merged directory
cd tmp/ma_sorted_3amplicon18S/
cp wang1141_13_F.Nem_0425_6_3_R/* merged_single18S/.
cd 27M_F_98_F.Klin0341_CR_18_R
for sample in ./*; do
    cat "$sample" >> "../merged_single18S/$sample"
done
cd ../515F_Y_118_F.806R_118_R
for sample in ./*; do
    cat "$sample" >> "../merged_single18S/$sample"
done
cd /SAN/Susanas_den/gitProj/LabQiime2/sorted/
#qiime2 starts here

#sample seem to be empty  too smaller reads and that makes errors when runnung vsearch
find . -name "*.gz" -type 'f' -size -2k -exec mv "{}" ../empty_samples18Smerged/ \;
mv tmp/ma_sorted_3amplicon18S/merged_single18S/S039-GKS_S39_L001_R2_001.fastq.gz tmp/ma_sorted_3amplicon18S/empty_samples18Smerged
mv tmp/ma_sorted_3amplicon18S/merged_single18S/S269-NEGATIVE-22_S269_L001_R1_001.fastq.gz tmp/ma_sorted_3amplicon18S/empty_samples18Smerged

#import artifact into qiime2 format
qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/ma_sorted_3amplicon18S/merged_single18S/ \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt \
      --output-path /SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/demux-18S.qza


#join read pairs
qiime vsearch join-pairs \
     --i-demultiplexed-seqs tmp/demux-18S.qza \
     --o-joined-sequences tmp/demux-18S-joined.qza

qiime demux summarize   \
      --i-data tmp/demux-18S-joined.qza \
      --o-visualization tmp/demux-18S-joined.qzv

# denoise with dada2
qiime dada2 denoise-paired \
      --i-demultiplexed-seqs tmp/demux-18S.qza \
      --p-trunc-len-f 200 \
      --p-trunc-len-r 200 \
      --p-n-threads 0 \
      --o-table tmp/table-dada2-18S.qza \
      --o-representative-sequences tmp/rep-seqs-dada2-18S.qza \
      --o-denoising-stats tmp/denoising-dada2-18S-stats.qza

# denoise with deblur
qiime quality-filter q-score \
      --i-demux tmp/demux-18S-joined.qza \
      --o-filtered-sequences  tmp/demux-18S-FILTERED.qza \
      --o-filter-stats  tmp/demux-18S-filter-STATS.qza

qiime deblur denoise-other \
      --i-reference-seqs /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURef_NR99_tax_silva.qza \
      --i-demultiplexed-seqs tmp/demux-18S-FILTERED.qza \
      --p-trim-length 200 \
      --p-sample-stats \
      --p-jobs-to-start 20 \
      --o-representative-sequences tmp/req-seqs-deblur-18S.qza \
      --o-table tmp/table-deblur-18S.qza \
      --o-stats tmp/deblur-stats-18S.qza

qiime vsearch cluster-features-de-novo \
      --i-table tmp/table-deblur-18S.qza \
      --i-sequences tmp/req-seqs-deblur-18S.qza \
      --p-perc-identity 0.99 \
      --o-clustered-table tmp/table-dn-99-deblur-18S.qza \
      --o-clustered-sequences tmp/rep-seqs-dn-99-deblur-18S.qza

qiime vsearch uchime-denovo \
      --i-table tmp/table-dn-99-deblur-18S.qza \
      --i-sequences tmp/rep-seqs-dn-99-deblur-18S.qza \
      --output-dir tmp/uchime-dn-out-deblur-18S

qiime feature-table filter-features \
      --i-table tmp/table-dn-99-deblur-18S.qza \
      --m-metadata-file tmp/uchime-dn-out-deblur-18S/nonchimeras.qza \
      --o-filtered-table tmp/uchime-dn-out-deblur-18S/table-nonchimeric-wo-borderline.qza
qiime feature-table filter-seqs \
      --i-data tmp/rep-seqs-dn-99-deblur-18S.qza \
      --m-metadata-file tmp/uchime-dn-out-deblur-18S/nonchimeras.qza \
      --o-filtered-data tmp/uchime-dn-out-deblur-18S/rep-seqs-nonchimeric-wo-borderline.qza

qiime feature-classifier classify-sklearn \
      --i-classifier /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/classifier-SILVA-138.1.qza \
      --i-reads tmp//uchime-dn-out-deblur-18S/rep-seqs-nonchimeric-wo-borderline.qza \
      --p-n-jobs -2 \
      --o-classification tmp/taxonomy-deblur-18S.qza

#export qiime to phyloseq
qiime tools export \
      --input-path tmp/uchime-dn-out-deblur-18S/table-nonchimeric-wo-borderline.qza \
      --output-path phyloseq

biom convert \
     -i phyloseq/feature-table.biom \
     -o phyloseq/otu_table.txt \
     --to-tsv

qiime tools export \
      --input-path tmp/taxonomy.qza \
      --output-path phyloseq/
