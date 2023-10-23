#!/bin/bash

# demultiplexing and importing using manifest 

qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path ./anemia_manifest_updated.txt \
  --output-path /data/anemia/demux_seqs_anemia.qza

qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path ./infant_manifest.txt \
  --output-path /data/anemia/demux_seqs_infant.qza


  # Create visualization of demultiplexed samples
qiime demux summarize \
  --i-data demux_seqs_anemia.qza \
  --o-visualization demux_seqs_anemia.qzv

qiime demux summarize \
  --i-data demux_seqs_infant.qza \
  --o-visualization demux_seqs_infant.qzv



# Determine ASVs with DADA2
# 247 is the shortest lengeth read, quality score is solid throughout
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_seqs_anemia.qza \
  --p-trim-left 0 \
  --p-trunc-len 247 \
  --o-representative-sequences rep-seqs_anemia.qza \
  --o-table table_anemia.qza \
  --o-denoising-stats stats_anemia.qza

# quality score is solid throughout
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_seqs_infant.qza \
  --p-trim-left 0 \
  --p-trunc-len 150 \
  --o-representative-sequences rep-seqs_infant.qza \
  --o-table table_infant.qza \
  --o-denoising-stats stats_infant.qza

# Visualize ASVs stats
# Anemia
qiime feature-table summarize \
  --i-table table_anemia.qza \
  --o-visualization table_anemia.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/anemia/anemia_metadata.txt
  
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_anemia.qza \
  --o-visualization rep-seqs_anemia.qzv

# Infant
qiime feature-table summarize \
  --i-table table_infant.qza \
  --o-visualization table_infant.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/infant/infant_metadata.txt
  
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_infant.qza \
  --o-visualization rep-seqs_infant.qzv


# filter the anemia table using the metadata, keeping only the 6 month non-anemic samples
# and filter for only children feed with bm
qiime feature-table filter-samples \
  --i-table table_anemia.qza \
  --m-metadata-file /mnt/datasets/project_2/anemia/anemia_metadata.txt \
  --p-where "age_months='6' AND anemia='normal' AND NOT diet IN('OtherMilk','OtherMilk.Solids', 'Missing: Not collected', 'Solids.OtherLiquids')" \
  --o-filtered-table filtered_table_anemia.qza

# filter the infant table using the metadata, keeping only the 6 month non-anemic samples
# and filter for only children feed with bm
qiime feature-table filter-samples \
  --i-table table_infant.qza \
  --m-metadata-file /mnt/datasets/project_2/infant/infant_metadata.txt \
  --p-where "age_category='6 months' AND feed IN ('breast', 'combined')" \
  --o-filtered-table filtered_table_infant.qza


# convert table.qza into readable .qzv files
# anemia
qiime feature-table summarize \
  --i-table filtered_table_anemia.qza \
  --o-visualization filtered_table_anemia.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/anemia/anemia_metadata.txt

# infant
qiime feature-table summarize \
  --i-table filtered_table_infant.qza \
  --o-visualization filtered_table_infant.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/infant/infant_metadata.txt

qiime metadata tabulate \
  --m-input-file stats_anemia.qza \
  --o-visualization stats_anemia.qzv

# Infant
# Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs_infant.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree_infant.qza \
  --o-rooted-tree rooted-tree_infant.qza 

# Alpha-rarefaction, max-depth set at 3rd quartile of frequency per sample
qiime diversity alpha-rarefaction \
  --i-table filtered_table_infant.qza \
  --i-phylogeny rooted-tree_infant.qza \
  --p-max-depth 42049 \
  --m-metadata-file /mnt/datasets/project_2/infant/infant_metadata.txt \
  --o-visualization alpha-rarefaction_infant.qzv

# Anemia
# Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs_anemia.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree_anemia.qza \
  --o-rooted-tree rooted-tree_anemia.qza 

# Alpha-rarefaction
qiime diversity alpha-rarefaction \
  --i-table filtered_table_anemia.qza \
  --i-phylogeny rooted-tree_anemia.qza \
  --p-max-depth 50000 \
  --m-metadata-file /mnt/datasets/project_2/anemia/anemia_metadata.txt \
  --o-visualization alpha-rarefaction_anemia.qzv