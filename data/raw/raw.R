qime2 codes

# on the server (qiime2-amplicon-2025.4)
cd /data

# make working directory + folders
mkdir -p project2_gastric_cancer_qiime2/{0_inputs,1_import,2_dada2,3_taxonomy,4_filtering,5_tree_rarefaction}
cd /data/project2_gastric_cancer_qiime2

# copy inputs
cp /datasets/project_2/gastric_cancer/gastric_cancer_manifest.tsv 0_inputs/manifest_raw.tsv
cp /datasets/project_2/gastric_cancer/gastric_cancer_metadata.tsv 0_inputs/metadata.tsv

# create single-end QIIME2 manifest (TSV) with direction=forward
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1 {print "sample-id","absolute-filepath","direction"; next}
{
  gsub("\r","",$1); gsub("\r","",$2);
  if($1!="" && $2!="") print $1,$2,"forward"
}' 0_inputs/manifest_raw.tsv > 0_inputs/manifest.tsv

# import into QIIME2
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path 0_inputs/manifest.tsv \
  --output-path 1_import/demux_seqs.qza

# summarize demux
qiime demux summarize \
  --i-data 1_import/demux_seqs.qza \
  --o-visualization 1_import/demux_summary.qzv



# DADA2 denoise (single-end; trunc-len=0 = no truncation)
# screen -S dada2_run2
# cd /data/project2_gastric_cancer_qiime2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs 1_import/demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --o-representative-sequences 2_dada2/rep-seqs.qza \
  --o-table 2_dada2/table.qza \
  --o-denoising-stats 2_dada2/stats.qza

# summarize feature table after denoising (sample counts + depth distribution)
qiime feature-table summarize \
  --i-table 2_dada2/table.qza \
  --o-visualization 2_dada2/table_summary.qzv \
  --m-sample-metadata-file 0_inputs/metadata.tsv

# view DADA2 denoising stats table
qiime metadata tabulate \
  --m-input-file 2_dada2/stats.qza \
  --o-visualization 2_dada2/denoising_stats.qzv

# view representative sequences (ASVs)
qiime feature-table tabulate-seqs \
  --i-data 2_dada2/rep-seqs.qza \
  --o-visualization 2_dada2/rep_seqs.qzv


