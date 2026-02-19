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
