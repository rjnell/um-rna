### Analyze how samples differ in their usage of splice junctions
### Input: SJ.out.tab files derived from STAR alignment (GRCh38)
### Output: table + images *****

# Cohort 1 SJ.out.tab files
dir_cohort1 = "F:/project-235-RNAseqPieterUvealMelanoma/analysis/11_TCGA/samples/"
splicing_files_cohort1 = list.files(path = dir_cohort1, pattern="SJ.out.tab", recursive = T)

# Merge cohorts
merged_files = c(paste0(dir_cohort1, splicing_files_cohort1))
merged_names = c(splicing_files_cohort1)

# Initialize the list of all splicing_loci
all_splicing_loci = NULL

# Initialize the list of all sample_names
sample_names = NULL

library(data.table)

# Iterate through all files
for (f in 1:length(merged_files)) {
  splicing_data = fread(merged_files[f], sep="\t", header = F, stringsAsFactors = F)
  splicing_loci = paste0(splicing_data$V1, "-", splicing_data$V2, "-", splicing_data$V3)
  all_splicing_loci = unique(c(all_splicing_loci, splicing_loci))
  sample_names = c(sample_names, strsplit(merged_names[f], "/")[[1]][1])
}

# Create a matrix
splicing_loci_usage = matrix(data = 0, nrow = length(all_splicing_loci), ncol = length(merged_files))
rownames(splicing_loci_usage) = all_splicing_loci
sample_names = substr(merged_names,3,38)
samplesheet = readxl::read_xlsx("data/data.xlsx")
colnames(splicing_loci_usage) = samplesheet$`Patient ID`[match(sample_names, stringr::str_replace_all(substr(samplesheet$`File Name`,0,36),"-","_"))]

# Iterate through all files
for (f in 1:length(merged_files)) {
  splicing_data = fread(merged_files[f], sep="\t", header = F, stringsAsFactors = F)
  splicing_loci = paste0(splicing_data$V1, "-", splicing_data$V2, "-", splicing_data$V3)
  splicing_loci_usage[splicing_loci, f] = splicing_data$V7
}

splicing_loci_usage[1:10,1:20]

###

# PCA on unnormalized data
combined_pca = prcomp(t(splicing_loci_usage))
plot(combined_pca$x, pch=16)
text(combined_pca$x, labels=rownames(combined_pca$x), pos=2, cex=0.5)

# Load DESeq2
library(DESeq2)

# Build samples_table
samples_table = as.matrix(factor(rep("cohort_1",80)),,drop=F)
rownames(samples_table) = colnames(splicing_loci_usage)
colnames(samples_table) = "cohort"

# Build DESeqData
DESeqData = DESeqDataSetFromMatrix(splicing_loci_usage,
                                   colData = samples_table,
                                   design = ~ 1)

# Perform vst, a variance stabilizing transformation 
vsd = vst(DESeqData, blind=F)
combined_normalized = assay(vsd)
combined_normalized_genes = rownames(combined_normalized)

# Illustrate outcome
head(combined_normalized, 3)

# Build PCA of combined normalized data
combined_normalized_pca = prcomp(t(combined_normalized))
plot(combined_normalized_pca$x, col=c(rep("#BE1E2D",35), rep("#2EAADE",75)), pch=16)
text(combined_normalized_pca$x, labels=rownames(combined_normalized_pca$x), pos=4, cex=0.5)


write.table(combined_normalized, "tcga-splicing-normalized.tsv", row.names = T, sep="\t")
saveRDS(combined_normalized, "tcga-splicing-normalized.RDS")

