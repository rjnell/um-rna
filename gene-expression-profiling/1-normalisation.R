###  
### STEP 1 - DESEQ2 on TCGA RNA-seq
###

# Set seed
set.seed(83539)

# Load fragments_per_gene data from tcga
tcga_fpg_path = "data/tcga.fragments_per_gene.tsv"
tcga_fpg = read.csv(tcga_fpg_path, stringsAsFactors = F, check.names = F, sep = "\t", row.names = 1)

# Plot PCA of raw data
tcga_pca = prcomp(t(tcga_fpg))
plot(tcga_pca$x, pch=16)

# Load samplesheet
samplesheet = read.csv("data/tcga-samplesheet.tsv", sep="\t", stringsAsFactors = F)

# Replace sample_IDs by sample_names
cases = rownames(tcga_pca$x)
cases = stringr::str_remove(cases,"s_")
cases = stringr::str_replace_all(cases,"_","-")
names = samplesheet$File.Name
names = stringr::str_remove(names,"_gdc_realn_rehead.bam")
cases = samplesheet$Case.ID[match(cases, names)]
colnames(tcga_fpg) = cases

# Install DESeq2
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# install.packages("colorspace")
# install.packages("tibble")

# Load DESeq2
library(DESeq2)

# Build samples_table
samples_table = matrix(nrow = 80, ncol = 1)
rownames(samples_table) = cases
colnames(samples_table) = "cohort"

# Build DESeqData
dds = DESeqDataSetFromMatrix(tcga_fpg,
                             colData = samples_table,
                             design = ~ 1)

# Only keep ENSG genes
keep = which(substr(rownames(dds),1,4) == "ENSG") 
dds = dds[keep,]

# Normalise and save
tcga_normalised = assay(vst(dds, blind=F))
write.table(tcga_normalised, "gene-expression-profiling/tcga-normalised.tsv", sep="\t")
