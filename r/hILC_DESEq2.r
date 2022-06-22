library("tximport")
library("readr")
library("tximportData")
library("tximeta")
library("DESeq2")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Use a meta data file to define the samples
samples <- read.table(file.path("/Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/RNAseq/meta/ILC3_meta.txt"), header = TRUE, sep="\t")
rownames(samples) <- samples$run

# Show the Samples
samples

# Get the file paths of all of the samples
files <- file.path("/Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/RNAseq/quants", samples$run, "quant.sf")

# Get the names of the samples
names(files) <- paste0(samples$run)

# Check to see if all files exist
all(file.exists(files))

files

# Change the row names to the sample names
rownames(samples) <- samples$run

# Import the tx to gene object
tx2gene <- read_csv(file.path("/Users/caz3so/scratch/20200620_hILC3_pchic/hg38_tx2gene.txt"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE, countsFromAbundance = "lengthScaledTPM")

# Ensure that all columns are the same order
all(colnames(txi$counts) == rownames(sample))

coldata <- samples
coldata$files <- files
coldata$names <- coldata$run

# Create the DESeq data set from Tximport data
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ CellType)
Gene_counts = counts(ddsTxi)

keep <- rowSums(counts(ddsTxi)) >= 10

ddsTxi <- ddsTxi[keep,]

Gene_counts = counts(ddsTxi)

# Plot the PCA for QC
plotPCA(vsd, intgroup=c("CellType", "Replicate"))

# Write the data to a file
write.table(Gene_counts, file="/Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/RNAseq/20200706_hILC3_CD4_GeneCounts.tsv", sep="\t", quote=F, col.names=NA)

ddsTxi <- DESeq(ddsTxi)

dds <- estimateSizeFactors(ddsTxi)

res <- results(ddsTxi)

vsd <- vst(ddsTxi, blind=FALSE)

vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)
pheatmap(vsd_cor, annotation = metadata)

normalized_counts <- counts(dds, normalized=TRUE)

write.table(vsd, file="/Users/caz3so/scratch/20200620_hILC3_pchic/20200611_function_map_pchic_2_features/20200622_hILC3_Colona_DESEQ2_vst.tsv", sep="\t", quote=F, col.names=NA)