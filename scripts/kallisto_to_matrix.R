#!/usr/bin/env Rscript
library(tximport)
library(readr)
library(tibble)

# Read the sample table
sample_table <- read_tsv(snakemake@input[["table"]], col_names = c("sample", "file"))

# Extract file paths
files <- sample_table$file
names(files) <- sample_table$sample

# Read the gene-to-transcript map
tx2gene <- read_tsv(snakemake@input[["gene_trans_map"]], col_names = c("GENEID", "TXNAME"))
tx2gene <- tx2gene[, c("TXNAME", "GENEID")]  # <- reordenar para tximport

# Import Kallisto results with tximport for gene-level counts
txi_gene <- tximport(files, type = "kallisto", txIn = TRUE, txOut = FALSE, tx2gene = tx2gene)

# Import Kallisto results with tximport for transcript-level TPM
txi_trans <- tximport(files, type = "kallisto", txIn = TRUE, txOut = TRUE)

# Save gene-level counts matrix
counts <- txi_gene$counts
colnames(counts) <- sample_table$sample
counts <- as.data.frame(counts)
counts <- rownames_to_column(counts, "gene_id")

write.csv(counts, snakemake@output[["gene_matrix"]], row.names = FALSE, quote = FALSE)


# Save transcript-level TPM matrix
tpm <- txi_trans$counts # counts a nivel de transcrito
colnames(tpm) <- sample_table$sample
tpm <- as.data.frame(tpm)
tpm <- rownames_to_column(tpm, "transcript_id")
write.csv(tpm, snakemake@output[["transcript_matrix"]], row.names = FALSE, quote = FALSE)