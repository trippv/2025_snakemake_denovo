#!/usr/bin/env Rscript

library(DESeq2)
library(ggplot2)
library(readr)
library(tibble)
library(RColorBrewer)

#### 1: Cargar datos ####
samples <- read.table(snakemake@input[["metadata"]], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(samples) <- samples$sample

counts <- read.csv(snakemake@input[["counts"]], header = TRUE, row.names = 1)
counts <- counts[, samples$sample, drop = FALSE]  # Asegura orden y correspondencia

# Redondear los valores de cuentas para DESeq2
counts <- round(counts, 0)

###2. Filtrar genes con poca representación
keep <- rowSums(counts >= 10) >= 2
counts_filtered <- counts[keep, ]

###3. Crear objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = samples, design = ~ group)
vsd <- varianceStabilizingTransformation(dds)

# ====== Análisis exploratorio ======

# Matriz de distancias
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- samples$sample
colnames(sampleDistMatrix) <- samples$sample
sampleDistMatrix <- rownames_to_column(as.data.frame(sampleDistMatrix), "sample")



# Guardar matriz de distancias en results/summary/qc
write.table(sampleDistMatrix, snakemake@output[["dist_matrix"]], 
            row.names = FALSE, quote = FALSE, sep = "\t")

### PCA =========================================================

# ========== Paleta de colores por grupo ==========
groups <- unique(samples$group)
n_groups <- length(groups)
colors <- setNames(brewer.pal(min(max(n_groups, 3), 8), "Set1")[1:n_groups], groups)


# Análisis de PCA usando vsd (mas rapido que log transformation)
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
pcaData$color <- colors[as.character(pcaData$group)]


pcaData_df <- data.frame(
  sample = pcaData$name,
  PC1 = pcaData$PC1,
  PC2 = pcaData$PC2,
  group = pcaData$group,
  color = pcaData$color
)

# Guardar datos de PCA en results/summary_qc
write.table(pcaData_df, snakemake@output[["pca"]],
            row.names = FALSE, quote = FALSE, sep = "\t")

# ====== Expresión diferencial =================================
# Leer umbrales desde el config de Snakemake
fdr_limit <- snakemake@config[["diff_exp"]][["fdr_threshold"]]
lfc_limit <- snakemake@config[["diff_exp"]][["lfc_threshold"]]


################Expresion por grupos pareados####
group_names <- unique(samples$group)
#generar todas las combinaciones de pares de grupos
comparisons <- combn(group_names, 2)

dds <- DESeq(dds)

# Crear directorio de salida si Snakemake no lo hizo
dir.create(snakemake@output[["de_dir"]], recursive = TRUE, showWarnings = FALSE)


# Vector para recolectar genes significativos de TODOS los contrastes
all_sig_genes <- c()

# Loop para realizar el análisis para cada par de grupos
for (i in 1:ncol(comparisons)) {
  g1 <- comparisons[1, i]
  g2 <- comparisons[2, i] 

# 1. Obtener resultados del contraste específico
  res_name <- paste0(g1, "_vs_", g2)
  res_curr <- results(dds, contrast = c("group", g1, g2))
  res_curr_df <- as.data.frame(res_curr)
  res_curr_df$gene <- rownames(res_curr_df)

# 2. Identificar significativos con los umbrales del config
  is_sig <- which(res_curr_df$padj < fdr_limit & abs(res_curr_df$log2FoldChange) >= lfc_limit)
  res_curr_df$significant <- "No significativo"
  res_curr_df$significant[is_sig] <- "Significativo"

# 3. Guardar el archivo individual de la comparación
  write.table(res_curr_df, 
              file = file.path(snakemake@output[["de_dir"]], 
              paste0("DE_", res_name, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)


# 4. Acumular los IDs de genes significativos para el heatmap global
  sig_genes_curr <- res_curr_df$gene[is_sig]
  all_sig_genes <- union(all_sig_genes, sig_genes_curr)
}

# Filtrar para obtener solo los 200 más significativos de la bolsa total
# Para esto, podemos re-extraer los p-valores mínimos de cualquiera de los resultados
# o simplemente tomar los primeros 200 si la lista es muy larga.
if (length(all_sig_genes) > 200) {
    # Una forma elegante: tomar los 200 con el padj promedio más bajo o el mínimo
    # Aquí, por simplicidad, tomamos los primeros 200 recolectados
    all_sig_genes <- all_sig_genes[1:200]
}

# Generar la matriz de abundancia con esta lista maestra
counts_de <- counts[all_sig_genes, ]
counts_de_mx <- as.matrix(counts_de)
counts_de_mx <- t(scale(t(counts_de_mx), center = TRUE, scale = TRUE))
counts_de_df <- rownames_to_column(as.data.frame(counts_de_mx), "gene")

write.table(counts_de_df, snakemake@output[["abundance"]], 
            row.names = FALSE, quote = FALSE, sep = "\t")


##### Volcano plot
# Usar el último contraste como ejemplo 
# si res_curr_df tiene mas de 800 filas, tomar solo los significativos

if(nrow(res_curr_df) > 800){
  res_df_volcano <- res_curr_df[res_curr_df$significant == "Significativo", ]
} else {
  res_df_volcano <- res_curr_df
} 


# Preparar tabla para volcano plot
res_df_volcano <- data.frame(
  gene = rownames(res_df_volcano),
  logFC = res_df_volcano$log2FoldChange,
  pval = -log(res_df_volcano$padj),
  significant = res_df_volcano$significant,
  contraste = res_name
)

# El color del volcano ahora respeta el lfc_limit del config
res_df_volcano$color <- ifelse(
  res_df_volcano$significant != "Significativo", "#b3aaaa",
  ifelse(res_df_volcano$logFC >= lfc_limit, "red", "blue")
)
# Guardar tabla para volcano plot en results/summary_qc
write.table(res_df_volcano, snakemake@output[["volcano"]], 
            row.names = FALSE, quote = FALSE, sep = "\t")

