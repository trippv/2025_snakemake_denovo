#!/usr/bin/env Rscript

library(yaml)

### --- Abundance heatmap YAML ---

abund <- read.delim(snakemake@input[["abundance"]], header = TRUE, sep = "\t")

abund_data <- lapply(seq_len(nrow(abund)), function(i) {
  gene_values <- abund[i, -1]
  as.list(setNames(as.numeric(gene_values), colnames(gene_values)))
})
names(abund_data) <- abund$gene

abund_yaml <- list(
  id = "mqc_abundance",
  section_name = "Mapa de calor de abundancias",
  description = "Mapa de calor (heatmap) con los valores de abundancia de los  genes/transcritos expresados diferencialmente (valores normalizados)",
  plot_type = "heatmap",
  pconfig = list(
    id = "mqc_abundance_heatmap",
    title = "Abundancias normalizadas",
    xlab = "Muestras",
    ylab = "Genes"
  ),
  data = abund_data
)

writeLines(as.yaml(abund_yaml), snakemake@output[["abundance_yaml"]])

### --- PCA YAML ---

pca <- read.delim(snakemake@input[["pca"]], header = TRUE, sep = "\t")

pca_points <- lapply(seq_len(nrow(pca)), function(i) {
  list(
    x = round(pca$PC1[i], 4),
    y = round(pca$PC2[i], 4),
    color = pca$color[i]
  )
})
names(pca_points) <- pca$sample

pca_yaml <- list(
  id = "mqc_pca",
  section_name = "Análisis PCA",
  description = "Este gráfico muestra los resultados del Análisis de Componentes Principales (PCA) con valores normalizados (VSD) con colores por grupo",
  plot_type = "scatter",
  pconfig = list(
    id = "mqc_pca_scatter_plot",
    title = "PCA",
    ylab = "Componente 2",
    xlab = "Componente 1"
  ),
  data = pca_points
)

writeLines(as.yaml(pca_yaml), snakemake@output[["pca_yaml"]])


###-----volcano plot ---------

# Leer el archivo volcano.txt
volcano <- read.delim(snakemake@input[["volcano"]], header = TRUE, sep = "\t")

volvano_title <- volcano$contraste[1]

# Crear lista de puntos con coordenadas y color
volcano_points <- lapply(seq_len(nrow(volcano)), function(i) {
  list(
    x = round(volcano$logFC[i], 4),
    y = round(volcano$pval[i], 4),
    color = volcano$color[i]
  )
})
names(volcano_points) <- volcano$gene

# Crear el objeto completo
volcano_yaml <- list(
  id = "mqc_volcano",
  section_name = "Gráfico Volcano (Solo se muestra un contraste)",
  description = "Gráfico de volcan los resultados del análisis de expresión diferencial. Cada punto representa un gen, con el eje X mostrando el log2 Fold Change y el eje Y mostrando el -log10(p-valor). Los colores indican la significancia y dirección del cambio en la expresión génica.",
  plot_type = "scatter",
  pconfig = list(
    id = "mqc_volcano_plot",
    title = paste("Volcano plot -", volvano_title),
    xlab = "log2 Fold Change",
    ylab = "-log10(p-valor)"  # Cambia esto si estás usando -log10(p)
  ),
  data = volcano_points
)

# Escribir archivo YAML
writeLines(as.yaml(volcano_yaml), snakemake@output[["volcano_yaml"]])
