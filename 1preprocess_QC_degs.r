# ================================================================================
# GSE209998 - Análisis Integrado: QC + DEGs con datos UQN
# Cáncer de Mama: Primary vs Metastatic
# Optimizado para construcción posterior de redes de coexpresión (GCN)
# ================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(RColorBrewer)
  library(pheatmap)
  library(ggrepel)
  library(data.table)
  library(cluster)
  library(vegan)
  library(tibble)
  library(GEOquery)
  library(limma)
  library(matrixStats)
  library(tibble)
})

set.seed(123)

# Configuración de directorios
setwd("C:/Users/marcr/OneDrive/Escritorio/DATASETS/NCBI/dataset_breast_UQN")
dir.create("./results", showWarnings = FALSE)
dir.create("./results/QC", showWarnings = FALSE)
dir.create("./results/DEGs", showWarnings = FALSE)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║  ANÁLISIS GSE209998: PRIMARY vs METASTATIC BREAST CANCER      ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# ================================================================================
# 1) CARGA DE DATOS Y METADATA
# ================================================================================

cat("═══ PASO 1: CARGA DE DATOS ═══\n")

# 1.1) Obtener metadata desde GEO
gse <- getGEO("GSE209998", GSEMatrix = TRUE)
eset <- gse[[1]]
pheno <- pData(eset)

# Extraer condición
cond <- rep(NA_character_, nrow(pheno))
cond[grepl("primary", pheno$source_name_ch1, ignore.case = TRUE)] <- "Primary"
cond[grepl("metast", pheno$source_name_ch1, ignore.case = TRUE)] <- "Metastatic"

cat("  ✓ Condiciones detectadas:\n")
print(table(cond, useNA = "ifany"))

gsm_ids <- pheno$geo_accession

# 1.2) Leer matriz UQN
uqn <- read.table(
  "GSE209998_UQN.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = NULL,
  comment.char = "",
  stringsAsFactors = FALSE  # Importante para evitar conversiones
)

gene_ids <- as.character(uqn[[1]])
uqn_mat <- as.matrix(uqn[, -1, drop = FALSE])

# CRÍTICO: Manejar duplicados causados por Excel
# Detectar genes convertidos a fechas (patrón: número-Mes)
date_pattern <- grepl("^\\d+-[A-Za-z]{3}$", gene_ids)
if (any(date_pattern)) {
  cat("  ⚠ ADVERTENCIA: Detectados", sum(date_pattern), 
      "genes convertidos a fechas por Excel\n")
  cat("    Ejemplos:", paste(head(unique(gene_ids[date_pattern]), 5), collapse = ", "), "\n")
}

# Hacer rownames únicos añadiendo sufijos a duplicados
gene_ids_unique <- make.unique(gene_ids, sep = "_dup")
n_duplicates <- sum(gene_ids != gene_ids_unique)
if (n_duplicates > 0) {
  cat("  ⚠ CORRECCIÓN:", n_duplicates, "genes duplicados renombrados\n")
}

rownames(uqn_mat) <- gene_ids_unique
colnames(uqn_mat) <- gsm_ids

cat("  ✓ Dimensiones originales:", nrow(uqn_mat), "genes x", ncol(uqn_mat), "muestras\n")
cat("  ✓ Genes únicos:", length(unique(gene_ids_unique)), "\n")

# 1.3) Construir metadata
meta <- data.frame(
  sample = gsm_ids,
  condition = cond,
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$sample

# Filtrar muestras sin condición
keep_samp <- !is.na(meta$condition)
meta <- meta[keep_samp, , drop = FALSE]
uqn_mat <- uqn_mat[, keep_samp, drop = FALSE]

# Renombrar muestras
primary_idx <- which(meta$condition == "Primary")
metastatic_idx <- which(meta$condition == "Metastatic")

new_names <- character(nrow(meta))
new_names[primary_idx] <- paste0("Primary_", seq_along(primary_idx))
new_names[metastatic_idx] <- paste0("Metastatic_", seq_along(metastatic_idx))

meta$sample_renamed <- new_names
rownames(meta) <- new_names
colnames(uqn_mat) <- new_names

cat("  ✓ Muestras finales:\n")
cat("    - Primary:", length(primary_idx), "\n")
cat("    - Metastatic:", length(metastatic_idx), "\n")
cat("    - Total:", ncol(uqn_mat), "\n\n")

# ================================================================================
# 2) TRANSFORMACIÓN LOGARÍTMICA (CRÍTICO PARA LIMMA)
# ================================================================================

cat("═══ PASO 2: TRANSFORMACIÓN LOGARÍTMICA ═══\n")

# IMPORTANTE: Los datos UQN están normalizados pero NO en escala log
# Para limma necesitamos log2(UQN + offset)
offset <- 1  # Evita log(0)

expr_log <- log2(uqn_mat + offset)

cat("  ✓ Transformación aplicada: log2(UQN + 1)\n")
cat("  ✓ Rango valores originales: [", round(min(uqn_mat), 2), ",", 
    round(max(uqn_mat), 2), "]\n")
cat("  ✓ Rango valores log2: [", round(min(expr_log), 2), ",", 
    round(max(expr_log), 2), "]\n\n")

# ================================================================================
# 3) FILTRADO DE GENES
# ================================================================================

cat("═══ PASO 3: FILTRADO DE GENES ═══\n")

# 3.1) Eliminar genes con varianza cero o NA
vars <- rowVars(expr_log, na.rm = TRUE)
keep_var <- !is.na(vars) & vars > 0
expr_filt <- expr_log[keep_var, , drop = FALSE]

cat("  ✓ Genes eliminados (var=0 o NA):", sum(!keep_var), "\n")

# 3.2) Filtrado por expresión mínima
# Mantener genes con expresión > percentil 25 en al menos 25% de muestras
min_samples <- ceiling(0.25 * ncol(expr_filt))
expr_threshold <- quantile(expr_filt, 0.25)
keep_expr <- rowSums(expr_filt > expr_threshold) >= min_samples
expr_filt <- expr_filt[keep_expr, , drop = FALSE]

cat("  ✓ Genes eliminados (baja expresión):", sum(!keep_expr), "\n")
cat("  ✓ Genes retenidos:", nrow(expr_filt), "\n\n")

# ================================================================================
# 4) CONTROL DE CALIDAD
# ================================================================================

cat("═══ PASO 4: CONTROL DE CALIDAD ═══\n")

group <- factor(meta$condition, levels = c("Primary", "Metastatic"))
names(group) <- rownames(meta)
colors_group <- c("Primary" = "#4ECDC4", "Metastatic" = "#FF6B6B")

# 4.1) Distribución de expresión
pdf("./results/QC/01_expression_distribution.pdf", width = 14, height = 6)
par(mfrow = c(1, 2), mar = c(10, 4, 3, 1))

# Boxplot
boxplot(expr_filt, las = 2, outline = FALSE, 
        col = colors_group[group],
        main = "Distribución log2(UQN+1) por muestra", 
        ylab = "log2(UQN + 1)", 
        cex.axis = 0.7, cex.names = 0.7)
legend("topleft", legend = names(colors_group), 
       fill = colors_group, bty = "n")

# Densidades
plot(density(expr_filt[, 1]), main = "Densidad de expresión", 
     xlab = "log2(UQN + 1)", lwd = 2, ylim = c(0, max(density(expr_filt[, 1])$y)))
for (i in 2:ncol(expr_filt)) {
  lines(density(expr_filt[, i]), lwd = 1, 
        col = colors_group[group[i]])
}
legend("topright", legend = names(colors_group), 
       col = colors_group, lwd = 2, bty = "n")
dev.off()

cat("  ✓ Gráfico 1: Distribución de expresión\n")

# 4.2) RLE (Relative Log Expression)
gene_medians <- rowMedians(expr_filt, na.rm = TRUE)
rle_mat <- sweep(expr_filt, 1, gene_medians, FUN = "-")
rle_summary <- data.frame(
  sample = colnames(rle_mat),
  median_RLE = apply(rle_mat, 2, median, na.rm = TRUE),
  IQR_RLE = apply(rle_mat, 2, IQR, na.rm = TRUE),
  condition = group
)
write.csv(rle_summary, "./results/QC/RLE_summary.csv", row.names = FALSE)

p_rle <- ggplot(rle_summary, aes(x = sample, y = median_RLE, color = condition)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "RLE: Mediana por muestra", 
       y = "Mediana RLE", x = "") +
  scale_color_manual(values = colors_group)

ggsave("./results/QC/02_RLE_median.png", p_rle, width = 12, height = 5, dpi = 300)
cat("  ✓ Gráfico 2: RLE analysis\n")

# 4.3) Identificación de HVGs (Highly Variable Genes)
# CRÍTICO: Estos genes se usarán para GCN según recomendación del supervisor
n_hvgs <- 5000
vars_all <- rowVars(expr_filt, na.rm = TRUE)
hvgs <- names(sort(vars_all, decreasing = TRUE))[1:min(n_hvgs, nrow(expr_filt))]
expr_hvg <- expr_filt[hvgs, ]

cat("  ✓ HVGs seleccionados para análisis:", length(hvgs), "\n")

# Guardar lista de HVGs
write.csv(
  data.frame(
    gene = hvgs,
    variance = vars_all[hvgs],
    mean_expr = rowMeans(expr_filt[hvgs, ])
  ),
  "./results/QC/HVGs_top5000.csv",
  row.names = FALSE
)

# 4.4) PCA con HVGs
pca <- prcomp(t(expr_hvg), scale. = TRUE, center = TRUE)
pca_var <- (pca$sdev^2) / sum(pca$sdev^2)
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = group,
  sample = colnames(expr_hvg)
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 2.8, max.overlaps = 30) +
  theme_minimal(base_size = 13) +
  labs(
    title = sprintf("PCA en %d HVGs", length(hvgs)),
    x = sprintf("PC1 (%.1f%%)", 100 * pca_var[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * pca_var[2])
  ) +
  scale_color_manual(values = colors_group)

ggsave("./results/QC/03_PCA_HVGs.png", p_pca, width = 10, height = 7, dpi = 300)
cat("  ✓ Gráfico 3: PCA\n")


# 4.5) MDS con Spearman y Pearson (todas las muestras)

# Spearman
cor_s <- cor(expr_filt, method = "spearman")
dist_s <- as.dist(1 - cor_s)
mds_s  <- cmdscale(dist_s, k = 2)

mds_df_s <- data.frame(
  Dim1      = mds_s[, 1],
  Dim2      = mds_s[, 2],
  condition = group,
  sample    = colnames(expr_filt)
)

p_mds_s <- ggplot(mds_df_s,
                  aes(Dim1, Dim2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 2.8, max.overlaps = 30) +
  theme_minimal(base_size = 13) +
  labs(title = "MDS (1 - Spearman)") +
  scale_color_manual(values = colors_group)

ggsave("./results/QC/04_MDS_Spearman_ALL.png",
       p_mds_s, width = 10, height = 7, dpi = 300)

# Pearson
cor_p <- cor(expr_filt, method = "pearson")
dist_p <- as.dist(1 - cor_p)
mds_p  <- cmdscale(dist_p, k = 2)

mds_df_p <- data.frame(
  Dim1      = mds_p[, 1],
  Dim2      = mds_p[, 2],
  condition = group,
  sample    = colnames(expr_filt)
)

p_mds_p <- ggplot(mds_df_p,
                  aes(Dim1, Dim2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 2.8, max.overlaps = 30) +
  theme_minimal(base_size = 13) +
  labs(title = "MDS (1 - Pearson)") +
  scale_color_manual(values = colors_group)

ggsave("./results/QC/04_MDS_Pearson_ALL.png",
       p_mds_p, width = 10, height = 7, dpi = 300)

cat("  ✓ Gráficos 4: MDS Spearman/Pearson (todas las muestras)\n")


# 4.6) PERMANOVA (sobre Spearman con todas las muestras)

set.seed(123)
adonis_result <- adonis2(dist_s ~ condition, data = meta, permutations = 999)
capture.output(adonis_result,
               file = "./results/QC/PERMANOVA_results_ALL.txt")
cat("  ✓ PERMANOVA completado (todas las muestras, p =",
    round(adonis_result$`Pr(>F)`[1], 4), ")\n")


# 4.7) Heatmaps de correlación (todas las muestras)

ann_col_all <- data.frame(condition = group,
                          row.names = colnames(expr_filt))

pheatmap(cor_s,
         annotation_col   = ann_col_all,
         annotation_colors = list(condition = colors_group),
         main            = "Correlación de Spearman (ALL)",
         show_colnames   = FALSE,
         show_rownames   = FALSE,
         filename        = "./results/QC/05_correlation_heatmap_Spearman_ALL.png",
         width           = 8, height = 7)

pheatmap(cor_p,
         annotation_col   = ann_col_all,
         annotation_colors = list(condition = colors_group),
         main            = "Correlación de Pearson (ALL)",
         show_colnames   = FALSE,
         show_rownames   = FALSE,
         filename        = "./results/QC/05_correlation_heatmap_Pearson_ALL.png",
         width           = 8, height = 7)

cat("  ✓ Gráficos 5: Heatmaps Spearman/Pearson (todas las muestras)\n\n")


# 4.8) Identificación y eliminación de muestras conflictivas
# (definidas a partir de los heatmaps/MDS previos)

samples_to_remove <- c(
  "Primary_24", "Primary_40", "Primary_19", "Primary_2",
  "Primary_36", "Primary_25", "Primary_15", "Primary_13",
  "Primary_29", "Primary_31", "Primary_44", "Primary_38",
  "Primary_16",
  "Metastatic_40", "Metastatic_20", "Metastatic_5", "Metastatic_16"
)

if (length(samples_to_remove) > 0) {
  write.csv(
    data.frame(sample = samples_to_remove),
    "./results/QC/samples_conflictive_manual_QC.csv",
    row.names = FALSE
  )
}

samples_to_keep <- setdiff(colnames(expr_filt), samples_to_remove)

expr_clean <- expr_filt[, samples_to_keep, drop = FALSE]
meta_clean <- meta[samples_to_keep, , drop = FALSE]

group_clean <- factor(meta_clean$condition, levels = c("Primary", "Metastatic"))
names(group_clean) <- rownames(meta_clean)


# 4.9) QC "CLEAN": Spearman y Pearson sin muestras conflictivas

# Spearman CLEAN
cor_s_clean  <- cor(expr_clean, method = "spearman")
dist_s_clean <- as.dist(1 - cor_s_clean)

set.seed(123)
adonis_result_clean <- adonis2(dist_s_clean ~ condition,
                               data = meta_clean,
                               permutations = 999)
capture.output(adonis_result_clean,
               file = "./results/QC/PERMANOVA_results_CLEAN.txt")

ann_col_clean <- data.frame(condition = group_clean,
                            row.names = colnames(expr_clean))

pheatmap(cor_s_clean,
         annotation_col   = ann_col_clean,
         annotation_colors = list(condition = colors_group),
         main            = "Correlación de Spearman (CLEAN)",
         show_colnames   = FALSE,
         show_rownames   = FALSE,
         filename        = "./results/QC/05_correlation_heatmap_Spearman_CLEAN.png",
         width           = 8, height = 7)

mds_s_clean <- cmdscale(dist_s_clean, k = 2)
mds_df_s_clean <- data.frame(
  Dim1      = mds_s_clean[, 1],
  Dim2      = mds_s_clean[, 2],
  condition = group_clean,
  sample    = colnames(expr_clean)
)

p_mds_s_clean <- ggplot(mds_df_s_clean,
                        aes(Dim1, Dim2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 2.8, max.overlaps = 30) +
  theme_minimal(base_size = 13) +
  labs(title = "MDS (1 - Spearman) - CLEAN") +
  scale_color_manual(values = colors_group)

ggsave("./results/QC/06_MDS_Spearman_CLEAN.png",
       p_mds_s_clean, width = 10, height = 7, dpi = 300)


# Pearson CLEAN

cor_p_clean  <- cor(expr_clean, method = "pearson")
dist_p_clean <- as.dist(1 - cor_p_clean)

pheatmap(cor_p_clean,
         annotation_col   = ann_col_clean,
         annotation_colors = list(condition = colors_group),
         main            = "Correlación de Pearson (CLEAN)",
         show_colnames   = FALSE,
         show_rownames   = FALSE,
         filename        = "./results/QC/05_correlation_heatmap_Pearson_CLEAN.png",
         width           = 8, height = 7)

mds_p_clean <- cmdscale(dist_p_clean, k = 2)
mds_df_p_clean <- data.frame(
  Dim1      = mds_p_clean[, 1],
  Dim2      = mds_p_clean[, 2],
  condition = group_clean,
  sample    = colnames(expr_clean)
)

p_mds_p_clean <- ggplot(mds_df_p_clean,
                        aes(Dim1, Dim2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 2.8, max.overlaps = 30) +
  theme_minimal(base_size = 13) +
  labs(title = "MDS (1 - Pearson) - CLEAN") +
  scale_color_manual(values = colors_group)

ggsave("./results/QC/06_MDS_Pearson_CLEAN.png",
       p_mds_p_clean, width = 10, height = 7, dpi = 300)

cat("  ✓ Versión CLEAN generada (Spearman/Pearson: heatmaps + MDS + PERMANOVA)\n\n")


# OPCIONAL: a partir de aquí usar expr_clean y meta_clean en el resto del pipeline
# Después de definir expr_clean / meta_clean en el QC:
group_clean <- factor(meta_clean$condition, levels = c("Primary", "Metastatic"))
names(group_clean) <- rownames(meta_clean)
# expr_filt <- expr_clean
meta      <- meta_clean
group     <- group_clean


# Después de crear expr_clean / meta_clean / group_clean

# 4.10) Recalcular HVGs sobre datos CLEAN
n_hvgs <- 5000
vars_all_clean <- rowVars(expr_clean, na.rm = TRUE)
hvgs <- names(sort(vars_all_clean, decreasing = TRUE))[1:min(n_hvgs, nrow(expr_clean))]
expr_hvg <- expr_clean[hvgs, ]

cat("  ✓ HVGs (CLEAN) seleccionados para análisis:", length(hvgs), "\n")

write.csv(
  data.frame(
    gene = hvgs,
    variance = vars_all_clean[hvgs],
    mean_expr = rowMeans(expr_clean[hvgs, ])
  ),
  "./results/QC/HVGs_top5000.csv",
  row.names = FALSE
)

# PCA con HVGs sobre datos CLEAN
pca <- prcomp(t(expr_hvg), scale. = TRUE, center = TRUE)
pca_var <- (pca$sdev^2) / sum(pca$sdev^2)
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = group_clean,
  sample = colnames(expr_hvg)
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 2.8, max.overlaps = 30) +
  theme_minimal(base_size = 13) +
  labs(
    title = sprintf("PCA en %d HVGs (CLEAN)", length(hvgs)),
    x = sprintf("PC1 (%.1f%%)", 100 * pca_var[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * pca_var[2])
  ) +
  scale_color_manual(values = colors_group)

ggsave("./results/QC/03_PCA_HVGs_CLEAN.png", p_pca, width = 10, height = 7, dpi = 300)
cat("  ✓ Gráfico 3: PCA (HVGs CLEAN)\n")

# === Exportar matriz de expresión CLEAN filtrada a HVGs (input para EnGNet) ===

# Asegurarse de que hvgs se ha calculado sobre expr_clean:
# vars_all_clean <- rowVars(expr_clean, na.rm = TRUE)
# hvgs <- names(sort(vars_all_clean, decreasing = TRUE))[1:min(5000, nrow(expr_clean))]
# expr_hvg <- expr_clean[hvgs, ]

# Guardar matriz de expresión: genes en filas, muestras en columnas
expr_hvg_df <- data.frame(
  gene = rownames(expr_hvg),
  expr_hvg,
  check.names = FALSE
)

write.table(
  expr_hvg_df,
  "./results/DEGs/expr_clean_HVGs5000.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("  ✓ Matriz expr_clean_HVGs5000.csv generada para EnGNet\n")



########## OBTENER CSV DOS CONDICIONES PRIMARY Y METASTASIC ##########
# Separar matrices por condición
primary_samples    <- rownames(meta_clean)[meta_clean$condition == "Primary"]
metastatic_samples <- rownames(meta_clean)[meta_clean$condition == "Metastatic"]

expr_primary    <- expr_clean[, primary_samples, drop = FALSE]
expr_metastatic <- expr_clean[, metastatic_samples, drop = FALSE]

cat("  ✓ Muestras Primary:", ncol(expr_primary), "\n")
cat("  ✓ Muestras Metastatic:", ncol(expr_metastatic), "\n")



# HVGs en Primary
vars_primary <- rowVars(expr_primary, na.rm = TRUE)
hvgs_primary <- names(sort(vars_primary, decreasing = TRUE))[1:min(n_hvgs, length(vars_primary))]
expr_hvg_primary <- expr_primary[hvgs_primary, ]

write.csv(
  data.frame(
    gene      = hvgs_primary,
    variance  = vars_primary[hvgs_primary],
    mean_expr = rowMeans(expr_primary[hvgs_primary, ])
  ),
  "./results/DEGs/HVGs_top5000_Primary.csv",
  row.names = FALSE
)

# HVGs en Metastatic
vars_metast <- rowVars(expr_metastatic, na.rm = TRUE)
hvgs_metast <- names(sort(vars_metast, decreasing = TRUE))[1:min(n_hvgs, length(vars_metast))]
expr_hvg_metast <- expr_metastatic[hvgs_metast, ]

write.csv(
  data.frame(
    gene      = hvgs_metast,
    variance  = vars_metast[hvgs_metast],
    mean_expr = rowMeans(expr_metastatic[hvgs_metast, ])
  ),
  "./results/DEGs/HVGs_top5000_Metastatic.csv",
  row.names = FALSE
)

# Exportar matrices de expresión HVGs por condición (input para EnGNet)

expr_hvg_primary_df <- data.frame(
  gene = rownames(expr_hvg_primary),
  expr_hvg_primary,
  check.names = FALSE
)

write.table(
  expr_hvg_primary_df,
  "./results/DEGs/expr_clean_HVGs5000_Primary.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

expr_hvg_metast_df <- data.frame(
  gene = rownames(expr_hvg_metast),
  expr_hvg_metast,
  check.names = FALSE
)

write.table(
  expr_hvg_metast_df,
  "./results/DEGs/expr_clean_HVGs5000_Metastatic.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("  ✓ Matrices para EnGNet generadas (Primary y Metastatic)\n")



# ================================================================================
# 5) ANÁLISIS DIFERENCIAL CON LIMMA
# ================================================================================

cat("═══ PASO 5: ANÁLISIS DIFERENCIAL (LIMMA) ═══\n")

# 5.1) Diseño experimental
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

cat("  ✓ Diseño experimental:\n")
print(head(design))

# 5.2) Ajuste del modelo lineal
# NOTA: Los datos ya están en log2, listos para limma
# CRÍTICO: Preservar rownames para evitar problemas posteriores
stopifnot(!is.null(rownames(expr_clean)))
stopifnot(length(rownames(expr_clean)) == nrow(expr_clean))

fit <- lmFit(expr_clean, design)

# Verificación de integridad
cat("  ✓ Genes en el modelo:", nrow(fit$coefficients), "\n")
cat("  ✓ Tipo de identificadores:", class(rownames(fit$coefficients)), "\n")
cat("  ✓ Ejemplos:", paste(head(rownames(fit$coefficients), 3), collapse = ", "), "\n")

# 5.3) Contraste
cont <- makeContrasts(
  Metastatic_vs_Primary = Metastatic - Primary,
  levels = design
)
fit2 <- contrasts.fit(fit, cont)

# ================================================================================
# 6) ESTRATEGIA DUAL: eBayes + TREAT
# ================================================================================

cat("\n═══ PASO 6: OBTENCIÓN DE DEGs ═══\n")

# 6.1) Método eBayes (para redes de coexpresión)
fit_ebayes <- eBayes(fit2)

# CRÍTICO: topTable reordena filas, debemos preservar los genes originales
deg_all <- topTable(fit_ebayes, number = Inf, adjust.method = "BH", 
                    coef = "Metastatic_vs_Primary", sort.by = "none")

# Asignar manualmente los rownames de expr_clean
rownames(deg_all) <- rownames(expr_clean)

# Ahora reordenar por p-valor para tener los más significativos primero
deg_all <- deg_all[order(deg_all$adj.P.Val), ]

# Verificación
cat("\n  ✓ Resultados eBayes obtenidos\n")
cat("    - Total genes:", nrow(deg_all), "\n")
cat("    - Primeros 3 genes:", paste(head(rownames(deg_all), 3), collapse = ", "), "\n")
cat("    - ¿Rownames únicos?:", length(unique(rownames(deg_all))) == nrow(deg_all), "\n")

# Diferentes umbrales
deg_fc05 <- subset(deg_all, adj.P.Val < 0.05 & abs(logFC) > 0.5)
deg_fc1 <- subset(deg_all, adj.P.Val < 0.05 & abs(logFC) > 1)
deg_fc2 <- subset(deg_all, adj.P.Val < 0.05 & abs(logFC) > 2)
deg_fc3 <- subset(deg_all, adj.P.Val < 0.05 & abs(logFC) > 3)

cat("\n--- eBayes (para construcción de GCN) ---\n")
cat("  Total genes significativos (FDR<0.05):", 
    sum(deg_all$adj.P.Val < 0.05), "\n")
cat("  DEGs |logFC|>0.5:", nrow(deg_fc05), "\n")
cat("  DEGs |logFC|>1.0:", nrow(deg_fc1), "\n")
cat("  DEGs |logFC|>2.0:", nrow(deg_fc2), "\n")
cat("  DEGs |logFC|>3.0:", nrow(deg_fc3), "\n")

# 6.2) Método TREAT (alta confianza)
fit_treat1 <- treat(fit2, lfc = 1)
deg_treat1 <- topTreat(fit_treat1, number = Inf, coef = "Metastatic_vs_Primary", 
                       sort.by = "none")
rownames(deg_treat1) <- rownames(expr_clean)
deg_treat1 <- deg_treat1[order(deg_treat1$adj.P.Val), ]
deg_treat1_sig <- subset(deg_treat1, adj.P.Val < 0.05)

fit_treat2 <- treat(fit2, lfc = 2)
deg_treat2 <- topTreat(fit_treat2, number = Inf, coef = "Metastatic_vs_Primary",
                       sort.by = "none")
rownames(deg_treat2) <- rownames(expr_clean)
deg_treat2 <- deg_treat2[order(deg_treat2$adj.P.Val), ]
deg_treat2_sig <- subset(deg_treat2, adj.P.Val < 0.05)

cat("\n--- TREAT (alta confianza estadística) ---\n")
cat("  Genes con |logFC|≥1 (FDR<0.05):", nrow(deg_treat1_sig), "\n")
cat("  Genes con |logFC|≥2 (FDR<0.05):", nrow(deg_treat2_sig), "\n")

# 6.3) Guardar resultados
write.csv(deg_all, "./results/DEGs/01_eBayes_all_genes.csv", row.names = TRUE)
write.csv(deg_fc05, "./results/DEGs/02_eBayes_FC05_for_GCN.csv", row.names = TRUE)
write.csv(deg_fc1, "./results/DEGs/03_eBayes_FC1.csv", row.names = TRUE)
write.csv(deg_fc2, "./results/DEGs/04_eBayes_FC2.csv", row.names = TRUE)
write.csv(deg_treat1, "./results/DEGs/05_TREAT_FC1_all.csv", row.names = TRUE)
write.csv(deg_treat1_sig, "./results/DEGs/06_TREAT_FC1_significant.csv", row.names = TRUE)
write.csv(deg_treat2, "./results/DEGs/07_TREAT_FC2_all.csv", row.names = TRUE)
write.csv(deg_treat2_sig, "./results/DEGs/08_TREAT_FC2_significant.csv", row.names = TRUE)

cat("\n  ✓ Tablas de DEGs guardadas en ./results/DEGs/\n")

# ================================================================================
# 7) VISUALIZACIONES
# ================================================================================

cat("\n═══ PASO 7: VISUALIZACIONES ═══\n")

# 7.1) Volcano plot
deg_all$regulation <- "No significant"
deg_all$regulation[deg_all$logFC >= 1 & deg_all$adj.P.Val < 0.05] <- "Up"
deg_all$regulation[deg_all$logFC <= -1 & deg_all$adj.P.Val < 0.05] <- "Down"
deg_all$regulation <- factor(deg_all$regulation, 
                             levels = c("Down", "No significant", "Up"))

top_up <- head(deg_all[deg_all$logFC > 0, ], 10)
top_down <- head(deg_all[deg_all$logFC < 0, ], 10)
top_genes <- rbind(top_up, top_down)

p_volcano <- ggplot(deg_all, aes(x = logFC, y = -log10(adj.P.Val), 
                                 color = regulation)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_text_repel(data = top_genes,
                  aes(label = rownames(top_genes)),
                  size = 2.5, max.overlaps = 20) +
  scale_color_manual(values = c("Down" = "#3498DB", 
                                "No significant" = "grey70", 
                                "Up" = "#E74C3C")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", 
             color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "black", linewidth = 0.5) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Volcano Plot: Metastatic vs Primary",
    subtitle = sprintf("Up: %d | Down: %d (FDR<0.05, |logFC|>1)",
                       sum(deg_all$regulation == "Up"),
                       sum(deg_all$regulation == "Down")),
    x = "log2 Fold Change",
    y = "-log10(FDR)",
    color = "Regulation"
  )

ggsave("./results/DEGs/09_volcano_plot.png", p_volcano, 
       width = 10, height = 8, dpi = 300)
cat("  ✓ Volcano plot guardado\n")

# 7.2) MA plot
pdf("./results/DEGs/10_MA_plot.pdf", width = 10, height = 6)
plotMA(fit_ebayes, main = "MA Plot: Metastatic vs Primary",
       status = deg_all$regulation[match(rownames(fit_ebayes), rownames(deg_all))],
       col = c("#3498DB", "grey70", "#E74C3C"))
abline(h = c(-1, 1), col = "black", lty = 2)
dev.off()
cat("  ✓ MA plot guardado\n")

# 7.3) Heatmap de top DEGs
if (nrow(deg_fc1) >= 20) {
  # Diagnóstico
  cat("\n  [DEBUG] Verificando compatibilidad de genes:\n")
  cat("    - Genes en deg_fc1:", nrow(deg_fc1), "\n")
  cat("    - Genes en expr_clean:", nrow(expr_clean), "\n")
  cat("    - Primeros 5 genes deg_fc1:", paste(head(rownames(deg_fc1), 5), collapse = ", "), "\n")
  cat("    - Primeros 5 genes expr_clean:", paste(head(rownames(expr_clean), 5), collapse = ", "), "\n")
  
  missing_genes <- setdiff(rownames(deg_fc1), rownames(expr_clean))
  cat("    - Genes faltantes:", length(missing_genes), "\n")
  
  # Solución: Usar genes que existan en ambos
  common_genes <- intersect(rownames(deg_fc1), rownames(expr_clean))
  cat("    - Genes comunes:", length(common_genes), "\n")
  
  if (length(common_genes) >= 10) {
    # Ordenar por logFC absoluto usando solo genes comunes
    deg_fc1_available <- deg_fc1[common_genes, ]
    deg_fc1_sorted <- deg_fc1_available[order(abs(deg_fc1_available$logFC), decreasing = TRUE), ]
    top_degs <- head(rownames(deg_fc1_sorted), min(50, nrow(deg_fc1_sorted)))
    
    ann_col_hm <- data.frame(
      Condition = group,
      row.names = colnames(expr_clean)
    )
    
    pheatmap(expr_clean[top_degs, ],
             scale = "row",
             show_rownames = TRUE,
             annotation_col = ann_col_hm,
             annotation_colors = list(Condition = colors_group),
             main = sprintf("Top %d DEGs (|logFC|>1)", length(top_degs)),
             fontsize_row = 6,
             filename = "./results/DEGs/11_heatmap_top50.png",
             width = 10, height = 12)
    cat("  ✓ Heatmap de top", length(top_degs), "DEGs guardado\n")
  } else {
    cat("  ⚠ ADVERTENCIA: Solo", length(common_genes), "genes comunes entre deg_fc1 y expr_clean\n")
    cat("  ⚠ Esto sugiere un problema en el pipeline. Revisar filtrado.\n")
  }
}

# 7.4) Distribución de logFC
p_fc <- ggplot(deg_all[deg_all$adj.P.Val < 0.05, ],
               aes(x = logFC, fill = logFC > 0)) +
  geom_histogram(bins = 50, alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"),
                    labels = c("Up", "Down")) +
  geom_vline(xintercept = c(-1, 0, 1), 
             linetype = c("dashed", "solid", "dashed"),
             color = "black", linewidth = 0.5) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Distribución de log2FC en DEGs significativos",
    x = "log2 Fold Change",
    y = "Número de genes",
    fill = "Direction"
  )

ggsave("./results/DEGs/12_logFC_distribution.png", p_fc, 
       width = 9, height = 6, dpi = 300)
cat("  ✓ Distribución de logFC guardada\n")

# ================================================================================
# 8) ANÁLISIS DE SUPERPOSICIÓN HVGs-DEGs
# ================================================================================

cat("\n═══ PASO 8: SUPERPOSICIÓN HVGs vs DEGs ═══\n")

# Análisis de overlap
deg_genes_fc05 <- rownames(deg_fc05)
deg_genes_fc1 <- rownames(deg_fc1)
deg_genes_fc2 <- rownames(deg_fc2)

overlap_fc05 <- intersect(hvgs, deg_genes_fc05)
overlap_fc1 <- intersect(hvgs, deg_genes_fc1)
overlap_fc2 <- intersect(hvgs, deg_genes_fc2)

overlap_summary <- data.frame(
  DEG_threshold = c("FC>0.5", "FC>1.0", "FC>2.0"),
  Total_DEGs = c(length(deg_genes_fc05), 
                 length(deg_genes_fc1), 
                 length(deg_genes_fc2)),
  DEGs_in_HVGs = c(length(overlap_fc05), 
                   length(overlap_fc1), 
                   length(overlap_fc2)),
  Percentage = c(
    round(100 * length(overlap_fc05) / length(deg_genes_fc05), 1),
    round(100 * length(overlap_fc1) / length(deg_genes_fc1), 1),
    round(100 * length(overlap_fc2) / length(deg_genes_fc2), 1)
  )
)

write.csv(overlap_summary, 
          "./results/QC/HVG_DEG_overlap_summary.csv", 
          row.names = FALSE)

cat("\n  Superposición HVGs (n=", length(hvgs), ") vs DEGs:\n")
print(overlap_summary)

# ================================================================================
# 9) RESUMEN Y RECOMENDACIONES
# ================================================================================

cat("\n\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║                     RESUMEN DEL ANÁLISIS                       ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("DATOS PROCESADOS:\n")
cat("  • Genes analizados:", nrow(expr_clean), "\n")
cat("  • Muestras Primary:", sum(group == "Primary"), "\n")
cat("  • Muestras Metastatic:", sum(group == "Metastatic"), "\n")
cat("  • HVGs identificados:", length(hvgs), "\n\n")

cat("RESULTADOS DEGs (eBayes):\n")
cat("  • Total significativos (FDR<0.05):", sum(deg_all$adj.P.Val < 0.05), "\n")
cat("  • |logFC|>0.5:", nrow(deg_fc05), 
    sprintf(" (%.1f%% en HVGs)\n", overlap_summary$Percentage[1]))
cat("  • |logFC|>1.0:", nrow(deg_fc1), 
    sprintf(" (%.1f%% en HVGs)\n", overlap_summary$Percentage[2]))
cat("  • |logFC|>2.0:", nrow(deg_fc2), 
    sprintf(" (%.1f%% en HVGs)\n", overlap_summary$Percentage[3]))

cat("\nRESULTADOS TREAT:\n")
cat("  • Evidencia robusta |logFC|≥1:", nrow(deg_treat1_sig), "\n")
cat("  • Evidencia robusta |logFC|≥2:", nrow(deg_treat2_sig), "\n")

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║            RECOMENDACIONES PARA ANÁLISIS GCN                   ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("✓ PARA CONSTRUCCIÓN DE RED (GCN):\n")
cat("  Archivo: ./results/QC/HVGs_top5000.csv\n")
cat("  Genes:", length(hvgs), "\n")
cat("  Razón: Mayor poder estadístico y detección de módulos\n\n")

cat("✓ PARA ENRIQUECIMIENTO DE MÓDULOS:\n")
cat("  Archivo primario: ./results/DEGs/02_eBayes_FC05_for_GCN.csv\n")
cat("  Genes:", nrow(deg_fc05), "\n")
cat("  Uso: Verificar si módulos están enriquecidos en DEGs\n\n")

cat("✓ PARA VALIDACIÓN EXPERIMENTAL:\n")
cat("  Archivo: ./results/DEGs/06_TREAT_FC1_significant.csv\n")
cat("  Genes:", nrow(deg_treat1_sig), "\n")
cat("  Razón: Alta confianza estadística\n\n")

cat("NOTA IMPORTANTE:\n")
cat("  ✓ Transformación log2 aplicada correctamente\n")
cat("  ✓ Datos listos para WGCNA u otros métodos de GCN\n")
cat("  ✓", round(overlap_summary$Percentage[2], 1), 
    "% de DEGs (FC>1) están en los HVGs\n\n")

# Guardar información de sesión
writeLines(capture.output(sessionInfo()), "./results/sessionInfo.txt")

cat("═══════════════════════════════════════════════════════════════\n")
cat("ANÁLISIS COMPLETADO:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Archivos generados:\n")
cat("  📁 ./results/QC/       → 6 figuras + tablas de control\n")
cat("  📁 ./results/DEGs/     → 8 tablas CSV + 4 visualizaciones\n")
cat("  📄 sessionInfo.txt     → Información de paquetes\n\n")




# ================================================================================
# ================================================================================
# POSTERIOR A LAS REDES Y CLUSTERING - RESUMEN Y RECOMENDACIONES
# ================================================================================
# ================================================================================

# Leer genes del cluster (una columna, sin cabecera)
clust2_genes <- read_tsv("Cluster_2.txt", col_names = FALSE)[[1]]
clust2_df <- data.frame(gene = clust2_genes, stringsAsFactors = FALSE)


# 2) Preparar tabla de DEGs
deg_tbl <- deg_treat1_sig %>%
  as.data.frame() %>%
  rownames_to_column("gene")   # convierte rownames en columna "gene"


# 3) Cruzar cluster con DEGs y resumir
clust2_annot <- clust2_df %>%
  left_join(deg_tbl, by = "gene")

# Resumen: cuántos genes y cuántos son DEGs
clust2_summary <- clust2_annot %>%
  summarise(
    n_genes    = n(),
    n_DEGs     = sum(!is.na(adj.P.Val) & adj.P.Val < 0.05),
    prop_DEGs  = round(100 * n_DEGs / n_genes, 1),
    mean_logFC = mean(logFC[!is.na(logFC)], na.rm = TRUE)
  )

clust2_summary


