# ================================================================================
# Análisis de Enriquecimiento Funcional (GO/KEGG)
# ================================================================================

# 1. Carga de librerías necesarias
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(hgu133plus2.db)
  library(enrichplot)
  library(ggplot2)
})

# 2. Configuración de rutas (Ruta proporcionada por el usuario)
base_path   <- "C:/Users/marcr/OneDrive/Escritorio/Master_Bioinf/Semestre_4/TFM/Resultados"
input_dir   <- file.path(base_path, "inputClusters")
output_dir  <- file.path(base_path, "Resultados_Enriquecimiento")

# Crear carpeta de resultados si no existe
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 3. Identificar todos los archivos Cluster_X.txt en la carpeta
# Esto evita tener que escribir los nombres uno por uno
cluster_files <- list.files(path = input_dir, pattern = "Cluster_.*\\.txt", full.names = TRUE)

cat("Se han encontrado", length(cluster_files), "clústeres para analizar.\n")

# 4. Bucle para procesar cada clúster automáticamente
for (file_path in cluster_files) {
  
  # Extraer el nombre del clúster para los archivos de salida
  cluster_name <- gsub("\\.txt", "", basename(file_path))
  cat("\n>>> Procesando:", cluster_name, "...\n")
  
  # Leer la lista de genes (Símbolos)
  gene_symbols <- readLines(file_path)
  
  # A. Mapeo de Gene Symbols a Entrez IDs
  # Es necesario para clusterProfiler
  entrez_ids <- mapIds(hgu133plus2.db, 
                       keys      = gene_symbols, 
                       column    = 'ENTREZID', 
                       keytype   = 'SYMBOL', 
                       multiVals = "first")
  
  # Limpiar IDs no mapeados
  entrez_ids_clean <- entrez_ids[!is.na(entrez_ids)]
  
  if (length(entrez_ids_clean) == 0) {
    cat("    ⚠ No se encontraron Entrez IDs para", cluster_name, ". Saltando...\n")
    next
  }
  
  # B. Análisis de Enriquecimiento GO (Gene Ontology)
  # Aplicando el método de corrección BH indicado en tu script original
  go_res <- enrichGO(gene          = entrez_ids_clean,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "all",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE)
  
  if (!is.null(go_res) && nrow(go_res@result) > 0) {
    write.csv(go_res, file.path(output_dir, paste0(cluster_name, "_GO_results.csv")), row.names = FALSE)
    
    # Generar Dotplot de GO
    p_go <- dotplot(go_res, split="ONTOLOGY", font.size = 10) + 
      facet_grid(ONTOLOGY~., scale ='free', space = 'free')
    ggsave(filename = file.path(output_dir, paste0(cluster_name, "_GO_plot.pdf")), 
           plot = p_go, device = "pdf", height = 10, width = 8)
  }
  
  # C. Análisis de Enriquecimiento KEGG (Pathways)
  kegg_res <- enrichKEGG(gene          = entrez_ids_clean,
                         organism      = 'hsa',
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05)
  
  if (!is.null(kegg_res) && any(kegg_res@result$p.adjust < 0.05)) {
    write.csv(kegg_res, file.path(output_dir, paste0(cluster_name, "_KEGG_results.csv")), row.names = FALSE)
    
    # Generar Dotplot de KEGG
    p_kegg <- dotplot(kegg_res) + ggtitle(paste("KEGG Enrichment -", cluster_name))
    ggsave(filename = file.path(output_dir, paste0(cluster_name, "_KEGG_plot.pdf")), 
           plot = p_kegg, device = "pdf", height = 7, width = 8)
  }
}

cat("\n✅ Proceso finalizado. Los resultados están en:", output_dir, "\n")