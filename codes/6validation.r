# ================================================================================
# Validación de Clústeres: Integración de DEGs (TREAT) y Hubs
# ================================================================================

library(dplyr)
library(readr)
library(stringr)

# 1. Rutas de archivos
path_base     <- "C:/Users/marcr/OneDrive/Escritorio/Master_Bioinf/Semestre_4/TFM/Resultados"
path_clusters <- file.path(path_base, "inputClusters")
path_degs     <- file.path(path_base, "./DEGs/05_TREAT_FC1_all.csv")
path_hubs     <- file.path(path_base, "infoHubs.txt")
output_dir    <- file.path(path_base, "Validacion_Final_Clusters")

if (!dir.exists(output_dir)) dir.create(output_dir)

# 2. Cargar y procesar la lista de Hubs
# El archivo contiene una cadena de texto con los nombres de los genes
hubs_raw <- readLines(path_hubs)
# Extraer solo los nombres de los genes (limpiando corchetes, comillas y espacios)
hubs_list <- str_extract_all(hubs_raw, "[A-Z0-9a-z\\.\\-]+")[[1]]
hubs_list <- hubs_list[!(hubs_list %in% c("Nodos", "hub", "c"))] # Limpiar etiquetas

# 3. Cargar tabla de DEGs
degs_all <- read_csv(path_degs) %>%
  rename(gene = 1) %>%
  mutate(is_DEG = adj.P.Val < 0.05)

# 4. Procesar Clústeres: 72, 54, 26, 42
cluster_ids <- c(72, 54, 26, 42)
resumen_biomarcadores <- data.frame()

for (id in cluster_ids) {
  file_name <- file.path(path_clusters, paste0("Cluster_", id, ".txt"))
  
  if (file.exists(file_name)) {
    cluster_genes <- readLines(file_name)
    
    # Crear tabla de validación
    valid_df <- data.frame(gene = cluster_genes) %>%
      left_join(degs_all, by = "gene") %>%
      mutate(is_Hub = gene %in% hubs_list) %>%
      mutate(Priority_Score = case_when(
        is_DEG & is_Hub  ~ "ALTA (Hub + DEG)",
        is_Hub           ~ "MEDIA (Solo Hub)",
        is_DEG           ~ "MEDIA (Solo DEG)",
        TRUE             ~ "BAJA"
      ))
    
    # Guardar reporte detallado del clúster
    write.csv(valid_df, file.path(output_dir, paste0("Validacion_Cluster_", id, ".csv")), row.names = FALSE)
    
    # Añadir al resumen global
    resumen_biomarcadores <- rbind(resumen_biomarcadores, data.frame(
      Cluster = id,
      Total_Genes = nrow(valid_df),
      Num_Hubs = sum(valid_df$is_Hub),
      Num_DEGs = sum(valid_df$is_DEG, na.rm = TRUE),
      Candidatos_Alta_Prioridad = sum(valid_df$Priority_Score == "ALTA (Hub + DEG)")
    ))
  }
}

print(resumen_biomarcadores)
write.csv(resumen_biomarcadores, file.path(output_dir, "Resumen_Priorizacion_Biomarcadores.csv"), row.names = FALSE)
