# Instalar (si es necesario) y cargar los paquetes de Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar clusterProfiler y la anotación para organismo humano
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"), ask = FALSE)

# Cargar las librerías
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("C:/Users/luisc/Desktop/Nueva carpeta/P2")

gene_data <- read.table("GSE261149_gene_count_matrix.txt",
                        header = FALSE, 
                        stringsAsFactors = FALSE)
genes_raw <- gene_data$V1

# Eliminar la parte de la versión y el conteo: quitamos todo desde el primer punto
genes_clean <- sapply(genes_raw, function(x) sub("\\..*", "", x))
# Mantener únicamente los IDs únicos
genes_unique <- unique(genes_clean)

# Análisis de enriquecimiento para Función Molecular (MF)
ego_MF <- enrichGO(gene          = genes_unique,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "MF",       # Cambiamos a Función Molecular
                   pAdjustMethod = "BH",       # Corrección de múltiples pruebas
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)       # Convierte a símbolos de gen para facilitar la interpretación

# Mostrar los primeros resultados de Función Molecular
head(as.data.frame(ego_MF))

# Crear un dotplot para visualizar los términos enriquecidos en Función Molecular
dotplot(ego_MF, showCategory = 10) + ggtitle("Enriquecimiento en Función Molecular")
