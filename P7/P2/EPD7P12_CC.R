# Instalar (si es necesario) y cargar los paquetes de Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar clusterProfiler y la anotación para organismo humano
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"), ask = FALSE)

# Cargar las librerías
install.packages("ggplot2")  # Ejecuta esto si ggplot2 no está instalado
library(ggplot2)
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


ego_CC <- enrichGO(gene          = genes_unique,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "CC",       # Analizamos componente celular
                   pAdjustMethod = "BH",       # Corrección de múltiples pruebas
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)       # Convierte a símbolos de gen para facilitar la interpretación

# Mostrar los primeros resultados
head(as.data.frame(ego_CC))


# Crear un dotplot para visualizar los términos enriquecidos en CC
dotplot(ego_CC, showCategory = 10) + ggtitle("Enriquecimiento en Componente Celular")
