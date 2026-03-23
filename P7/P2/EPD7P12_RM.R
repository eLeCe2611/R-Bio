# Instalar (si es necesario) y cargar los paquetes requeridos
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar clusterProfiler y la anotación para organismo humano
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"), ask = FALSE)

# Cargar las librerías
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("C:/Users/luisc/Desktop/Nueva carpeta/P2")

# Leer el archivo con la lista de géneros
gene_data <- read.table("GSE261149_gene_count_matrix.txt",
                        header = FALSE, 
                        stringsAsFactors = FALSE)
genes_raw <- gene_data$V1

# Limpiar los IDs: eliminar la parte de la versión y, en algunos casos, el conteo
genes_clean <- sapply(genes_raw, function(x) sub("\\..*", "", x))
genes_unique <- unique(genes_clean)

# Convertir los IDs de ENSEMBL a ENTREZID para el análisis KEGG
gene_conversion <- bitr(genes_unique,
                        fromType = "ENSEMBL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

entrez_genes <- unique(gene_conversion$ENTREZID)

# Análisis de enriquecimiento con KEGG (ruta metabólica)
ekegg <- enrichKEGG(gene          = entrez_genes,
                    organism      = "hsa",         # hsa para humano
                    pAdjustMethod = "BH",          # Corrección de pruebas múltiples
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)

# Mostrar los primeros resultados de la enriquecimiento en KEGG
head(as.data.frame(ekegg))

# Crear un dotplot para visualizar las rutas KEGG enriquecidas
dotplot(ekegg, showCategory = 10) + ggtitle("Enriquecimiento en rutas KEGG")
