# 1. Instalación y carga de paquetes necesarios
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Instalar los paquetes necesarios: GEOquery, org.Hs.eg.db, topGO
BiocManager::install("GEOquery", ask = FALSE)
BiocManager::install("org.Hs.eg.db", ask = FALSE)
BiocManager::install("topGO", ask = FALSE)
BiocManager::install("Rgraphviz", ask = FALSE)  # Opcional para visualización

library(GEOquery)
library(org.Hs.eg.db)
library(topGO)
library(Rgraphviz)  # Opcional

# 2. Descarga y exploración del dataset GSE227018
gse <- getGEO("GSE227018", GSEMatrix = TRUE)
cat("Número de ExpressionSet descargados:", length(gse), "\n")

# Selecciona el primer ExpressionSet (en este ejemplo, suponemos que ese es el que corresponde)
eset <- gse[[1]]

# Extrae la matriz de expresión y la información fenotípica
expression_matrix <- exprs(eset)
sample_info <- pData(eset)
cat("Primeras filas de la matriz de expresión:\n")
head(expression_matrix)
cat("\nInformación de las muestras:\n")
head(sample_info)

# 3. Preparación del vector de genes para el análisis de enriquecimiento
# Notar que en este dataset la columna "ID_REF" (los nombres de fila) contiene distintos formatos;
# aprovechamos aquellos que comienzan con "ENSG" (IDs de Ensembl) para poder usar el paquete org.Hs.eg.db.
probeIDs <- rownames(expression_matrix)
ensgIDs <- probeIDs[grep("^ENSG", probeIDs)]
# Muchas veces los IDs vienen con sufijos (ej. "_st" o "_s_st"), así que extraemos la parte principal:
ensgIDs_clean <- sapply(ensgIDs, function(x) unlist(strsplit(x, "_"))[1])
ensgIDs_clean <- unique(ensgIDs_clean)  # Nos quedamos con IDs únicos

# Para demostrar, simulamos un vector de “p‑valores” (en tu análisis real estos provendrían de, por ejemplo, limma)
set.seed(123)
geneList <- runif(length(ensgIDs_clean), min = 0, max = 1)
names(geneList) <- ensgIDs_clean

# Definimos la función que selecciona los genes “diferenciales”
topDiffGenes <- function(pval) {
  return(pval < 0.01)
}

# 4. Creación del objeto topGOdata para la ontología "CC" (Componente Celular)
# Usamos annFUN.org con mapping = "org.Hs.eg.db" y especificamos que los IDs son de tipo "ensembl".
GOdata_CC <- new("topGOdata",
                 description = "Enriquecimiento en Componente Celular - GSE227018",
                 ontology = "CC", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 nodeSize = 10,
                 annot = annFUN.org,
                 mapping = "org.Hs.eg.db",
                 ID = "ensembl")

# 5. Ejecución del test de enriquecimiento (test de Fisher clásico)
resultFisher_CC <- runTest(GOdata_CC, algorithm = "classic", statistic = "fisher")

# 6. Generación de una tabla con los términos GO más significativos en CC
allRes_CC <- GenTable(GOdata_CC,
                      classicFisher = resultFisher_CC,
                      orderBy = "classicFisher",
                      topNodes = 10)
cat("\nTabla de resultados (top 10 términos) para Componente Celular:\n")
print(allRes_CC)

# 7. (Opcional) Visualización del árbol ontológico
cat("\nVisualización del grafo GO para los términos significativos:\n")
showSigOfNodes(GOdata_CC, score(resultFisher_CC), firstSigNodes = 5, useInfo = "all")
