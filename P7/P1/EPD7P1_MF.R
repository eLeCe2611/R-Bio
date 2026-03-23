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

# Seleccionamos el primer ExpressionSet
eset <- gse[[1]]

# Extraemos la matriz de expresión y la información fenotípica
expression_matrix <- exprs(eset)
sample_info <- pData(eset)
cat("Primeras filas de la matriz de expresión:\n")
head(expression_matrix)
cat("\nInformación de las muestras:\n")
head(sample_info)

# 3. Preparación del vector de genes para el análisis de enriquecimiento
# En este dataset los nombres de fila tienen distintos formatos;
# filtramos aquellos que comienzan con "ENSG" (IDs de Ensembl).
probeIDs <- rownames(expression_matrix)
ensgIDs <- probeIDs[grep("^ENSG", probeIDs)]
# Muchos IDs pueden venir con sufijos (ej. "_st" o "_s_st"), extraemos la parte principal:
ensgIDs_clean <- sapply(ensgIDs, function(x) unlist(strsplit(x, "_"))[1])
ensgIDs_clean <- unique(ensgIDs_clean)

# Filtramos solo aquellos IDs que se encuentran en la base de datos org.Hs.eg.db
keys_db <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
validIDs <- intersect(ensgIDs_clean, keys_db)
cat("Número de IDs de Ensembl válidos:", length(validIDs), "\n")

result <- head(select(org.Hs.eg.db,
                      keys = validIDs[1:5],
                      keytype = "ENSEMBL",
                      columns = c("GO", "ONTOLOGY")))
print(result)

# Para el ejemplo, simulamos un vector de “p‑valores” para cada ID válido.
# En un análisis real estos p‑valores provendrían de un análisis diferencial.
set.seed(123)
geneList <- runif(length(validIDs), min = 0, max = 1)
names(geneList) <- validIDs

# Para la demostración, definimos como “genes diferenciales” aquellos con p-valor < 0.5.
# (Recuerda: en tu estudio real usarías el umbral que considere adecuado.)
topDiffGenes <- function(pval) {
  return(pval < 0.5)
}

# 4. Creación del objeto topGOdata para Función Molecular (MF)
# Se usa la ontología "MF" y se especifica que los IDs son de tipo "ENSEMBL" (en mayúsculas).
GOdata_MF <- new("topGOdata",
                 description = "Enriquecimiento en Función Molecular - GSE227018",
                 ontology = "MF", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 nodeSize = 10,
                 annot = annFUN.org,
                 mapping = "org.Hs.eg.db",
                 ID = "ENSEMBL")

# 5. Ejecución del test de enriquecimiento (test de Fisher clásico)
resultFisher_MF <- runTest(GOdata_MF, algorithm = "classic", statistic = "fisher")

# 6. Generación de una tabla con los términos GO (top 10) más significativos para Función Molecular (MF)
allRes_MF <- GenTable(GOdata_MF,
                      classicFisher = resultFisher_MF,
                      orderBy = "classicFisher",
                      topNodes = 10)
cat("\nTabla de resultados (top 10 términos) para Función Molecular:\n")
print(allRes_MF)

# 7. (Opcional) Visualización del árbol ontológico
cat("\nVisualización del grafo GO para los términos significativos en MF:\n")
showSigOfNodes(GOdata_MF, score(resultFisher_MF), firstSigNodes = 5, useInfo = "all")
