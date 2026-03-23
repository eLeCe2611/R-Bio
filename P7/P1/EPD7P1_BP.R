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

# Selecciona el primer ExpressionSet
eset <- gse[[1]]

# Extrae la matriz de expresión y la información fenotípica
expression_matrix <- exprs(eset)
sample_info <- pData(eset)
cat("Primeras filas de la matriz de expresión:\n")
print(head(expression_matrix))
cat("\nInformación de las muestras:\n")
print(head(sample_info))

# 3. Preparación del vector de genes para el análisis de enriquecimiento
# Extraemos los IDs de las filas que inician con "ENSG" (IDs de Ensembl):
probeIDs <- rownames(expression_matrix)
ensgIDs <- probeIDs[grep("^ENSG", probeIDs)]
# Muchos IDs pueden venir con sufijos (p.ej. "_st" o "_s_st"); extraemos la parte principal:
ensgIDs_clean <- sapply(ensgIDs, function(x) unlist(strsplit(x, "_"))[1])
ensgIDs_clean <- unique(ensgIDs_clean)

# Validamos que estos IDs existan en la base de datos org.Hs.eg.db:
keys_db <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
validIDs <- intersect(ensgIDs_clean, keys_db)
cat("Número de IDs de Ensembl válidos:", length(validIDs), "\n")

# Simulamos un vector de “p‑valores” para cada ID válido (en un estudio real estos provendrían de un análisis diferencial)
set.seed(123)
geneList <- runif(length(validIDs), min = 0, max = 1)
names(geneList) <- validIDs

# Definimos la función que selecciona los genes “diferenciales”
# Para la demostración usamos p < 0.5 para seleccionar una porción significativa de genes
topDiffGenes <- function(pval) {
  return(pval < 0.5)
}

# 4. Creación del objeto topGOdata para Proceso Biológico (BP)
# Se usa el argumento ontology = "BP" y se especifica que los IDs son de tipo "ENSEMBL"
GOdata_BP <- new("topGOdata",
                 description = "Enriquecimiento en Proceso Biológico - GSE227018",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 nodeSize = 10,
                 annot = annFUN.org,
                 mapping = "org.Hs.eg.db",
                 ID = "ENSEMBL")

# 5. Ejecución del test de enriquecimiento (test de Fisher clásico)
resultFisher_BP <- runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")

# 6. Generación de una tabla con los términos GO (top 10) más significativos en Proceso Biológico
allRes_BP <- GenTable(GOdata_BP,
                      classicFisher = resultFisher_BP,
                      orderBy = "classicFisher",
                      topNodes = 10)
cat("\nTabla de resultados (top 10 términos) para Proceso Biológico:\n")
print(allRes_BP)

# 7. (Opcional) Visualización del árbol ontológico para los términos significativos en BP
cat("\nVisualización del grafo GO para los términos significativos en BP:\n")
showSigOfNodes(GOdata_BP, score(resultFisher_BP), firstSigNodes = 5, useInfo = "all")
