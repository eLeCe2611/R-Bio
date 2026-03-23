# ===========================================
# INSTALACIÓN Y CARGA DE PAQUETES NECESARIOS
# ===========================================

# Instalar paquetes de CRAN
install.packages(c("cluster", "factoextra", "biclust", "caret", "NbClust"))

# Instalar Bioconductor si no está presente
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar paquetes desde Bioconductor
BiocManager::install(c("STRINGdb", "Biobase", "GEOquery"), ask = FALSE)

# Cargar librerías necesarias para el análisis
library(STRINGdb)
library(igraph)
library(Biobase)
library(caret)
library(cluster)
library(factoextra)
library(pROC)
library(GEOquery)

# ===========================================
# PASO 1: DESCARGA Y PREPROCESADO DE DATOS
# ===========================================

# Descargar el dataset GDS6063 desde GEO
gds_obj <- getGEO("GDS6063")

# Extraer la tabla de expresión
expr_data <- Table(gds_obj)

# Verificar estructura de los datos
head(expr_data)

# Seleccionar solo las columnas con datos de expresión (GSM)
expr_matrix <- expr_data[, grep("^GSM", colnames(expr_data))]

# Convertir a matriz numérica
expr_matrix <- as.matrix(sapply(expr_matrix, as.numeric))

# Asignar nombres de fila con identificadores de genes
rownames(expr_matrix) <- expr_data$IDENTIFIER

# Eliminar genes con valores faltantes
expr_matrix <- na.omit(expr_matrix)

# Reducir a los primeros 300 genes para evitar sobrecarga computacional
expr_matrix <- expr_matrix[1:300, ]

# Confirmar estructura final
head(expr_matrix)

# ===========================================
# PASO 2: CONSTRUCCIÓN DE RED DE CORRELACIÓN
# ===========================================

# Calcular matriz de correlación entre genes (Pearson)
gene_cor <- cor(t(expr_matrix), method = "pearson")

# Crear grafo a partir de la matriz de correlación
g_net <- graph.adjacency(gene_cor, weighted = TRUE, mode = "lower")
V(g_net)$name <- rownames(expr_matrix)

# Eliminar enlaces con correlación menor a 0.80
g_net <- delete.edges(g_net, E(g_net)[weight < 0.80])

# Simplificar la red: eliminar bucles y enlaces múltiples
g_net <- simplify(g_net, remove.multiple = TRUE, remove.loops = TRUE)

# Eliminar nodos sin conexiones
g_net <- delete.vertices(g_net, which(degree(g_net) < 1))

# Visualizar la red
plot(g_net, layout = layout_with_fr, vertex.size = 5, vertex.label = NA)

# ===========================================
# PASO 3: VALIDACIÓN BIOLÓGICA CON STRINGdb
# ===========================================

# Preparar identificadores de genes para STRING
genes_ids <- data.frame(gene = V(g_net)$name)

# Inicializar el objeto STRINGdb para Mus musculus (TaxID 10090)
string_db <- STRINGdb$new(
  version = "11",
  species = 10090,
  score_threshold = 0,
  input_directory = ""
)

# Mapear los identificadores a STRING
mapped_ids <- string_db$map(genes_ids, "gene", removeUnmappedRows = TRUE)

# Extraer identificadores mapeados correctamente
hits <- mapped_ids$STRING_id

# Visualizar la red de interacción conocida según STRING
string_db$plot_network(hits)

# ===========================================
# PASO 4: ENRIQUECIMIENTO FUNCIONAL
# ===========================================

# Obtener enriquecimiento en procesos biológicos (GO)
enrichment_go <- string_db$get_enrichment(hits, category = "Process")

# Obtener enriquecimiento en rutas KEGG
enrichment_kegg <- string_db$get_enrichment(hits, category = "KEGG")

# Mostrar resultados principales
head(enrichment_go)
head(enrichment_kegg)

# ===========================================
# PASO 5: ANÁLISIS DE CLÚSTERES FUNCIONALES
# ===========================================

# Obtener clústeres funcionales en la red
clusters_list <- string_db$get_clusters(hits)

# Visualizar los 4 primeros clústeres
par(mfrow = c(2, 2))
for (i in 1:4) {
  string_db$plot_network(clusters_list[[i]])
}
