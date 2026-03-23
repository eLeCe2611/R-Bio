## 1. Definir la ruta completa al archivo
file_path <- "C:/Users/luisc/Desktop/Nueva carpeta/GSE227018_series_matrix.txt"

# Leer todas las líneas del archivo
lines <- readLines(file_path)

# Identificar las líneas donde comienza y termina la tabla de expresión
start_line <- grep("!series_matrix_table_begin", lines)
end_line   <- grep("!series_matrix_table_end", lines)

# Extraer únicamente el bloque que contiene la matriz de expresión
data_text <- lines[(start_line + 1):(end_line - 1)]

# Convertir el bloque de texto en una tabla de datos
exp_data <- read.table(text = data_text, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# Visualizar las primeras filas para confirmar que se ha leido correctamente
head(exp_data)


## 2. Preparar la matriz numérica de expresión y calcular correlaciones

# Separamos los identificadores y la matriz de expresión (convirtiéndola a numérica)
genes <- exp_data[, 1]  # Identificadores de genes
expr_matrix <- as.matrix(exp_data[, -1])
mode(expr_matrix) <- "numeric"  # Forzamos a numérico, en caso de que sea necesario

# Para calcular la correlación entre genes (cada gen representado por su perfil de expresión en las muestras),
# se debe transponer la matriz: filas=genes, columnas=muestras -> t(expr_matrix) da filas = muestras.
# Así cor(t(expr_matrix)) genera la correlación entre genes.
cor_pearson <- cor(t(expr_matrix), method = "pearson")
cor_spearman <- cor(t(expr_matrix), method = "spearman")

## 3. Umbralización de las correlaciones y construcción de matrices de adyacencia

# Definir un umbral de interés. En este ejemplo, usamos 0.8.
threshold <- 0.8

# Convertir las matrices de correlación en matrices de adyacencia.
# Se toman valores absolutos para considerar tanto correlaciones positivas como negativas.
adj_pearson <- (abs(cor_pearson) >= threshold) * 1
diag(adj_pearson) <- 0   # Eliminar autorrelaciones (valor 1 en la diagonal)

adj_spearman <- (abs(cor_spearman) >= threshold) * 1
diag(adj_spearman) <- 0

## 4. Creación de las redes de co-expresión utilizando igraph

# Instalar y cargar el paquete si aún no lo tienes:
# install.packages("igraph")
library(igraph)

# Generar grafos a partir de las matrices de adyacencia.
# El parámetro mode="undirected" es adecuado ya que la correlación es simétrica.
g_pearson <- graph_from_adjacency_matrix(adj_pearson, mode = "undirected", weighted = TRUE, diag = FALSE)
g_spearman <- graph_from_adjacency_matrix(adj_spearman, mode = "undirected", weighted = TRUE, diag = FALSE)

# Asignar los nombres de los nodos usando los identificadores de genes
V(g_pearson)$name <- genes
V(g_spearman)$name <- genes


## 5. Análisis y comparación de las redes

# Resumen básico de los grafos
cat("Red basada en Pearson:\n")
print(summary(g_pearson))
cat("\nRed basada en Spearman:\n")
print(summary(g_spearman))

# Número de nodos y enlaces:
cat("\nNodos y enlaces (Pearson):", vcount(g_pearson), "nodos,", ecount(g_pearson), "enlaces\n")
cat("Nodos y enlaces (Spearman):", vcount(g_spearman), "nodos,", ecount(g_spearman), "enlaces\n")

# Distribución de grados
deg_pearson <- degree(g_pearson)
deg_spearman <- degree(g_spearman)

cat("\nGrado promedio (Pearson):", mean(deg_pearson), "\n")
cat("Grado promedio (Spearman):", mean(deg_spearman), "\n")

# Coeficiente de clustering promedio (indicador de cohesión local)
clust_pearson <- transitivity(g_pearson, type = "average")
clust_spearman <- transitivity(g_spearman, type = "average")

cat("\nCoeficiente de clustering (Pearson):", clust_pearson, "\n")
cat("Coeficiente de clustering (Spearman):", clust_spearman, "\n")

# Opcional: visualizar las redes (para datasets pequeños)
# plot(g_pearson, main = "Red de Co-expresión (Pearson)", vertex.size = 5, vertex.label.cex = 0.7)
plot(g_spearman, main = "Red de Co-expresión (Spearman)", vertex.size = 5, vertex.label.cex = 0.7)

# Exportar la red basada en Pearson en formato GraphML al directorio "C:/Users/luisc/Desktop/Cosas"
write_graph(g_pearson, file = "C:/Users/luisc/Desktop/Cosas/Red_pearson.sif", format = "graphml")

