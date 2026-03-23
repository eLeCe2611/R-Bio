# Cargar librerías necesarias
library(data.table)  # Para carga eficiente de datos grandes
library(cluster)     # Para clustering
library(factoextra)  # Para visualización

# Definir rutas de archivo
input_file <- "C:/Users/luisc/Desktop/Cosas/Universidad/3º Tercero/Segundo Cuatrismestre/Bioinformática/Expander/MicroArray/GSE27272_non-normalized_cord_blood.txt"
output_file <- "C:/Users/luisc/Desktop/MicroArray_Procesado.txt"

# Cargar los datos
df <- fread(input_file, header = TRUE, sep = "\t", na.strings = c("", "NA"))

# Revisar estructura
print(head(df))

# ------- LIMPIEZA DE DATOS --------

# Eliminar la primera columna si solo contiene IDs o valores vacíos
if (all(is.na(df[[1]]))) {
  df <- df[, -1, with = FALSE]
}

# Asegurar que todas las columnas sean numéricas
numeric_cols <- names(df)[sapply(df, is.character)]
df[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]

# Eliminar columnas con más del 30% de valores NA
threshold <- 0.3 * nrow(df)
df <- df[, which(colSums(is.na(df)) < threshold), with = FALSE]

# Eliminar filas con valores NA restantes
df <- na.omit(df)

# -------- NORMALIZACIÓN --------

# Normalizar los datos (Z-score)
df_scaled <- scale(df)

# -------- CLUSTERING (si es necesario) --------

# Determinar número óptimo de clusters con el método del codo
fviz_nbclust(df_scaled, kmeans, method = "wss")

# Aplicar k-means clustering con el número óptimo de clusters (ejemplo: k=3)
set.seed(123)
kmeans_result <- kmeans(df_scaled, centers = 3, nstart = 25)

# Agregar el cluster al dataset
df$Cluster <- factor(kmeans_result$cluster)

# Guardar los datos procesados en la ruta especificada
write.table(df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

# Mensaje de éxito
cat("✅ Archivo guardado en:", output_file, "\n")
