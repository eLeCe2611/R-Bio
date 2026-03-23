# Cargar librerías necesarias
library(caret)       # Preprocesamiento de datos
library(cluster)     # Algoritmos de clusterización
library(factoextra)  # Visualización de clustering
library(data.table)  # Manejo eficiente de datos

# Definir la ruta del archivo
ruta_archivo <- "C:/Users/luisc/Desktop/Cosas/Universidad/3º Tercero/Segundo Cuatrismestre/Bioinformática/Expander/RNA-Seq/GSE100797_ProcessedData.txt"

# Cargar los datos
dataGDS <- fread(ruta_archivo, header = TRUE, sep = "\t")

# Ver las primeras filas para inspección
head(dataGDS)

# --- LIMPIEZA DE DATOS ---
# Eliminar valores nulos
dataGDS <- na.omit(dataGDS)

# Identificar y reemplazar valores atípicos por la media de la columna
samplesGDS <- colnames(dataGDS)[-1]  # Excluir primera columna si es identificador de genes

for (col in samplesGDS) {
  Q1 <- quantile(dataGDS[[col]], 0.25, na.rm = TRUE)
  Q3 <- quantile(dataGDS[[col]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  limitDown <- Q1 - 1.5 * IQR
  limitUp <- Q3 + 1.5 * IQR
  media <- mean(dataGDS[[col]], na.rm = TRUE)

  dataGDS[[col]][dataGDS[[col]] >= limitUp] <- media
  dataGDS[[col]][dataGDS[[col]] <= limitDown] <- media
}

# --- PREPROCESAMIENTO ---
# Preprocesamiento: centrar y escalar los datos
preprocessParams <- preProcess(dataGDS[, -1, with = FALSE], method = c("center", "scale"))

# Convertir a data.frame antes de la asignación para evitar problemas de compatibilidad
dataGDS_scaled <- as.data.frame(predict(preprocessParams, dataGDS[, -1, with = FALSE]))

# Reintegrar la columna de identificadores (asumiendo que la primera columna es identificador)
dataGDS <- cbind(dataGDS[, 1, with = FALSE], dataGDS_scaled)

# Verificar estructura final
str(dataGDS)

# --- CLUSTERIZACIÓN ---
# Determinar el número óptimo de clusters con el método del codo
wss <- numeric(15)
for (i in 1:15) {
  kmeans_result <- kmeans(dataGDS[, -1, with = FALSE], centers = i, nstart = 10)
  wss[i] <- kmeans_result$tot.withinss
}

# Graficar la selección del número óptimo de clusters
plot(1:15, wss, type = "b", xlab = "Número de Clusters", ylab = "Suma de cuadrados dentro del grupo",
     main = "Método del codo para selección de K")

# Aplicar K-Means con el número óptimo (cambiar este número según el gráfico)
numClusters <- 6  # Ajustar según el gráfico
set.seed(100)
clusters <- kmeans(dataGDS[, -1, with = FALSE], centers = numClusters, nstart = 10)

# Agregar los clusters al dataset
dataGDS$Cluster <- as.factor(clusters$cluster)

# Extraer solo las columnas numéricas (excluyendo identificadores y la columna de clusters)
data_numeric <- dataGDS[, -c(1, ncol(dataGDS)), with = FALSE]  

# Verificar que todas las columnas sean numéricas
str(data_numeric)

# Visualizar los clusters con PCA
fviz_cluster(clusters, data = data_numeric, geom = "point", 
             ellipse = TRUE, main = "Visualización de Clusters")

# Guardar el dataset procesado con clusters
fwrite(dataGDS, "C:/Users/luisc/Desktop/RNA_Procesado.txt", sep = "\t")

# Mostrar resumen de los clusters
table(dataGDS$Cluster)
