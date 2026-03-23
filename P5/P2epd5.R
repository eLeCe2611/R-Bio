# Instalar y cargar librerías necesarias
install.packages(c("factoextra", "cluster", "pheatmap", "biclust"))
library(factoextra)
library(cluster)
library(pheatmap)
library(biclust)

# Leer los datos sin establecer row.names al principio
data <- read.table("C:/Users/luisc/Desktop/Cosas/Universidad/3º Tercero/Segundo Cuatrismestre/Bioinformática/Expander/RNA-Seq/GSE261149_gene_count_matrix.txt", header = TRUE, row.names=1)

# Eliminar genes sin expresión
data_filtered <- data[rowSums(data) > 0, ]

# Normalización log2(1+valor) para corregir escala
data_log <- log2(data_filtered + 1)

# Para reproducibilidad
set.seed(123)  

# Seleccionar una muestra aleatoria de 500 genes (ajusta según sea necesario)
sample_data <- data_log[sample(nrow(data_log), 500), ]

# Aplicar método "Gap Statistic" en la muestra
fviz_nbclust(sample_data, kmeans, method = "gap_stat")

# Aplicar k-means con k = 4 (como ejemplo)
set.seed(123)
clusters_km <- kmeans(sample_data, centers=4, nstart=25)

# Visualizar clustering
fviz_cluster(clusters_km, data = sample_data)

# Calcular matriz de distancias euclidianas
dist_matrix <- dist(sample_data, method = "euclidean")

# Aplicar clustering PAM con 4 clusters
clusters_pam <- pam(dist_matrix, k = 4)

# Visualizar los clusters PAM
fviz_cluster(list(data = sample_data, cluster = clusters_pam$clustering))

# Convertir a matriz numérica
sample_data_matrix <- as.matrix(sample_data)

# Aplicar biclustering Spectral
biclusters_spectral <- biclust(sample_data_matrix, method = BCSpectral())

# Visualizar los biclusters (esta parte debe ser adaptada si quieres hacer un heatmap)
pheatmap(sample_data_matrix)




