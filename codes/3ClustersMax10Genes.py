import pandas as pd
import os


# Load the CSV file Marc
df = pd.read_csv('C:/Users/marcr/OneDrive/Escritorio/Master_Bioinf/Semestre_4/TFM/Resultados/NetworkClusters06.csv') # Marc Git
output_path = 'C:/Users/marcr/OneDrive/Escritorio/Master_Bioinf/Semestre_4/TFM/Resultados/Filtered_Clusters.csv' # Marc Git
output_txt_dir = 'C:/Users/marcr/OneDrive/Escritorio/Master_Bioinf/Semestre_4/TFM/Resultados/inputClusters' # Marc Git


# Display the first few rows of the dataframe to understand its structure
df.head()

# Agrupar por la columna '__glayCluster' y eliminar grupos con menos de 10 filas
filtered_df = df.groupby('__glayCluster').filter(lambda x: len(x) >= 10)

# Guardar el dataframe filtrado en un nuevo archivo CSV

filtered_df.to_csv(output_path, index=False)

print(f"Archivo guardado en: {output_path}")

# Guardar un archivo .txt por cada cluster, con los valores de la columna 'name'
for cluster_value, group in filtered_df.groupby('__glayCluster'):
    # Obtener los valores de la columna 'name'
    cluster_names = group['name']
    
    # Crear el nombre del archivo .txt
    txt_file_path = os.path.join(output_txt_dir, f'Cluster_{cluster_value}.txt')
    
    # Guardar los nombres en el archivo .txt sin cabecera
    cluster_names.to_csv(txt_file_path, index=False, header=False)

print(f"Archivos guardados en: {output_txt_dir}")