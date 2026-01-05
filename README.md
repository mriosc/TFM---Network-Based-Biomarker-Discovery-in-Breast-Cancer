# Network-Based Biomarker Discovery in Breast Cancer

Este repositorio contiene el flujo de trabajo bioinformático desarrollado para el Trabajo de Fin de Máster (TFM) enfocado en la identificación de biomarcadores para la metástasis en cáncer de mama mediante el análisis de redes de coexpresión diferencial.

## Descripción del Proyecto
El proyecto emplea una metodología innovadora basada en la **sustracción de grafos**. A partir de datos transcriptómicos de tumores primarios y metastásicos, se construyen redes de coexpresión y se aíslan las interacciones exclusivas del estado metastásico para identificar reguladores maestros (Hubs) integrados en módulos biológicamente coherentes.

## Flujo de Trabajo
1. **Preprocesamiento y DEGs**: Descarga de datos de GEO, normalización y análisis de expresión diferencial (método TREAT).
2. **Construcción de Redes**: Generación de redes de coexpresión basadas en 5000 HVGs (Highly Variable Genes).
3. **Sustracción de Redes**: Generación de la red exclusiva de tumor mediante la resta de interacciones comunes.
4. **Análisis Modular**: Clustering (GLay) y filtrado de comunidades por masa crítica (>10 genes).
5. **Identificación de Hubs**: Análisis topológico mediante centralidad de grado ($Q3 + 1.5 \times IQR$).
6. **Enriquecimiento Funcional**: Anotación biológica mediante Gene Ontology (GO) y KEGG.

## Estructura de Archivos
- `1preprocess_QC_degs.r`: Scripts en R para control de calidad y expresión diferencial.
- `2restaGrafos.py`: Script en Python para la sustracción lógica de redes.
- `3ClustersMax10Genes.py`: Filtrado de clústeres por tamaño.
- `4hubExploration.py`: Identificación estadística de nodos Hub.
- `5enrichment_auto.R`: Automatización de análisis funcional.

## Requisitos e Instalación
Para replicar este análisis, clona el repositorio y asegura tener instalados los requisitos de Python y R:

```bash
git clone [https://github.com/mriosc/TFM---Network-Based-Biomarker-Discovery-in-Breast-Cancer.git](https://github.com/mriosc/TFM---Network-Based-Biomarker-Discovery-in-Breast-Cancer.git)
pip install -r requirements.txt