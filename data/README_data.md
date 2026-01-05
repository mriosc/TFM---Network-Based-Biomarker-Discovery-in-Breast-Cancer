# Instrucciones de Datos

Para ejecutar el flujo de trabajo, los datos de expresión deben obtenerse siguiendo estos pasos:

1. **Fuente**: Los datos primarios se descargan automáticamente mediante el script `1preprocess_QC_degs.r` desde el repositorio **Gene Expression Omnibus (GEO)**.
2. **Series Utilizadas**: El dataset final "CLEAN" es una integración de múltiples series de cáncer de mama (asegúrate de referenciar las GSE IDs aquí).
3. **Archivos Generados**:
   - `expr_clean_HVGs5000.txt`: Matriz de expresión de los 5000 genes con mayor variabilidad.
   - `06_TREAT_FC1_significant.csv`: Resultados del análisis diferencial.
4. **Estructura**: Los scripts esperan encontrar las matrices de expresión en formato tabular (.txt o .csv) con los símbolos de los genes en la primera columna.