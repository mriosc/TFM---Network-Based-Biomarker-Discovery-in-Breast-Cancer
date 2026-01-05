import networkx as nx
import csv

# Carga de datos (Marc)
normalSarcoma = nx.read_edgelist('C:/Users/marcr/OneDrive/Escritorio/Master_Bioinf/Semestre_4/TFM/Resultados/06PrimaryNetwork_NameInteraction.csv', delimiter=(';'), nodetype=str)
tumorSarcoma = nx.read_edgelist('C:/Users/marcr/OneDrive/Escritorio/Master_Bioinf/Semestre_4/TFM/Resultados/06MetastasicNetwork_NameInteraction.csv', delimiter=(';'), nodetype=str)

# Invertir nodos del normal
invertedNormal = [(target, source) for source, target in normalSarcoma.edges()]
invertedNormal = nx.Graph(invertedNormal)

# Restar tumor menos el normal
interactionTumorSarcomaExclussive = set(tumorSarcoma.edges()) - set(normalSarcoma.edges()) # interacciones exclusivas tumorales
tumorSarcomaExclussive = nx.Graph(list(interactionTumorSarcomaExclussive))

# Restar resultado anterior menos el invertido
interactionTumorSarcomaExclussive = set(tumorSarcomaExclussive.edges()) - set(invertedNormal.edges()) # interacciones exclusivas tumorales
tumorSarcomaExclussive = nx.Graph(list(interactionTumorSarcomaExclussive))

# Extraer las aristas de tumorSarcomaExclussive
edges_tumor_exclusive = list(tumorSarcomaExclussive.edges())


# Guardar las aristas en un fichero CSV (Marc)
with open('C:/Users/marcr/OneDrive/Escritorio/Master_Bioinf/Semestre_4/TFM/Resultados/NetworktumorExclussive.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["GEN1", "GEN2"])
    writer.writerows(edges_tumor_exclusive)
