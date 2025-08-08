import pandas as pd
import networkx as nx
import numpy as np
from scipy.linalg import expm

# === Cargar datos ===
filepath = 'all_web_interactions.csv'
df = pd.read_csv(filepath)

# === Filtrar por sitio y sin conflictos ===
site = 'Penhale_Sands'
df_site = df[(df['Site'] == site) & (df['Conflict'].isna())]

# === Crear grafo ===
G = nx.Graph()
for _, row in df_site.iterrows():
    planta = row['Lower_Taxon']
    insecto = row['Upper_Taxon']
    G.add_node(planta)
    G.add_node(insecto)
    G.add_edge(planta, insecto)

# === Usar componente conexa más grande si no es conexo ===
if not nx.is_connected(G):
    G_lcc = G.subgraph(max(nx.connected_components(G), key=len)).copy()
else:
    G_lcc = G

# === Métricas estructurales (solo LCC para lo necesario) ===
shortest_paths = dict(nx.all_pairs_shortest_path_length(G_lcc))
distances = [dist for sp in shortest_paths.values() for dist in sp.values()]
mean_path_length = np.mean(distances)
std_dev = np.std(distances, ddof=1)
sem = std_dev / np.sqrt(len(distances))

A = nx.adjacency_matrix(G).todense()
A_lcc = nx.adjacency_matrix(G_lcc).todense()

def compute_be(A):
    exp_A = expm(A)
    exp_neg_A = expm(-A)
    return np.trace(exp_neg_A) / np.trace(exp_A)

# === Centralidad (se calculan sobre todo G) ===
def top_n_central_nodes(centrality_dict, n=10):
    return sorted(centrality_dict.items(), key=lambda x: x[1], reverse=True)[:n]

degree_centrality = nx.degree_centrality(G)
closeness_centrality = nx.closeness_centrality(G)
betweenness_centrality = nx.betweenness_centrality(G)
eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000)
katz_centrality = nx.katz_centrality(G)
pagerank = nx.pagerank(G)
subgraph_centrality = nx.subgraph_centrality(G)

# === Guardar resultados ===
filename = f"{site}_structure.txt".replace(" ", "_")

with open(filename, 'w', encoding='utf-8') as f:
    f.write("=== Métricas estructurales ===\n")
    f.write(f"Density: {nx.density(G):.4f}\n")
    f.write(f"Average clustering: {nx.average_clustering(G):.4f}\n")
    f.write(f"Transitivity: {nx.transitivity(G):.4f}\n")
    f.write(f"Average path length (LCC): {mean_path_length:.4f}\n")
    f.write(f"Standard deviation of path length (LCC): {std_dev:.4f}\n")
    f.write(f"Diameter (LCC): {nx.diameter(G_lcc)}\n")
    f.write(f"Degree assortativity: {nx.degree_pearson_correlation_coefficient(G):.4f}\n")
    f.write(f"Bipartivity index (LCC): {compute_be(A_lcc):.4f}\n\n")

    f.write("=== Top 10 por centralidad ===\n\n")
    def write_top_centrality(name, centrality_dict):
        f.write(f"{name}:\n")
        for node, value in top_n_central_nodes(centrality_dict):
            f.write(f"  {node}: {value:.4f}\n")
        f.write("\n")

    write_top_centrality("Degree Centrality", degree_centrality)
    write_top_centrality("Closeness Centrality", closeness_centrality)
    write_top_centrality("Betweenness Centrality", betweenness_centrality)
    write_top_centrality("Eigenvector Centrality", eigenvector_centrality)
    write_top_centrality("Katz Centrality", katz_centrality)
    write_top_centrality("PageRank", pagerank)
    write_top_centrality("Subgraph Centrality", subgraph_centrality)

print(f"✅ Resultados guardados en: {filename}")
