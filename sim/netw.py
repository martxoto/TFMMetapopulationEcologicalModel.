import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

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
    if G.has_edge(planta, insecto):
        G[planta][insecto]['weight'] += 1
    else:
        G.add_edge(planta, insecto, weight=1)

# === Normalizar pesos ===
weights = [G[u][v]['weight'] for u, v in G.edges()]
max_weight = max(weights)
for u, v in G.edges():
    G[u][v]['weight'] /= max_weight

# === Guardar archivo CSV con interacciones ===
interaction_data = []
for u, v in G.edges():
    planta, insecto = (u, v) if u in df_site['Lower_Taxon'].values else (v, u)
    interaction_data.append({
        'Planta': planta,
        'Insecto': insecto,
        'Peso_Normalizado': G[planta][insecto]['weight'],
        'Peso_Bruto': int(G[planta][insecto]['weight'] * max_weight)
    })

df_out = pd.DataFrame(interaction_data)
df_out.to_csv(f'interactions{site}.csv', index=False)

# === Dibujar grafo ===
plt.figure(figsize=(12, 12))
layout = nx.spring_layout(G, seed=42)
normalized_weights = [G[u][v]['weight'] for u, v in G.edges()]

nx.draw_networkx_nodes(G, pos=layout, node_color='skyblue', node_size=300)
nx.draw_networkx_labels(G, pos=layout, font_size=8)
nx.draw_networkx_edges(G, pos=layout, width=[30 * np.sqrt(w) for w in normalized_weights], edge_color='gray')

plt.title(f"Red de interacciones normalizada en {site}", fontsize=14)
plt.axis('off')
plt.tight_layout()
plt.show()
