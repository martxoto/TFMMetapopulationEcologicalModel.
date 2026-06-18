import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

filepath = 'all_web_interactions.csv'
df = pd.read_csv(filepath)

site = 'Hangman_Hill'
df_site = df[(df['Site'] == site) & (df['Conflict'].isna())]

habitats = df_site['Habitat'].unique()
habitats.sort()
numpatch = len(habitats)

patch_map = {habitat_name: i for i, habitat_name in enumerate(habitats)}

print(patch_map)

G_list = [nx.Graph() for _ in range(numpatch)]


for _, row in df_site.iterrows():
    
    plant = row['Lower_Taxon']
    insect = row['Upper_Taxon']
    habitat_name = row['Habitat']
	
    patch_id = patch_map[habitat_name]
	
    G = G_list[patch_id] 
	
    G.add_node(plant)
    G.add_node(insect)
    if G.has_edge(plant, insect):
        G[plant][insect]['weight'] += 1
    else:
        G.add_edge(plant, insect, weight=1)

for G in G_list:
	if G.number_of_edges() > 0:
	    weights = [G[u][v]['weight'] for u, v in G.edges()]
	    max_weight = max(weights)
	    if max_weight > 0:
		    for u, v in G.edges():
		        G[u][v]['weight'] /= max_weight

# === 6. Guardar TODO en un solo archivo de 4 columnas ===
#    (Usaremos un .txt separado por espacios, es más fácil para C++)

all_interaction_data = []
all_plants = set() # Para identificar plantas vs insectos
all_insects = set()

for patch_id, G in enumerate(G_list):
    # Necesitamos saber quién es planta y quién insecto
    # (Tu método de comprobar 'u in df_site['Lower_Taxon'].values' es lento)
    # Vamos a pre-calcularlos (simplificado)
    for node in G.nodes():
        if df_site[df_site['Lower_Taxon'] == node].shape[0] > 0:
            all_plants.add(node)
        else:
            all_insects.add(node)

for patch_id, G in enumerate(G_list):
    for u, v in G.edges():
        # Determinar planta e insecto
        planta, insecto = (u, v) if u in all_plants else (v, u)
        
        all_interaction_data.append({
            'patch_id': patch_id,
            'Planta': planta,
            'Insecto': insecto,
            'Peso_Normalizado': G[planta][insecto]['weight']
        })

df_out = pd.DataFrame(all_interaction_data)

# Guardar como .txt separado por espacios y sin cabecera
output_filename = f'interactions_{site}_patches.txt'
df_out.to_csv(output_filename, sep=' ', index=False, header=False,
              columns=['patch_id', 'Planta', 'Insecto', 'Peso_Normalizado'])

print(f"Fichero de interacciones guardado en: {output_filename}")
'''
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
'''