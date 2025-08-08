import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

filepath = 'all_web_interactions.csv'
df = pd.read_csv(filepath)

site = 'Penhale_Sands'
df_site = df[(df['Site'] == site) & (df['Conflict'].isna())]

G = nx.Graph()
for _, row in df_site.iterrows():
    planta = row['Lower_Taxon']
    insecto = row['Upper_Taxon']
    G.add_node(planta, type='plant')
    G.add_node(insecto, type='insect')
    G.add_edge(planta, insecto)

plantas = [n for n, d in G.nodes(data=True) if d['type'] == 'plant']
insectos = [n for n, d in G.nodes(data=True) if d['type'] == 'insect']

plant_degrees = [G.degree(p) for p in plantas]
insect_degrees = [G.degree(i) for i in insectos]

max_k = max(max(plant_degrees), max(insect_degrees))

plant_hist = np.zeros(max_k + 1)
insect_hist = np.zeros(max_k + 1)

for k in plant_degrees:
    plant_hist[k] += 1
for k in insect_degrees:
    insect_hist[k] += 1

plant_hist /= len(plantas)
insect_hist /= len(insectos)

degrees = np.arange(max_k + 1)

plt.figure(figsize=(10, 6))
plt.bar(degrees - 0.2, plant_hist, width=0.4, label='Plants', color='green', alpha=0.6, edgecolor='black')
plt.bar(degrees + 0.2, insect_hist, width=0.4, label='Insects', color='orange', alpha=0.6, edgecolor='black')
plt.xlabel('k')
plt.ylabel('P(k)')
plt.title(f'Degree distribution in {site}')
plt.legend()
plt.gca().spines['top'].set_color('none')
plt.gca().spines['right'].set_color('none')
plt.gca().spines['left'].set_color('none')
plt.gca().spines['bottom'].set_color('none')
plt.tight_layout()
plt.show()
