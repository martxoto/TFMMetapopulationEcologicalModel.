import pandas as pd
import networkx as nx

# === CONFIGURACIÓN ===
filepath = '/home/martxoto/Escritorio/TFM/sim/interactions/all_web_interactions.csv'
site = 'Woolacombe' # Cambia esto por la red que quieras analizar (ej: Walborough, Hangman_Hill)

# 1. Cargamos y filtramos los datos
df = pd.read_csv(filepath)
df_site = df[(df['Site'] == site) & (df['Conflict'].isna())]

print(f"Procesando {site}... ({len(df_site)} avistamientos empíricos encontrados)")

# 2. Creamos un ÚNICO grafo para representar el ecosistema colapsado (Mónada)
G = nx.Graph()

# 3. Conteo global (ignorando el hábitat/parche original)
for _, row in df_site.iterrows():
    plant = row['Lower_Taxon']
    insect = row['Upper_Taxon']
    
    G.add_node(plant)
    G.add_node(insect)
    
    if G.has_edge(plant, insect):
        G[plant][insect]['weight'] += 1
    else:
        G.add_edge(plant, insect, weight=1)

# 4. Normalización global
if G.number_of_edges() > 0:
    # Buscamos el máximo de TODOS los avistamientos juntos
    weights = [G[u][v]['weight'] for u, v in G.edges()]
    max_weight = max(weights)
    print(f"El peso máximo global encontrado es: {max_weight} observaciones.")
    
    if max_weight > 0:
        for u, v in G.edges():
            G[u][v]['weight'] /= max_weight

# 5. Formatear la salida para el simulador C++
# Usamos un 'set' para sacar todas las plantas rápido y saber quién es quién
all_plants = set(df_site['Lower_Taxon'].unique())
all_interaction_data = []

for u, v in G.edges():
    # Determinar quién es la planta y quién el insecto
    planta, insecto = (u, v) if u in all_plants else (v, u)
    
    all_interaction_data.append({
        'patch_id': 0,  # <-- TODO SE ASIGNA AL PARCHE 0 (Mónada colapsada)
        'Planta': planta,
        'Insecto': insecto,
        'Peso_Normalizado': G[planta][insecto]['weight']
    })

df_out = pd.DataFrame(all_interaction_data)

# 6. Guardar archivo .txt
output_filename = f'interactions/collapsedInteractions/3interactions_{site}_COLLAPSED.txt'
df_out.to_csv(output_filename, sep=' ', index=False, header=False,
              columns=['patch_id', 'Planta', 'Insecto', 'Peso_Normalizado'])

print(f"\n✅ Fichero colapsado guardado en: {output_filename}")
print(f"Total de interacciones únicas en la red: {len(df_out)}")