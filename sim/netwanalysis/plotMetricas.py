import pandas as pd
import networkx as nx
import numpy as np
from scipy.stats import spearmanr
import itertools
import os
import matplotlib.pyplot as plt
import seaborn as sns

# === 1. CONFIGURACIÓN ===
directorio_script = os.path.dirname(os.path.abspath(__file__))
filepath = '/home/martxoto/Escritorio/TFM/sim/interactions/all_web_interactions.csv'

# 🌟 Cambia esto si tu columna de parches se llama distinto (ej. 'Habitat')
COLUMNA_PARCHE = 'Habitat' 

# === 2. FUNCIÓN NODF ===
def calcular_nodf(G):
    plantas = {n for n, d in G.nodes(data=True) if d.get('bipartite') == 0}
    insectos = {n for n, d in G.nodes(data=True) if d.get('bipartite') == 1}
    if len(plantas) < 2 or len(insectos) < 2: return 0.0
        
    plantas_ord = sorted(list(plantas), key=lambda x: G.degree(x), reverse=True)
    insectos_ord = sorted(list(insectos), key=lambda x: G.degree(x), reverse=True)
    
    try:
        matriz = nx.bipartite.biadjacency_matrix(G, row_order=plantas_ord, column_order=insectos_ord, weight=None).toarray()
    except AttributeError:
        matriz = nx.algorithms.bipartite.matrix.biadjacency_matrix(G, row_order=plantas_ord, column_order=insectos_ord, weight=None).toarray()
    
    matriz = np.where(matriz > 0, 1, 0)
    filas, cols = matriz.shape
    nodf_filas, nodf_cols = 0.0, 0.0
    pares_filas, pares_cols = 0, 0
    
    for i in range(filas - 1):
        for j in range(i + 1, filas):
            grado_i, grado_j = np.sum(matriz[i]), np.sum(matriz[j])
            if grado_i > grado_j and grado_j > 0: 
                nodf_filas += (np.sum(matriz[i] * matriz[j]) / grado_j) * 100
            pares_filas += 1
            
    for i in range(cols - 1):
        for j in range(i + 1, cols):
            grado_i, grado_j = np.sum(matriz[:, i]), np.sum(matriz[:, j])
            if grado_i > grado_j and grado_j > 0: 
                nodf_cols += (np.sum(matriz[:, i] * matriz[:, j]) / grado_j) * 100
            pares_cols += 1
            
    pares_totales = pares_filas + pares_cols
    if pares_totales == 0: return 0.0
    return (nodf_filas + nodf_cols) / pares_totales

# === 3. CARGA DE DATOS Y EXTRACCIÓN AUTOMÁTICA DE SITIOS ===
try:
    df_completo = pd.read_csv(filepath)
except FileNotFoundError:
    print(f"Error: No se encontró el archivo {filepath}")
    exit()

# 🌟 MAGIA: Extraemos todos los sitios únicos del CSV automáticamente
SITIOS_A_COMPARAR = df_completo['Site'].dropna().unique()
print(f"🌍 Se han detectado {len(SITIOS_A_COMPARAR)} ecosistemas en el archivo.")
print("Calculando métricas para todos ellos. Esto puede tardar un poco...")

resultados = []

# === 4. PROCESAMIENTO DE TODOS LOS SITIOS ===
for site in SITIOS_A_COMPARAR:
    df_site = df_completo[(df_completo['Site'] == site) & (df_completo['Conflict'].isna())].copy()
    if df_site.empty:
        continue

    parches = df_site[COLUMNA_PARCHE].unique()
    redes_parches = {}
    
    for p in parches:
        datos_parche = df_site[df_site[COLUMNA_PARCHE] == p]
        G = nx.Graph()
        for _, row in datos_parche.iterrows():
            G.add_node(row['Lower_Taxon'], bipartite=0)
            G.add_node(row['Upper_Taxon'], bipartite=1)
            G.add_edge(row['Lower_Taxon'], row['Upper_Taxon'])
        redes_parches[p] = G
        
    G_global = nx.compose_all(list(redes_parches.values()))
    
    # 1. Conectancia Global
    P = len({n for n, d in G_global.nodes(data=True) if d.get('bipartite') == 0})
    V = len({n for n, d in G_global.nodes(data=True) if d.get('bipartite') == 1})
    C_global = G_global.number_of_edges() / (P * V) if (P * V) > 0 else 0
    
    # 2. NODF Global
    nodf_global = calcular_nodf(G_global)
    
    # Métricas Espaciales (Requieren >= 2 parches)
    cv_asimetria, jaccard_medio, spearman_medio = np.nan, np.nan, np.nan
    
    if len(parches) >= 2:
        # Asimetría
        tamaños = [G.number_of_edges() for G in redes_parches.values()]
        cv_asimetria = np.std(tamaños) / np.mean(tamaños) if np.mean(tamaños) > 0 else 0
        
        # Jaccard y Spearman
        pares = list(itertools.combinations(parches, 2))
        jaccard_vals, spearman_vals = [], []
        
        plantas_glob = {n for n, d in G_global.nodes(data=True) if d.get('bipartite') == 0}
        
        for p1, p2 in pares:
            # Jaccard
            e1 = {tuple(sorted(e)) for e in redes_parches[p1].edges()}
            e2 = {tuple(sorted(e)) for e in redes_parches[p2].edges()}
            union = len(e1.union(e2))
            if union > 0: jaccard_vals.append(len(e1.intersection(e2)) / union)
            
            # Spearman
            G1, G2 = redes_parches[p1], redes_parches[p2]
            p_comunes = list(plantas_glob.intersection(set(G1.nodes())).intersection(set(G2.nodes())))
            if len(p_comunes) > 2:
                corr, _ = spearmanr([G1.degree(n) for n in p_comunes], [G2.degree(n) for n in p_comunes])
                if not np.isnan(corr): spearman_vals.append(corr)
                
        if jaccard_vals: jaccard_medio = np.mean(jaccard_vals)
        if spearman_vals: spearman_medio = np.mean(spearman_vals)

    resultados.append({
        'Ecosistema': site,
        'Conectancia Global': C_global,
        'Anidamiento (NODF)': nodf_global,
        'Asimetría (CV)': cv_asimetria,
        'Redundancia (Jaccard)': jaccard_medio,
        'Estabilidad Hubs (Spearman)': spearman_medio
    })

df_resultados = pd.DataFrame(resultados)

# === 5. CREACIÓN DEL DASHBOARD VISUAL ===
sns.set_theme(style="whitegrid")
# Hacemos la figura un poco más ancha por si hay muchos sitios
fig, axes = plt.subplots(2, 3, figsize=(20, 12)) 
fig.suptitle('Análisis Topológico y Espacial de Todos los Ecosistemas', fontsize=22, fontweight='bold', y=0.98)

metricas = [
    ('Redundancia (Jaccard)', axes[0, 0], '#2ca02c', 'Solapamiento espacial de interacciones (0 a 1)'),
    ('Estabilidad Hubs (Spearman)', axes[0, 1], '#9467bd', 'Correlación de Hubs entre parches (-1 a 1)'),
    ('Asimetría (CV)', axes[0, 2], '#ff7f0e', 'Desequilibrio de tamaño (CV > 0.5 es asimétrico)'),
    ('Conectancia Global', axes[1, 0], '#1f77b4', 'Densidad de la red total (0 a 1)'),
    ('Anidamiento (NODF)', axes[1, 1], '#d62728', 'Estructura Nested / Especialización (0 a 100)')
]

for col, ax, color_hex, desc in metricas:
    # Dibujamos las barras
    sns.barplot(data=df_resultados, x='Ecosistema', y=col, ax=ax, color=color_hex)
    ax.set_ylabel('Valor', fontsize=12)
    ax.set_xlabel('')
    ax.set_title(f"{col}\n{desc}", fontsize=14, fontweight='bold')
    
    # 🌟 Rotamos las etiquetas 45 grados y las alineamos a la derecha para que no se pisen
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)

# Eliminamos el último subplot vacío (2 filas x 3 columnas = 6 huecos, usamos 5)
fig.delaxes(axes[1, 2])

plt.tight_layout()
plt.subplots_adjust(top=0.88) # Dejar espacio para el título principal

# Guardar y mostrar
nombre_imagen = "dashboard_metricas_TODOS.png"
plt.savefig(nombre_imagen, dpi=300, bbox_inches='tight')
print(f"\n✅ ¡Completado! Se han procesado las {len(metricas)} métricas para {len(SITIOS_A_COMPARAR)} sitios.")
print(f"🖼️ Dashboard de alta resolución guardado como: {nombre_imagen}")
plt.show()