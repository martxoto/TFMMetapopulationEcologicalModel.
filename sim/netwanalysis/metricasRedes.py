import pandas as pd
import networkx as nx
import numpy as np
from scipy.stats import spearmanr
import itertools
import os

# === 1. CONFIGURACIÓN ===
directorio_script = os.path.dirname(os.path.abspath(__file__))
filepath = '/home/martxoto/Escritorio/TFM/sim/interactions/all_web_interactions.csv'

SITE = 'Studland'
COLUMNA_PLANTA = 'Lower_Taxon'
COLUMNA_INSECTO = 'Upper_Taxon'
COLUMNA_PARCHE = 'Habitat'

# === 2. CARGA Y FILTRADO ===
try:
    df = pd.read_csv(filepath)
except FileNotFoundError:
    print(f"Error: No se encontró el archivo {filepath}")
    exit()

df_site = df[(df['Site'] == SITE) & (df['Conflict'].isna())].copy()

if df_site.empty:
    print(f"Error: No se encontraron datos para el sitio '{SITE}'.")
    exit()

parches = df_site[COLUMNA_PARCHE].unique()
print(f"\n{'='*50}")
print(f"🌍 PERFIL ECOLÓGICO DE: {SITE}")
print(f"Número de parches locales: {len(parches)} ({', '.join(map(str, parches))})")
print(f"{'='*50}\n")

# === 3. CONSTRUCCIÓN DE REDES POR PARCHE ===
redes_parches = {}
for p in parches:
    datos_parche = df_site[df_site[COLUMNA_PARCHE] == p]
    G = nx.Graph()
    for _, row in datos_parche.iterrows():
        planta = row[COLUMNA_PLANTA]
        insecto = row[COLUMNA_INSECTO]
        G.add_node(planta, bipartite=0)
        G.add_node(insecto, bipartite=1)
        G.add_edge(planta, insecto)
    redes_parches[p] = G

# === 4. CÁLCULO DE MÉTRICAS ===

# A. Conectancia (Local y Global)
print("--- 1. CONECTANCIA (C) ---")
# Global
G_global = nx.compose_all(list(redes_parches.values()))
plantas_glob = {n for n, d in G_global.nodes(data=True) if d.get('bipartite') == 0}
insectos_glob = {n for n, d in G_global.nodes(data=True) if d.get('bipartite') == 1}
C_global = G_global.number_of_edges() / (len(plantas_glob) * len(insectos_glob))
print(f"Global (Mónada):\t {C_global:.4f}")

for p, G in redes_parches.items():
    P = len({n for n, d in G.nodes(data=True) if d.get('bipartite') == 0})
    V = len({n for n, d in G.nodes(data=True) if d.get('bipartite') == 1})
    E = G.number_of_edges()
    C_local = E / (P * V) if (P * V) > 0 else 0
    print(f"Parche {p}:\t\t {C_local:.4f}")

# B. Asimetría de Tamaño (Efecto Fuente-Sumidero)
print("\n--- 2. ASIMETRÍA DE TAMAÑO ---")
tamaños = [G.number_of_edges() for G in redes_parches.values()]
if len(tamaños) > 1:
    desviacion = np.std(tamaños)
    media = np.mean(tamaños)
    cv = desviacion / media # Coeficiente de variación
    print(f"Interacciones por parche: {tamaños}")
    print(f"Asimetría (CV):\t\t {cv:.4f} (0 = Simetría perfecta, >0.5 = Muy asimétrico)")
else:
    print("Solo hay un parche, no hay asimetría espacial.")

# C. Beta-Diversidad de Interacciones (Índice de Jaccard)
print("\n--- 3. REDUNDANCIA ESPACIAL (JACCARD) ---")
if len(parches) >= 2:
    pares = list(itertools.combinations(parches, 2))
    jaccard_vals = []
    for p1, p2 in pares:
        edges1 = set(redes_parches[p1].edges())
        edges2 = set(redes_parches[p2].edges())
        # Normalizar el orden de las tuplas para que A-B y B-A coincidan
        edges1 = {tuple(sorted(e)) for e in edges1}
        edges2 = {tuple(sorted(e)) for e in edges2}
        
        interseccion = len(edges1.intersection(edges2))
        union = len(edges1.union(edges2))
        J = interseccion / union if union > 0 else 0
        jaccard_vals.append(J)
        print(f"Jaccard {p1} vs {p2}:\t {J:.4f} ({interseccion} enlaces compartidos de {union})")
    
    print(f"Jaccard Medio:\t\t {np.mean(jaccard_vals):.4f} (0 = Aislamiento total, 1 = Redundancia total)")
else:
    print("Se requieren al menos 2 parches para medir la beta-diversidad.")

# D. Correlación del Grado (El Rol del Hub)
print("\n--- 4. ESTABILIDAD DE LOS HUBS (CORRELACIÓN DE SPEARMAN) ---")
if len(parches) >= 2:
    for p1, p2 in pares:
        G1, G2 = redes_parches[p1], redes_parches[p2]
        plantas_comunes = list(plantas_glob.intersection(set(G1.nodes())).intersection(set(G2.nodes())))
        
        if len(plantas_comunes) > 2:
            grados_1 = [G1.degree(n) for n in plantas_comunes]
            grados_2 = [G2.degree(n) for n in plantas_comunes]
            corr, p_val = spearmanr(grados_1, grados_2)
            print(f"Spearman {p1} vs {p2}:\t {corr:.4f} (p-value: {p_val:.3f})")
            if corr > 0.7 and p_val < 0.05:
                print(" -> El Hub es el mismo en ambos parches (Estructura centralizada).")
            elif corr < 0.3:
                print(" -> Los Hubs cambian según el parche (Estructura descentralizada).")
        else:
            print(f"Spearman {p1} vs {p2}:\t No hay suficientes plantas en común para correlacionar.")
else:
    print("Se requieren al menos 2 parches para medir la correlación espacial.")

# E. Anidamiento (NODF - Nestedness metric based on Overlap and Decreasing Fill)
print("\n--- 5. ANIDAMIENTO (NODF) ---")

def calcular_nodf(G):
    # Extraemos plantas e insectos
    plantas = {n for n, d in G.nodes(data=True) if d.get('bipartite') == 0}
    insectos = {n for n, d in G.nodes(data=True) if d.get('bipartite') == 1}
    
    if len(plantas) < 2 or len(insectos) < 2:
        return 0.0 # Redes muy pequeñas no tienen anidamiento medible
        
    # Ordenamos nodos de mayor a menor grado (Decreasing Fill)
    plantas_ord = sorted(list(plantas), key=lambda x: G.degree(x), reverse=True)
    insectos_ord = sorted(list(insectos), key=lambda x: G.degree(x), reverse=True)
    
    # Creamos la matriz de adyacencia
    try:
        matriz = nx.bipartite.biadjacency_matrix(G, row_order=plantas_ord, column_order=insectos_ord, weight=None).toarray()
    except AttributeError:
        # Por compatibilidad con distintas versiones de NetworkX
        matriz = nx.algorithms.bipartite.matrix.biadjacency_matrix(G, row_order=plantas_ord, column_order=insectos_ord, weight=None).toarray()
    
    # MUY IMPORTANTE: Nos aseguramos de que sea estrictamente binaria (1 o 0)
    matriz = np.where(matriz > 0, 1, 0)
    
    filas, cols = matriz.shape
    nodf_filas = 0.0
    pares_filas = 0
    
    # Solapamiento de filas (Plantas)
    for i in range(filas - 1):
        for j in range(i + 1, filas):
            grado_i, grado_j = np.sum(matriz[i]), np.sum(matriz[j])
            
            # FIX: Solo calculamos si el grado de j es estrictamente mayor que 0
            if grado_i > grado_j and grado_j > 0: 
                interseccion = np.sum(matriz[i] * matriz[j])
                nodf_filas += (interseccion / grado_j) * 100
            pares_filas += 1
            
    nodf_cols = 0.0
    pares_cols = 0
    
    # Solapamiento de columnas (Insectos)
    for i in range(cols - 1):
        for j in range(i + 1, cols):
            grado_i, grado_j = np.sum(matriz[:, i]), np.sum(matriz[:, j])
            
            # FIX: Solo calculamos si el grado de j es estrictamente mayor que 0
            if grado_i > grado_j and grado_j > 0: 
                interseccion = np.sum(matriz[:, i] * matriz[:, j])
                nodf_cols += (interseccion / grado_j) * 100
            pares_cols += 1
            
    pares_totales = pares_filas + pares_cols
    if pares_totales == 0: return 0.0
    
    return (nodf_filas + nodf_cols) / pares_totales

# Calculamos NODF para la red global y por parche
nodf_global = calcular_nodf(G_global)
print(f"Global (Mónada):\t {nodf_global:.2f} / 100")

for p, G in redes_parches.items():
    nodf_local = calcular_nodf(G)
    print(f"Parche {p}:\t\t {nodf_local:.2f} / 100")

print(f"\n{'='*50}")