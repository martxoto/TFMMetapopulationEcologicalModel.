import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# === CONFIGURACIÓN ===
# Ajusta la ruta a tu archivo concreto
INTERACTIONS_FILE = 'interactions/interactions_Penhale_Sands_patches.txt' 
EVO_P_FILE = 'evolutionp.txt' 
EVO_V_FILE = 'evolutionv.txt'
VIABILITY_THRESHOLD = 1.0e-6  # <--- Umbral de corte

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.size': 12})

def get_degrees(filename):
    """Calcula el grado (k) de cada especie desde el archivo de interacciones"""
    try:
        df = pd.read_csv(filename, sep=' ', header=None, names=['Patch', 'Plant', 'Insect', 'Weight'])
        # Contamos conexiones únicas
        k_plants = df.groupby('Plant')['Insect'].nunique()
        k_insects = df.groupby('Insect')['Plant'].nunique()
        return k_plants, k_insects
    except FileNotFoundError:
        print(f"Error: No encuentro {filename}")
        return None, None

def get_abundances(filename):
    """Lee la última línea del fichero de evolución temporal"""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            if not lines: return None
            # Última línea = Estado Estacionario
            data = [float(x) for x in lines[-1].strip().split()]
            # La primera columna es el tiempo (t), la quitamos
            return np.array(data[1:])
    except FileNotFoundError:
        print(f"Error: No encuentro {filename}")
        return None

# 1. Cargar Datos
print("Cargando datos...")

k_p_series, k_v_series = get_degrees(INTERACTIONS_FILE)
raw_p = get_abundances(EVO_P_FILE)
raw_v = get_abundances(EVO_V_FILE)

if k_p_series is None or raw_p is None or raw_v is None:
    print("Faltan archivos. Abortando.")
    exit()

# Reconstruir biomasa total por especie sumando parches
num_plants = len(k_p_series)
num_patches_p = len(raw_p) // num_plants
biomass_p = raw_p.reshape((num_plants, num_patches_p)).sum(axis=1)

num_insects = len(k_v_series)
num_patches_v = len(raw_v) // num_insects
biomass_v = raw_v.reshape((num_insects, num_patches_v)).sum(axis=1)

# Alineamos los grados con el orden del vector (asumimos orden alfabético/índice)
k_p_sorted = k_p_series.sort_index().values
k_v_sorted = k_v_series.sort_index().values

# === 2. FILTRADO (LA PARTE NUEVA) ===
print(f"\nFiltrando especies con biomasa < {VIABILITY_THRESHOLD}...")

# Máscara Plantas
mask_p = biomass_p > VIABILITY_THRESHOLD
biomass_p_clean = biomass_p[mask_p]
k_p_clean = k_p_sorted[mask_p]
print(f" -> Plantas Vivas: {len(biomass_p_clean)} / {len(biomass_p)}")

# Máscara Insectos
mask_v = biomass_v > VIABILITY_THRESHOLD
biomass_v_clean = biomass_v[mask_v]
k_v_clean = k_v_sorted[mask_v]
print(f" -> Insectos Vivos: {len(biomass_v_clean)} / {len(biomass_v)}")

if len(biomass_p_clean) == 0 or len(biomass_v_clean) == 0:
    print("ADVERTENCIA: Todo está muerto. No se pueden pintar gráficas logarítmicas.")
    exit()

# === 3. GRÁFICAS ===
fig, axs = plt.subplots(2, 2, figsize=(14, 10))

# A) Rank-Abundance (Plantas) - LogScale
sorted_p = np.sort(biomass_p_clean)[::-1] # Ordenamos de mayor a menor SOLO las vivas
axs[0, 0].plot(range(1, len(sorted_p)+1), sorted_p, 'o-', color='forestgreen')
axs[0, 0].set_yscale('log')
axs[0, 0].set_xscale('log')
axs[0, 0].set_title('Rank-Abundance (Plantas Vivas)')
axs[0, 0].set_ylabel('Biomasa (log)')
axs[0, 0].set_xlabel('Rango (log)')

# B) Rank-Abundance (Insectos) - LogScale
sorted_v = np.sort(biomass_v_clean)[::-1]
axs[0, 1].plot(range(1, len(sorted_v)+1), sorted_v, 'o-', color='darkorange')
axs[0, 1].set_yscale('log')
axs[0, 1].set_xscale('log')
axs[0, 1].set_title('Rank-Abundance (Insectos Vivos)')
axs[0, 1].set_ylabel('Biomasa (log)')
axs[0, 1].set_xlabel('Rango (log)')

# C) Correlación Grado-Abundancia (Plantas)
# Usamos las versiones _clean para mantener la correspondencia k <-> biomasa
axs[1, 0].scatter(k_p_clean, biomass_p_clean, color='forestgreen', alpha=0.6, edgecolors='k')
axs[1, 0].set_yscale('log')
axs[1, 0].set_xscale('log')
axs[1, 0].set_title('Correlación Grado-Abundancia (Plantas)')
axs[1, 0].set_ylabel('Biomasa (log)')
axs[1, 0].set_xlabel('Grado k (log)')

# D) Correlación Grado-Abundancia (Insectos)
axs[1, 1].scatter(k_v_clean, biomass_v_clean, color='darkorange', alpha=0.6, edgecolors='k')
axs[1, 1].set_yscale('log')
axs[1, 1].set_xscale('log')
axs[1, 1].set_title('Correlación Grado-Abundancia (Insectos)')
axs[1, 1].set_ylabel('Biomasa (log)')
axs[1, 1].set_xlabel('Grado k (log)')

plt.tight_layout()
plt.show()