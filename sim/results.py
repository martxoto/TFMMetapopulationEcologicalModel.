import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === 1. CONFIGURACIÓN ===
FILENAME = "resultsTarget.txt" 
OUTPUT_IMAGE = "resultados.png"

# === 2. CARGA DE DATOS ===
try:
    data = np.loadtxt(FILENAME)
    if data.size == 0:
        print(f"Error: El archivo '{FILENAME}' está vacío.")
        exit()
    print(f"Datos cargados correctamente de '{FILENAME}'. {len(data)} pasos detectados.")
except Exception as e:
    print(f"Error al leer los datos: {e}")
    exit()

# === 3. ASIGNACIÓN DE COLUMNAS A VARIABLES ===
k_eff = data[:, 0]
robustness = data[:, 1]
surv_plants = data[:, 2]
surv_insects = data[:, 3]
pollination = data[:, 4]
shannon_p = data[:, 5]
shannon_v = data[:, 6]
norm_biomass = data[:, 7] 

# === 4. GENERACIÓN DEL PANEL (2 filas x 2 columnas) ===
sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10), sharex=True)
colors = sns.color_palette("muted")

# --- GRÁFICA [0, 0]: ROBUSTEZ Y BIODIVERSIDAD (DOBLE EJE) ---
ax1 = axes[0, 0]
color_rob = colors[0]

# Pintamos la robustez en el eje Y principal (izquierdo)
line1 = ax1.plot(k_eff, robustness, '-', color=color_rob, linewidth=2.5, label="Robustez (R)")
ax1.set_ylabel("Radio de Robustez ($R$)", color=color_rob, fontsize=11, fontweight='bold')
ax1.tick_params(axis='y', labelcolor=color_rob)
ax1.set_ylim([0, 1.05])

# Creamos el segundo eje Y que comparte el mismo eje X
ax2 = ax1.twinx()  
color_plants = colors[2] # Típicamente verde en seaborn muted
color_insects = colors[1] # Típicamente naranja

# Pintamos plantas e insectos en el eje Y secundario (derecho)
line2 = ax2.plot(k_eff, surv_plants, '--', color=color_plants, linewidth=2, label="Plantas Vivas")
line3 = ax2.plot(k_eff, surv_insects, ':', color=color_insects, linewidth=2, label="Insectos Vivos")
ax2.set_ylabel("Número de Especies", color='gray', fontsize=11, fontweight='bold')
ax2.tick_params(axis='y', labelcolor='gray')
ax2.grid(False) # Quitamos la rejilla secundaria para evitar cruces de líneas raros

# Unificamos las leyendas de ambos ejes en una sola cajita
lines = line1 + line2 + line3
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='best', fontsize=10)

ax1.set_title("Robustez y Biodiversidad", fontsize=13, fontweight='bold')

# --- GRÁFICA [0, 1]: ÍNDICES DE GINI ---
ax_gini = axes[0, 1]
ax_gini.plot(k_eff, shannon_p, '-', color=colors[2], linewidth=2, label='Shannon Plantas')
ax_gini.plot(k_eff, shannon_v, '-', color=colors[1], linewidth=2, label='Shannon Insectos')
ax_gini.set_title("Entropía de Shannon", fontsize=13, fontweight='bold')
ax_gini.set_ylabel("Entropía de Shannon", fontsize=11)
ax_gini.legend(loc='best', fontsize=10)
ax_gini.set_ylim([0, 4])

# --- GRÁFICA [1, 0]: SERVICIO DE POLINIZACIÓN ---
ax_pol = axes[1, 0]
ax_pol.plot(k_eff, pollination, '-', color=colors[3], linewidth=2.5)
ax_pol.set_title("Servicio de Polinización", fontsize=13, fontweight='bold')
ax_pol.set_xlabel("Extinciones Efectivas ($k_{eff}$)", fontsize=11)
ax_pol.set_ylabel("Biomasa Insectos ($P_S$)", fontsize=11)

# --- GRÁFICA [1, 1]: BIOMASA NORMALIZADA ---
ax_bio = axes[1, 1]
ax_bio.plot(k_eff, norm_biomass, '-', color=colors[4], linewidth=2.5)
ax_bio.set_title("Biomasa Total Normalizada", fontsize=13, fontweight='bold')
ax_bio.set_xlabel("Extinciones Efectivas ($k_{eff}$)", fontsize=11)
ax_bio.set_ylabel("Biomasa Total ($B_{norm}$)", fontsize=11)
ax_bio.set_ylim([0, 1.05])

# === 5. AJUSTES FINALES Y GUARDADO ===
plt.tight_layout()
plt.savefig(OUTPUT_IMAGE, dpi=300)
print(f"Panel guardado como '{OUTPUT_IMAGE}'.")
plt.show()