import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# === CONFIGURACIÓN ===
INPUT_CSV = 'datos_barrido_D.csv'
OUTPUT_PDF = 'curvas_robustez_macroscopica.pdf'
COLOR_TRAYECTORIAS = 'gray'
COLOR_MEDIA = 'darkred'

# 1. Cargar los datos
df = pd.read_csv(INPUT_CSV)

# Separar el vector de dispersión (D) de las ejecuciones (Run_X)
D_vals = df['D']
trayectorias = df.drop(columns=['D'])

# 2. Configuración de estilo académico
plt.style.use('seaborn-v0_8-paper')
sns.set_context("paper", font_scale=1.5)
fig, ax = plt.subplots(figsize=(10, 6))

# 3. Dibujar el fondo estocástico (Spaghetti plot)
# Iteramos sobre cada columna Run y la pintamos fina y transparente
for col in trayectorias.columns:
    ax.plot(D_vals, trayectorias[col], color=COLOR_TRAYECTORIAS, alpha=0.3, linewidth=1.5)

# 4. Calcular el comportamiento macroscópico (Media y Desviación)
media_R = trayectorias.mean(axis=1)
std_R = trayectorias.std(axis=1)

# Dibujar la curva maestra
ax.plot(D_vals, media_R, color=COLOR_MEDIA, linewidth=3.5, marker='o', markersize=6, label='Media del Ensamble')

# Sombrear la desviación estándar (el "ruido" topológico)
ax.fill_between(D_vals, media_R - std_R, media_R + std_R, color=COLOR_MEDIA, alpha=0.15, label='Desviación Estándar ($\pm 1\sigma$)')

# 5. Formato físico y matemático del lienzo
# Usamos symlog porque tenemos el valor D=0.0 que no existe en logaritmo puro
ax.set_xscale('symlog', linthresh=0.01) 
ax.set_xlabel('Tasa de Dispersión ($D$)', fontsize=14)
ax.set_ylabel('Robustez del Ecosistema ($R$)', fontsize=14)
ax.set_title('Transición de la Robustez frente al Acoplamiento Espacial', fontsize=16, pad=15)

# Forzar ticks en el eje X para que coincidan con tus valores importantes
ticks_importantes = [0.0, 0.01, 0.1, 1.0, 10.0, 100.0]
ax.set_xticks(ticks_importantes)
ax.set_xticklabels([str(t) for t in ticks_importantes])

ax.legend(loc='lower right', frameon=True, shadow=True)
ax.grid(True, which="major", ls="-", alpha=0.3)
ax.grid(True, which="minor", ls="--", alpha=0.1)

# 6. Guardar y mostrar
plt.tight_layout()
plt.savefig(OUTPUT_PDF, dpi=300, format='pdf')
print(f"Gráfica exportada con calidad de publicación a: {OUTPUT_PDF}")
plt.show()