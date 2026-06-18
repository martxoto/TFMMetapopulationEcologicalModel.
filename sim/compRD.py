import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === CONFIGURACIÓN DE ARCHIVOS ===
FILE_DIADAS = 'RDdiadas.csv'
FILE_TRIADAS = 'RDtriadas.csv'
OUTPUT_FIG = 'Comparativa_Robustez_Diadas_vs_Triadas.pdf'

# 1. Cargar los datos experimentales
df_diadas = pd.read_csv(FILE_DIADAS)
df_triadas = pd.read_csv(FILE_TRIADAS)

# Extraer el vector de dispersión D (es el eje X, idéntico para ambos)
D_vals = df_diadas['D']

# Aislar puramente los datos de robustez eliminando la columna D
trayectorias_diadas = df_diadas.drop(columns=['D'])
trayectorias_triadas = df_triadas.drop(columns=['D'])

# 2. Análisis Estadístico del Ensamble
# Pandas automáticamente ignora los valores NaN (vacíos) al calcular medias y desviaciones
media_diadas = trayectorias_diadas.mean(axis=1)
std_diadas = trayectorias_diadas.std(axis=1)

media_triadas = trayectorias_triadas.mean(axis=1)
std_triadas = trayectorias_triadas.std(axis=1)

# 3. Configuración del lienzo (Calidad Publicación/TFM)
plt.style.use('seaborn-v0_8-paper')
sns.set_context("paper", font_scale=1.5)
fig, ax = plt.subplots(figsize=(10, 6))

# 4. Trazar régimen de Díadas (Azul)
ax.plot(D_vals, media_diadas, color='#1f77b4', linewidth=3, marker='o', markersize=7, 
        label=f'Díadas (Media, N={trayectorias_diadas.shape[1]})')
ax.fill_between(D_vals, media_diadas - std_diadas, media_diadas + std_diadas, 
                color='#1f77b4', alpha=0.2)

# 5. Trazar régimen de Tríadas (Rojo)
ax.plot(D_vals, media_triadas, color='#d62728', linewidth=3, marker='s', markersize=7, 
        label=f'Tríadas (Media, N={trayectorias_triadas.shape[1]})')
ax.fill_between(D_vals, media_triadas - std_triadas, media_triadas + std_triadas, 
                color='#d62728', alpha=0.2)

# 6. Formato Físico y Matemático
ax.set_xscale('symlog', linthresh=0.01) # symlog permite graficar el 0.0 sin romper el logaritmo
ax.set_xlabel('Tasa de Dispersión (D)', fontsize=14, fontweight='bold')
ax.set_ylabel('Robustez Media del Ensamble (R)', fontsize=14, fontweight='bold')

# Forzar marcadores legibles en el eje X
ticks_importantes = [0.0, 0.01, 0.1, 1.0, 10.0, 100.0]
ax.set_xticks(ticks_importantes)
ax.set_xticklabels([str(t) for t in ticks_importantes])

# Elementos estéticos auxiliares
ax.legend(loc='lower right', frameon=True, shadow=True, fontsize=12)
ax.grid(True, which="major", ls="-", alpha=0.4)
ax.grid(True, which="minor", ls="--", alpha=0.1)

# 7. Exportación
plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=300, format='pdf')
print(f"[+] Gráfica comparativa exportada exitosamente a: {OUTPUT_FIG}")
plt.show()