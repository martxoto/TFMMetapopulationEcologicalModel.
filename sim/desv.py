import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === CONFIGURACIÓN ===
FILE_DIADAS = 'RDdiadas.csv'
FILE_TRIADAS = 'RDtriadas.csv'

# 1. Cargar datos
df_diadas = pd.read_csv(FILE_DIADAS)
df_triadas = pd.read_csv(FILE_TRIADAS)

D_vals = df_diadas['D']

# 2. Calcular exclusivamente las Desviaciones Estándar (sigma)
std_diadas = df_diadas.drop(columns=['D']).std(axis=1)
std_triadas = df_triadas.drop(columns=['D']).std(axis=1)

# 3. Configurar el lienzo
plt.style.use('seaborn-v0_8-paper')
sns.set_context("paper", font_scale=1.5)
fig, ax = plt.subplots(figsize=(8, 5))

# 4. Trazar la comparativa de dispersión
ax.plot(D_vals, std_diadas, color='#1f77b4', linewidth=3, marker='o', label='Variabilidad Díadas ($\sigma$)')
ax.plot(D_vals, std_triadas, color='#d62728', linewidth=3, marker='s', label='Variabilidad Tríadas ($\sigma$)')

# 5. Formato
ax.set_xscale('symlog', linthresh=0.01)
ax.set_xlabel('Tasa de Dispersión ($D$)', fontweight='bold')
ax.set_ylabel('Desviación Estándar de la Robustez ($\sigma$)', fontweight='bold')

# Rellenar el área entre las dos curvas para destacar la "ganancia" en estabilidad
ax.fill_between(D_vals, std_triadas, std_diadas, where=(std_diadas > std_triadas), 
                interpolate=True, color='green', alpha=0.1, label='Reducción de Varianza')

ax.legend(loc='upper right', frameon=True)
ax.grid(True, which="major", ls="-", alpha=0.3)

plt.tight_layout()
plt.savefig('Analisis_Varianza_Topologica.pdf', dpi=300)
plt.show()