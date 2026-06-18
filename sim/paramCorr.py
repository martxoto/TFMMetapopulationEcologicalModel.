import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Cargar los datos
try:
    df = pd.read_csv("parametros_optimos.csv")
    print(f"Cargadas {len(df)} configuraciones exitosas.")
except FileNotFoundError:
    print("Error: No se encuentra 'parametros_optimos.csv'.")
    exit()

# 2. Configuración de estilo
sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Título global
#fig.suptitle("Análisis del Espacio de Parámetros Viables (Pattern-Oriented Modeling)", fontsize=16, y=1.05)

# --- GRÁFICO 1: Compensación rp vs alphap ---
# Añadimos legend=False
sns.scatterplot(data=df, x='rp', y='alphap', 
                hue='Plantas_Ini', size='Total_Ext_Sec_Plantas', 
                sizes=(50, 300), palette='viridis', alpha=0.8, ax=axes[0], legend=False)
#axes[0].set_title("1. Compensación: Mortalidad vs Beneficio", fontsize=12)
axes[0].set_xlabel("rp", fontsize=11)
axes[0].set_ylabel("alphap", fontsize=11)

# --- GRÁFICO 2: La Frontera de Supervivencia (rp vs rv) ---
# Añadimos legend=False
sns.scatterplot(data=df, x='rp', y='rv', 
                hue='Plantas_Ini', size='Total_Ext_Sec_Plantas', 
                sizes=(50, 300), palette='viridis', alpha=0.8, ax=axes[1], legend=False)
#axes[1].set_title("2. Trade-off de Mortalidad Cruzada", fontsize=12)
axes[1].set_xlabel("rp", fontsize=11)
axes[1].set_ylabel("rv", fontsize=11)

# --- GRÁFICO 3: Asimetría del Mutualismo (alphav vs alphap) ---
# Aquí SÍ dejamos que se genere la leyenda
sns.scatterplot(data=df, x='alphav', y='alphap', 
                hue='Plantas_Ini', size='Total_Ext_Sec_Plantas', 
                sizes=(50, 300), palette='viridis', alpha=0.8, ax=axes[2])
axes[2].plot([1, 10], [1, 10], color='red', linestyle='--', label='alphap = alphav')
#axes[2].set_title("3. Asimetría del Mutualismo", fontsize=12)
axes[2].set_xlabel("alphav", fontsize=11)
axes[2].set_ylabel("alphap", fontsize=11)

# EL TRUCO: Movemos la leyenda fuera del gráfico a la derecha
axes[2].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

# Ajustes estéticos finales
plt.tight_layout()
# bbox_inches='tight' asegura que la leyenda no se recorte al guardar la imagen
plt.savefig("analisis_parametros_lhs.png", dpi=300, bbox_inches='tight')
print("¡Gráfico guardado como 'analisis_parametros_lhs.png'!")
plt.show()