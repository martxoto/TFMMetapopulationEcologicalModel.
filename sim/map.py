import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import pandas as pd

# 1. DEFINIR EL BARRIDO DE PARÁMETROS (Mortalidad rp y rv)
# Rango suave/medio para darles oportunidad de sobrevivir al inicio
puntos = 15
rp_vals = np.linspace(-0.1, -1, puntos)
rv_vals = np.linspace(-0.1, -1, puntos)

# Matrices para los dos mapas
matriz_ext_plantas = np.zeros((len(rp_vals), len(rv_vals)))
matriz_plantas_ini = np.zeros((len(rp_vals), len(rv_vals)))

print("Iniciando barrido de rp y rv...")
total_sims = puntos * puntos
sim_actual = 0

for i, rp in enumerate(rp_vals):
    for j, rv in enumerate(rv_vals):
        sim_actual += 1
        
        # EJECUTAR C++ (Asegúrate de que tu main recibe rp y rv como argv[1] y argv[2])
        comando = ["./exp1map", str(rp), str(rv)]
        subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        # LEER RESULTADOS
        try:
            df = pd.read_csv("results.txt", sep='\s+', comment='#', header=None)
            fila_k0 = df[df[0] == 0]
            fila_k1 = df[df[0] == 1]
            
            if not fila_k0.empty and not fila_k1.empty:
                plantas_ini = fila_k0.iloc[0][2]
                plantas_k1 = fila_k1.iloc[0][2]
                
                ext_sec_plantas = (plantas_ini - 1) - plantas_k1
                
                # Guardamos ambos datos
                matriz_plantas_ini[i, j] = plantas_ini
                matriz_ext_plantas[i, j] = max(0, ext_sec_plantas)
            else:
                matriz_plantas_ini[i, j] = 0
                matriz_ext_plantas[i, j] = -1
                
        except Exception as e:
            matriz_plantas_ini[i, j] = 0
            matriz_ext_plantas[i, j] = -1
            
        if sim_actual % 10 == 0:
            print(f"Progreso: {sim_actual}/{total_sims} simulaciones...")

# 3. GRAFICAR LOS DOS MAPAS
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

x_labels = np.round(rv_vals, 2)
y_labels = np.round(rp_vals, 2)

# MAPA 1: Biodiversidad Inicial
sns.heatmap(matriz_plantas_ini, xticklabels=x_labels, yticklabels=y_labels, 
            cmap="viridis", annot=True, ax=ax1)
ax1.set_title("Biodiversidad Inicial (Plantas vivas en k=0)")
ax1.set_xlabel("Tasa de Insectos (rv)")
ax1.set_ylabel("Tasa de Plantas (rp)")
ax1.invert_yaxis()

# MAPA 2: Cascadas de Extinción
sns.heatmap(matriz_ext_plantas, xticklabels=x_labels, yticklabels=y_labels, 
            cmap="YlOrRd", annot=True, ax=ax2)
ax2.set_title("Extinciones Secundarias de Plantas (Tras 1 corte)")
ax2.set_xlabel("Tasa de Insectos (rv)")
ax2.set_ylabel("Tasa de Plantas (rp)")
ax2.invert_yaxis()

plt.tight_layout()
plt.savefig("mapa_doble_rp_rv.png", dpi=300)
plt.show()