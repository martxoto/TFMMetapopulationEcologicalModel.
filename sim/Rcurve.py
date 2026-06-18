import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === CONFIGURACIÓN ===
EXECUTABLE = './bin/expbacklocal'  
METRIC_FILE = "robustnessR.txt" 

# Parámetros que enviaremos a C++
rp, rv, ap, av = -0.5, -0.8, 5.0, 9.0
D_values = [0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 5.0, 10.0, 20.0, 50.0, 100.0]
#D_values = [0.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0]  
results_R = []

print(f"Iniciando barrido de D para {len(D_values)} valores...")
print("-" * 40)
print(f"{'D':<10} | {'R (Robustez)':<15}")
print("-" * 40)

for D in D_values:
    # 1. Le pasamos los 5 parámetros en fila india
    comando = [EXECUTABLE, str(rp), str(rv), str(ap), str(av), str(D)]
    
    if os.path.exists(METRIC_FILE):
        os.remove(METRIC_FILE)
        
    try:
        # 2. Lanzamos el comando silenciando la salida de C++ para no ensuciar la terminal
        subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        
        # 3. Leer el resultado
        if os.path.exists(METRIC_FILE):
            with open(METRIC_FILE, 'r') as f:
                r_val = float(f.read().strip())
                results_R.append(r_val)
                print(f"{D:<10} | {r_val:.5f}")
        else:
            print(f"Error: No se generó {METRIC_FILE} para D={D}")
            results_R.append(0)
            
    except subprocess.CalledProcessError as e:
        print(f"Error ejecutando C++ para D={D}: {e}")
        results_R.append(0)

# === GRAFICAR ===
sns.set_theme(style="whitegrid")
plt.figure(figsize=(10, 6))

plt.plot(D_values, results_R, 'o-', color='navy', linewidth=2, markersize=8)

max_R = max(results_R)
best_D = D_values[results_R.index(max_R)]
plt.plot(best_D, max_R, 'ro', markersize=12, label=f'Óptimo (D={best_D}, R={max_R:.3f})')

plt.title("Sensibilidad a la Dispersión (Target Attack en Díada)", fontsize=14, fontweight='bold')
plt.xlabel("Tasa de Dispersión ($D$)", fontsize=12)
plt.ylabel("Robustez del Ecosistema ($R_{int}$)", fontsize=12)
plt.xscale('symlog', linthresh=0.1) 
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()

plt.tight_layout()
plt.savefig('optimizacion_D.png', dpi=300)
print("-" * 40)
print(f"Gráfica guardada como 'optimizacion_D.png'. Mejor D = {best_D}")
plt.show()