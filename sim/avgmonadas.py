import os
import shutil
import subprocess
import numpy as np
import pandas as pd

# --- 1. CONFIGURACIÓN DE PARÁMETROS ---
rp, rv, ap, av, D= -0.1, -0.1, 8.0, 8.0, 1.0

# --- 2. RUTAS SEGÚN TU ESTRUCTURA DE CARPETAS ---
ejecutable = "./bin/expbackw" 
carpeta_redes = "interactions"
archivo_output_cpp = "robustnessR.txt" 

# La ruta donde C++ espera encontrar el archivo copiado
ruta_destino_cpp = os.path.join(carpeta_redes, "redActual.txt")

# --- 3. BUSCAR ARCHIVOS ---
archivos_monadas = [f for f in os.listdir(carpeta_redes) if f.startswith('3') and f.endswith('.txt')]

if not archivos_monadas:
    print(f"Error: No se encontraron archivos que empiecen por '1' en la carpeta '{carpeta_redes}'.")
    exit()

print(f"Se encontraron {len(archivos_monadas)} mónadas para analizar.\n")

# --- 4. BUCLE DE SIMULACIÓN ---
resultados = []

for red in archivos_monadas:
    print(f"Simulando {red}...", end=" ", flush=True)
    
    # SEGURO: Borramos el output viejo antes de empezar
    if os.path.exists(archivo_output_cpp):
        os.remove(archivo_output_cpp)
    
    # Copiamos la red directamente DENTRO de la carpeta interactions
    ruta_origen = os.path.join(carpeta_redes, red)
    shutil.copy(ruta_origen, ruta_destino_cpp)
    
    # Ejecutamos el simulador
    comando = [ejecutable, str(rp), str(rv), str(ap), str(av), str(D)]
    subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Leemos el resultado
    try:
        with open(archivo_output_cpp, "r") as f:
            rint = float(f.read().strip())
            resultados.append({'Red': red, 'Robustez_Rint': rint})
            print(f"Rint = {rint:.4f}")
    except FileNotFoundError:
        print("❌ Error: C++ falló y no generó resultados.")

# Limpiamos el archivo falso de la carpeta interactions
if os.path.exists(ruta_destino_cpp):
    os.remove(ruta_destino_cpp)

# --- 5. ESTADÍSTICA FINAL ---
if resultados:
    df = pd.DataFrame(resultados)
    media = np.mean(df['Robustez_Rint'])
    desviacion = np.std(df['Robustez_Rint'], ddof=1)
    
    print("\n" + "="*40)
    print("      RESULTADOS GLOBALES (MÓNADAS)")
    print("="*40)
    print(f"Redes analizadas : {len(df)}")
    print(f"Media (Rint)     : {media:.4f}")
    print(f"Desv. Típica     : {desviacion:.4f}")
    print("="*40)
    
    df.to_csv("resultados_monadas.csv", index=False)