import subprocess
import pandas as pd
import os

# === CONFIGURACIÓN ===
EXECUTABLE = './bin/expbacklocal'  
METRIC_FILE = "robustnessR.txt" 
OUTPUT_CSV = "datos_barrido_D.csv"

# Parámetros fijos (Escenario Canónico)
rp, rv, ap, av = -0.5, -0.8, 5.0, 9.0
D_values = [0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 5.0, 10.0, 20.0, 50.0, 100.0]

results_R = []

print(f"Iniciando simulación para {len(D_values)} valores de D...")

# === EJECUCIÓN Y RECOLECCIÓN ===
for D in D_values:
    comando = [EXECUTABLE, str(rp), str(rv), str(ap), str(av), str(D)]
    
    if os.path.exists(METRIC_FILE):
        os.remove(METRIC_FILE)
        
    try:
        subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        
        if os.path.exists(METRIC_FILE):
            with open(METRIC_FILE, 'r') as f:
                r_val = float(f.read().strip())
                results_R.append(r_val)
                print(f"D={D:<6} | R={r_val:.5f}")
        else:   
            print(f"Error: No se generó {METRIC_FILE} para D={D}")
            results_R.append(float('nan'))
            
    except subprocess.CalledProcessError as e:
        print(f"Error ejecutando C++ para D={D}")
        results_R.append(float('nan'))

# === GESTIÓN DEL ARCHIVO CSV (AÑADIR COLUMNAS) ===
if os.path.exists(OUTPUT_CSV):
    # Leer archivo existente
    df = pd.read_csv(OUTPUT_CSV)
    
    # Control de errores: evitar mezclar vectores de D de distintos tamaños
    if len(df) != len(D_values):
        print("¡ERROR FATAL! El número de filas no coincide con el archivo existente.")
        exit(1)
        
    # Calcular el nombre de la nueva columna (Ej: si hay D y Run_1, el nuevo es Run_2)
    numero_ejecucion = len(df.columns) # La primera columna es D, así que la longitud equivale al nuevo índice
    nombre_columna = f"Run_{numero_ejecucion}"
    
    # Añadir los datos y guardar
    df[nombre_columna] = results_R
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\n[+] Datos añadidos a '{OUTPUT_CSV}' en la columna: {nombre_columna}")

else:
    # Crear archivo nuevo si no existe
    nombre_columna = "Run_1"
    df = pd.DataFrame({
        'D': D_values,
        nombre_columna: results_R
    })
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\n[+] Archivo '{OUTPUT_CSV}' creado con la columna: {nombre_columna}")