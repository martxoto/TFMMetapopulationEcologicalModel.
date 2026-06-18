import numpy as np
import pandas as pd
import subprocess
import csv
import os
from scipy.stats import qmc

# 1. DEFINIR LÍMITES DE EXPLORACIÓN PARA LOS 4 PARÁMETROS
# Formato: [Mínimo, Máximo]
limites = {
    'rp': [-1.0, -0.01],     # Mortalidad plantas (suave a media)
    'rv': [-1.0, -0.01],     # Mortalidad insectos
    'alphap': [1.0, 10.0],   # Beneficio para plantas
    'alphav': [1.0, 10.0]    # Beneficio para insectos
}

num_simulaciones = 150  # Puedes subirlo a 5000 si lo dejas por la noche

print(f"Generando {num_simulaciones} combinaciones con LHS...")
sampler = qmc.LatinHypercube(d=4)
sample = sampler.random(n=num_simulaciones)

# Escalar las muestras a nuestros límites reales
l_bounds = [limites[k][0] for k in limites]
u_bounds = [limites[k][1] for k in limites]
parametros_lhs = qmc.scale(sample, l_bounds, u_bounds)

# Preparar archivo CSV de salida (Añadimos Total Plantas y Total Insectos)
archivo_salida = "parametros_optimos.csv"
with open(archivo_salida, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['rp', 'rv', 'alphap', 'alphav', 'Plantas_Ini', 'Insectos_Ini', 'Total_Ext_Sec_Plantas', 'Total_Ext_Sec_Insectos'])

print("Iniciando criba de parámetros...")

exitos = 0

for i, row in enumerate(parametros_lhs):
    rp, rv, ap, av = row
    
    # 2. EJECUTAR EL SIMULADOR DE RAREZA
    # Asegúrate de poner el nombre correcto del ejecutable que lee argumentos
    comando = ["./bin/expbackw", str(rp), str(rv), str(ap), str(av), str(5.0)]
    subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # 3. EVALUAR EL PATRÓN (Pattern-Oriented Modeling)
    try:
        # Leemos el experimento de rareza (results.txt)
        df = pd.read_csv("resultsTarget.txt", sep=r'\s+', comment='#', header=None)
        
        # Fila inicial (k=0)
        fila_k0 = df[df[0] == 0]
        if fila_k0.empty: continue
        
        plantas_ini = fila_k0.iloc[0][2]
        insectos_ini = fila_k0.iloc[0][3]
        
        # CONDICIÓN 1: Supervivencia inicial digna (>= 2 para el toy model o redes pequeñas)
        if plantas_ini >= 6 and insectos_ini >= 2:
            
            total_ext_sec_plantas = 0
            total_ext_sec_insectos = 0
            
            # Recorremos los pasos para ver si hubo extinciones secundarias
            for step in range(1, len(df)):
                plantas_eliminadas_manual = df.iloc[step][0]
                surv_plants = df.iloc[step][2]
                surv_insects = df.iloc[step][3]
                
                # CÁLCULO DEL TOTAL ACUMULADO
                # Plantas: (Las que había al inicio - las que he matado yo) - las que quedan vivas
                ext_p = (plantas_ini - plantas_eliminadas_manual) - surv_plants
                
                # Insectos: Como yo no mato insectos a mano, todas sus muertes son secundarias
                ext_i = insectos_ini - surv_insects
                
                # Actualizamos el total (el máximo valor acumulado de la resta es el total de caídas)
                if ext_p > total_ext_sec_plantas:
                    total_ext_sec_plantas = ext_p
                if ext_i > total_ext_sec_insectos:
                    total_ext_sec_insectos = ext_i
            
            # CONDICIÓN 2: Tienen que haber cascadas secundarias de plantas
            if total_ext_sec_plantas > 0:
                exitos += 1
                print(f"🎯 ¡ÉXITO {exitos}! rp={rp:.3f}, rv={rv:.3f}, ap={ap:.2f}, av={av:.2f} | Plantas Ini:{plantas_ini} -> Ext. Plantas:{total_ext_sec_plantas} | Ext. Insectos:{total_ext_sec_insectos}")
                
                # Guardamos en el CSV inmediatamente por si se corta la ejecución
                with open(archivo_salida, mode='a', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow([rp, rv, ap, av, plantas_ini, insectos_ini, total_ext_sec_plantas, total_ext_sec_insectos])

    except Exception as e:
        pass # Ignoramos ejecuciones fallidas o archivos vacíos

    # Barra de progreso
    if (i + 1) % 100 == 0:
        print(f"Progreso: {i + 1}/{num_simulaciones} completados. Éxitos hasta ahora: {exitos}")

print(f"\n¡Búsqueda terminada! Se encontraron {exitos} configuraciones óptimas.")
print(f"Revisa el archivo '{archivo_salida}'.")