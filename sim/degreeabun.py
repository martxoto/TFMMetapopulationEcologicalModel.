import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# === CONFIGURACIÓN ===
INTERACTION_FILE = 'interactions/interactions_Penhale_Sands_patches.txt' 
EVO_P_FILE = 'evolutionp.txt'
EVO_V_FILE = 'evolutionv.txt'

sns.set_theme(style="whitegrid")

def parse_interactions_and_get_degrees(filename):
    """
    Replica la lógica de carga de C++ para asegurar que los IDs coinciden
    y calcula el Grado (número de parejas únicas) de cada especie.
    """
    plant_index = {}
    insect_index = {}
    plant_count = 0
    insect_count = 0
    
    # Estructuras para guardar conexiones (Sets para contar únicos)
    plant_connections = {} # ID -> Set of insect IDs
    insect_connections = {} # ID -> Set of plant IDs

    print(f"Leyendo red desde {filename}...")
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) < 4: continue
                
                # Formato: site patch plant insect weight (o similar según tu parser C++)
                # Tu parser C++: site >> plant >> insect >> weight
                # Python split() separa por espacios
                # Asumimos que el archivo tiene 4 columnas: site_id plant_name insect_name weight
                
                # Ajusta esto si tu archivo tiene otro formato
                # En tu C++ loadGamma: site >> plant >> insect >> weight
                try:
                    site = parts[0]
                    plant = parts[1]
                    insect = parts[2]
                    weight = float(parts[3])
                except ValueError:
                    continue 

                # Asignar IDs (Lógica idéntica a C++)
                if plant not in plant_index:
                    plant_index[plant] = plant_count
                    plant_connections[plant_count] = set()
                    plant_count += 1
                
                if insect not in insect_index:
                    insect_index[insect] = insect_count
                    insect_connections[insect_count] = set()
                    insect_count += 1
                
                p_id = plant_index[plant]
                v_id = insect_index[insect]
                
                # Guardar conexión si el peso es positivo
                if weight > 0:
                    plant_connections[p_id].add(v_id)
                    insect_connections[v_id].add(p_id)
                    
    except FileNotFoundError:
        print(f"Error: No encuentro {filename}")
        return None, None

    # Calcular grados
    k_plants = [len(plant_connections[i]) for i in range(plant_count)]
    k_insects = [len(insect_connections[i]) for i in range(insect_count)]
    
    return k_plants, k_insects

def get_abundances(filename, num_species):
    """Lee la última línea del archivo de evolución para obtener el estado estacionario"""
    try:
        # Leemos solo la última línea para ahorrar memoria
        with open(filename, 'r') as f:
            lines = f.readlines()
            if not lines: return None
            last_line = lines[-1].strip().split()
            
            # La primera columna es el tiempo, el resto son poblaciones
            # Formato: t p0_s0 p0_s1 ... p1_s0 ...
            # Necesitamos sumar la población de cada especie en todos sus parches
            
            data = np.array([float(x) for x in last_line[1:]])
            
            # El array data tiene longitud: num_species * num_patches
            # Deducimos num_patches
            if len(data) % num_species != 0:
                print(f"Error dimensiones: Datos {len(data)} no divisible por especies {num_species}")
                return None
            
            num_patches = int(len(data) / num_species)
            
            # Reshape a (num_species, num_patches) y sumar parches
            abundances_matrix = data.reshape((num_species, num_patches))
            total_abundances = np.sum(abundances_matrix, axis=1)
            
            return total_abundances
            
    except FileNotFoundError:
        print(f"Error: No encuentro {filename}")
        return None

# === EJECUCIÓN ===

# 1. Obtener Estructura (Grados)
k_p, k_v = parse_interactions_and_get_degrees(INTERACTION_FILE)

if k_p is not None:
    # 2. Obtener Dinámica (Abundancias)
    n_p = get_abundances(EVO_P_FILE, len(k_p))
    n_v = get_abundances(EVO_V_FILE, len(k_v))

    if n_p is not None and n_v is not None:
        
        # === GRÁFICA 1: CORRELACIÓN GRADO-ABUNDANCIA ===
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Plantas
        axes[0].scatter(k_p, n_p, color='forestgreen', alpha=0.7, s=80, edgecolors='k')
        axes[0].set_title("Plantas: Grado vs Abundancia")
        axes[0].set_xlabel("Grado (Número de polinizadores)")
        axes[0].set_ylabel("Abundancia de Equilibrio (log)")
        axes[0].set_yscale('log') # Escala logarítmica suele verse mejor
        axes[0].set_xscale('log')
        axes[0].grid(True, which="both", ls="--", alpha=0.3)
        
        # Correlación
        if len(k_p) > 1:
            m, b = np.polyfit(np.log1p(k_p), np.log1p(n_p), 1)
            axes[0].plot(k_p, np.exp(b) * np.power(k_p, m), 'r--', alpha=0.5, label=f'Tendencia (pendiente={m:.2f})')
            axes[0].legend()

        # Insectos
        axes[1].scatter(k_v, n_v, color='darkorange', alpha=0.7, s=80, edgecolors='k')
        axes[1].set_title("Insectos: Grado vs Abundancia")
        axes[1].set_xlabel("Grado (Número de plantas)")
        axes[1].set_ylabel("Abundancia de Equilibrio (log)")
        axes[1].set_yscale('log')
        axes[1].set_xscale('log')
        axes[1].grid(True, which="both", ls="--", alpha=0.3)

        plt.tight_layout()
        plt.savefig('correlacion_grado_abundancia.png', dpi=300)
        plt.show()

        # === GRÁFICA 2: DISTRIBUCIÓN DE ABUNDANCIAS (Curva de Rango-Abundancia) ===
        plt.figure(figsize=(10, 6))
        
        # Ordenamos de mayor a menor
        sorted_np = np.sort(n_p)[::-1]
        sorted_nv = np.sort(n_v)[::-1]
        
        plt.plot(range(1, len(sorted_np)+1), sorted_np, 'o-', color='forestgreen', label='Plantas')
        plt.plot(range(1, len(sorted_nv)+1), sorted_nv, 'o-', color='darkorange', label='Insectos')
        
        plt.yscale('log')
        plt.title("Curvas de Rango-Abundancia (Whittaker plots)")
        plt.xlabel("Rango de Especie (1 = Más abundante)")
        plt.ylabel("Abundancia (log)")
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.3)
        
        plt.savefig('distribucion_abundancias.png', dpi=300)
        plt.show()