import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

# === CONFIGURACIÓN ===
# Nombre de tus archivos de salida C++
FILES_TO_PLOT = ['evolutionp.txt', 'evolutionv.txt'] 

# ¿Quieres escala logarítmica? (Útil para ver si se van a cero o explotan)
USE_LOG_SCALE = False 

def plot_file(filename):
    if not os.path.exists(filename):
        print(f"AVISO: El archivo '{filename}' no existe. Saltando...")
        return

    print(f"Procesando {filename}...")

    try:
        # Leemos el archivo. 
        # sep='\s+' maneja espacios múltiples generados por C++
        # header=None porque tu C++ no escribe cabecera "t p1 p2..."
        df = pd.read_csv(filename, sep='\s+', header=None, engine='python')
        
        # Eliminamos columnas vacías (a veces C++ deja un espacio al final de la línea)
        df = df.dropna(axis=1, how='all')

        # La columna 0 es el TIEMPO (t)
        t = df.iloc[:, 0]
        
        # El resto de columnas son las POBLACIONES
        # (Si tienes parches, cada columna es una especie en un parche)
        populations = df.iloc[:, 1:]

        # === GRAFICAR ===
        plt.figure(figsize=(10, 6))
        
        # Pintamos todas las columnas
        # El bucle automático asigna colores distintos
        for col in populations.columns:
            plt.plot(t, populations[col], label=f'Columna {col}')

        plt.title(f'Evolución Temporal: {filename}')
        plt.xlabel('Tiempo (t)')
        plt.ylabel('Biomasa / Población')
        
        if USE_LOG_SCALE:
            plt.yscale('log')
            plt.ylabel('Biomasa (Escala Log)')

        # Si hay demasiadas leyendas (muchas especies), ocultamos la leyenda para que no tape
        if len(populations.columns) <= 10:
            plt.legend(loc='best')
        else:
            print(f"  -> Hay {len(populations.columns)} series. Ocultando leyenda para claridad.")

        plt.grid(True, which='both', linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        # Guardar imagen
        img_name = filename.replace('.txt', '.png')
        plt.savefig(img_name)
        print(f"  -> Gráfica guardada como: {img_name}")
        
        # Mostrar en pantalla
        plt.show()

    except Exception as e:
        print(f"ERROR leyendo {filename}: {e}")
        print("Asegúrate de que el archivo no esté vacío y tenga formato numérico.")

# === MAIN ===
if __name__ == "__main__":
    # Si le pasas un archivo por argumento, usa ese. Si no, usa la lista por defecto.
    if len(sys.argv) > 1:
        for f in sys.argv[1:]:
            plot_file(f)
    else:
        for f in FILES_TO_PLOT:
            plot_file(f)