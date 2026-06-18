import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 1. CONFIGURACIÓN DE ESTILO (Para que quede de artículo científico)
sns.set_theme(style="whitegrid")
plt.figure(figsize=(10, 6))

# 2. DEFINIR ARCHIVOS Y SUS PROPIEDADES VISUALES
archivos = {
    "results.txt": {"label": "Rareza (Ataca plantas raras)", "color": "forestgreen", "linestyle": "-"},
    "resultsRandom.txt": {"label": "Aleatorio (Random)", "color": "royalblue", "linestyle": "--"},
    "resultsTarget.txt": {"label": "Dirigido (Ataca Hubs)", "color": "firebrick", "linestyle": "-."}
}

# 3. LEER Y PINTAR CADA ARCHIVO
for archivo, props in archivos.items():
    if os.path.exists(archivo):
        try:
            # Leemos el archivo. 
            # sep='\s+' separa por espacios. comment='#' ignora la cabecera.
            df = pd.read_csv(archivo, sep=r'\s+', comment='#', header=None)
            
            # Columna 0: Número de plantas eliminadas manualmente (kEffective)
            # Columna 1: Ratio de Robustez
            eje_x = df[0]
            eje_y = df[1]
            
            plt.plot(eje_x, eje_y, 
                     label=props["label"], 
                     color=props["color"], 
                     linestyle=props["linestyle"], 
                     linewidth=2.5)
            print(f"✅ Cargado y procesado: {archivo}")
        except Exception as e:
            print(f"❌ Error al procesar {archivo}: {e}")
    else:
        print(f"⚠️ OJO: No se encontró el archivo {archivo}. ¿Lo has ejecutado?")

# 4. DETALLES ESTÉTICOS DEL GRÁFICO
plt.title("Tolerancia del Ecosistema ante Estrategias de Extinción", fontsize=15, pad=15)
plt.xlabel("Nº de Plantas Eliminadas Directamente", fontsize=12)
plt.ylabel("Robustez (Fracción de especies supervivientes)", fontsize=12)

# Fijamos los límites del eje Y de 0 a 1 (100% de especies iniciales)
plt.ylim(-0.05, 1.05)
plt.xlim(left=0)

plt.legend(fontsize=11, frameon=True, shadow=True)
plt.tight_layout()

# 5. GUARDAR Y MOSTRAR
plt.savefig("comparativa_robustez.png", dpi=300)
print("Gráfica guardada como 'comparativa_robustez.png'")
plt.show()