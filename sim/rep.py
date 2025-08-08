import pandas as pd
import matplotlib.pyplot as plt

def graficar_columnas_txt(fichero, columnas_y):
    """
    Grafica las columnas dadas contra la primera columna (que se asume es el tiempo).
    
    Args:
    - fichero: str, ruta al archivo .txt.
    - columnas_y: lista de índices (desde 1) o nombres de columnas a graficar.
    """
    # Leer el archivo (asumiendo separación por espacios o tabs)
    df = pd.read_csv(fichero, sep='\s+', engine='python')

    print(f"Columnas detectadas: {list(df.columns)}")

    # Convertir los índices a nombres si hace falta
    nombres_y = []
    for col in columnas_y:
        if isinstance(col, int):
            idx = col - 1  # Convertir a base 0
            if idx < 1 or idx >= df.shape[1]:
                raise ValueError(f"Número de columna inválido: {col}")
            nombres_y.append(df.columns[idx])
        elif isinstance(col, str):
            if col not in df.columns:
                raise ValueError(f"Nombre de columna '{col}' no encontrado.")
            nombres_y.append(col)

    # Graficar
    plt.figure(figsize=(10, 6))
    for name in nombres_y:
        plt.plot(df.iloc[:, 0], df[name], label=name)
    plt.xlabel("t")
    plt.ylabel("Abundances")
    plt.legend(loc='upper left')
    ax = plt.gca()
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(False)
    plt.tight_layout()
    plt.show()

# Ejemplo de uso
graficar_columnas_txt(
    'evolutionp.txt',
    list(range(2, 38))  # columnas 2 a 37
)
