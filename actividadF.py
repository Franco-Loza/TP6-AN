import numpy as np
import matplotlib.pyplot as plt
import time
from actividadA import resolucion_ecuacion_difusion_linea_1D
from actividadC import resolucion_ecuacion_difusion_2D

def analisis_costo_computacional():
    """Retorna un diccionario con los datos de tiempos."""
    resultados = {}
    print(">>> Midiendo tiempos (Actividad F)...")

    # A: 1D Lineal
    M_values_A = [49, 99, 199, 399]
    N_values_A = [500, 1000, 2000, 4000]
    tiempos_A = []
    for M, N in zip(M_values_A, N_values_A):
        start = time.time()
        resolucion_ecuacion_difusion_linea_1D(0.01, 1.0, 0.5, M, N)
        tiempos_A.append(time.time() - start)
    resultados['A'] = {'M': M_values_A, 'N': N_values_A, 'tiempos': tiempos_A}

    # C: 2D Lineal
    M_values_C = [9, 19, 39, 79]
    tiempos_C = []
    for M in M_values_C:
        start = time.time()
        resolucion_ecuacion_difusion_2D(0.01, 1.0, 1.0, 0.1, M, M, 100)
        tiempos_C.append(time.time() - start)
    resultados['C'] = {'M': M_values_C, 'tiempos': tiempos_C}

    return resultados


def graficar_escalamiento_computacional(resultados):
    plt.figure(figsize=(10, 6))

    # A
    M_A, N_A = resultados['A']['M'], resultados['A']['N']
    ops_A = [m * n for m, n in zip(M_A, N_A)]
    plt.loglog(ops_A, resultados['A']['tiempos'], 'o-', label='Actividad A (1D Lineal)')

    # C
    M_C = resultados['C']['M']
    ops_C = [100 * m ** 2 for m in M_C]
    plt.loglog(ops_C, resultados['C']['tiempos'], 's-', label='Actividad C (2D Lineal ADI)')

    # Referencias
    x_ref = np.logspace(4, 7, 50)
    plt.loglog(x_ref, x_ref * (resultados['A']['tiempos'][-1] / ops_A[-1]), 'k--', alpha=0.5, label='O(N)')

    plt.xlabel('Operaciones')
    plt.ylabel('Tiempo [s]')
    plt.title('Actividad F: Escalamiento')
    plt.legend()
    plt.grid(True, which="both", alpha=0.2)
    plt.savefig('actividadF_escalamiento.png')
    print("F: Gráfico 'actividadF_escalamiento.png' generado.")