import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

def resolucion_ecuacion_difusion_linea_1D(D0, L, T_final, M, N):
    """
    Resuelve la ecuación de difusión lineal 1D (Actividad a)
    usando el esquema de Diferencias Finitas de Euler Implícito.

    Args:
        D0 (float): Coeficiente de difusividad constante.
        L (float): Longitud del dominio espacial.
        T_final (float): Tiempo final de la simulación.
        M (int): Número de puntos internos de la malla espacial (M+1 nodos).
        N (int): Número de pasos de tiempo.

    Returns:
        tuple: (x, theta_num, theta_an, error_L2, r)
            x (np.array): Coordenadas espaciales.
            theta_num (np.array): Solución numérica al tiempo T_final.
            theta_an (np.array): Solución analítica al tiempo T_final.
            error_L2 (float): Error L2 entre la solución numérica y analítica.
            r (float): El parámetro de estabilidad/precisión r = D0*dt/dx^2.
    """
    # 1. Parámetros de Discretización
    dx = L / (M + 1)
    dt = T_final / N
    r = D0 * dt / (dx ** 2)

    x = np.linspace(0, L, M + 2)
    theta = np.sin(np.pi * x / L)

    theta_interior = theta[1:-1]

    lambda_sq = (np.pi / L) ** 2
    theta_an = np.sin(np.pi * x / L) * np.exp(-D0 * lambda_sq * T_final)

    main_diag = (1 + 2 * r) * np.ones(M)
    off_diag = -r * np.ones(M - 1)

    A_banded = np.vstack([np.insert(off_diag, 0, 0),
                          main_diag,
                          np.append(off_diag, 0)])

    # 5. Bucle Temporal (Euler Implícito)
    for _ in range(N):
        # Lado derecho b = theta^n
        b = theta_interior.copy()

        # Resolvemos el sistema: A * theta^(n+1) = b
        theta_interior_new = solve_banded((1, 1), A_banded, b)

        theta_interior = theta_interior_new


    theta_num = np.zeros(M + 2)
    theta_num[1:-1] = theta_interior

    error_L2 = np.sqrt(np.sum((theta_num - theta_an) ** 2 * dx))

    return x, theta_num, theta_an, error_L2, r

def get_computational_cost_a(M, N):
    """
    Retorna el orden de complejidad computacional del DF 1D Implícito.
    """
    # El algoritmo de Thomas (resolución tridiagonal) es O(M).
    return f"O(N * M) = O({N} * {M})", N * M

def graficar_imagen_A(x, theta_num, theta_an, T_final, r, error_L2, D0):
    """Genera el gráfico comparando la solución numérica y analítica y lo guarda."""
    plt.figure(figsize=(10, 6))

    theta_init = np.sin(np.pi * x)
    plt.plot(x, theta_init, ':', color='gray', label='Condición Inicial ($t=0$)')

    plt.plot(x, theta_an, '-', color='red', linewidth=2, label=f'Solución Analítica ($t={T_final}$)')
    plt.plot(x, theta_num, 'o', color='blue', markersize=4, label=f'Solución Numérica ($t={T_final}$)')

    plt.xlabel('Posición $x$ (m)')
    plt.ylabel('Saturación $\\theta$')
    plt.title(f'Actividad a) - Difusión 1D Lineal ($D_0={D0}$)')

    text_info = (f'$r = D_0 \\Delta t/\\Delta x^2 = {r:.4f}$\n'
                 f'Error $L^2$: {error_L2:.6e}')
    plt.gca().text(0.05, 0.95, text_info, transform=plt.gca().transAxes,
                   verticalalignment='top', bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.7))

    plt.legend()
    plt.grid(True)

    plt.savefig('actividadA.png')
    plt.close()