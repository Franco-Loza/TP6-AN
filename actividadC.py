# actividadC.py
import numpy as np
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt

def analytical_solution_2d(x, y, t, D0, Lx, Ly):
    """
    Solución analítica para el caso 2D con condición inicial sin(pi x/Lx)*sin(pi y/Ly)
    """
    X, Y = np.meshgrid(x, y, indexing='ij')
    lambda_x = np.pi / Lx
    lambda_y = np.pi / Ly
    return np.sin(lambda_x * X) * np.sin(lambda_y * Y) * np.exp(-D0 * (lambda_x**2 + lambda_y**2) * t)

def resolucion_ecuacion_difusion_2D(D0, Lx, Ly, T_final, Nx, Ny, N):
    """
    Método ADI (Alternating Direction Implicit) para la ecuación de difusión 2D:
        ∂θ/∂t = D0 * (∂²θ/∂x² + ∂²θ/∂y²)
    Condiciones de borde: Dirichlet homogéneas (θ=0 en el contorno)
    """
    dx = Lx / (Nx + 1)
    dy = Ly / (Ny + 1)
    dt = T_final / N

    x = np.linspace(0, Lx, Nx + 2)
    y = np.linspace(0, Ly, Ny + 2)

    # Condición inicial
    X, Y = np.meshgrid(x, y, indexing='ij')
    theta = np.sin(np.pi * X / Lx) * np.sin(np.pi * Y / Ly)

    # Coeficientes del método ADI
    rx = D0 * dt / (2 * dx * dx)
    ry = D0 * dt / (2 * dy * dy)

    # Matrices tridiagonales constantes (banded form para solve_banded)
    Ax_ab = np.zeros((3, Nx))
    Ax_ab[1, :] = 1 + 2 * rx
    Ax_ab[0, 1:] = -rx
    Ax_ab[2, :-1] = -rx

    Ay_ab = np.zeros((3, Ny))
    Ay_ab[1, :] = 1 + 2 * ry
    Ay_ab[0, 1:] = -ry
    Ay_ab[2, :-1] = -ry

    # Bucle temporal
    theta_int = theta.copy()
    for _ in range(N):
        # Paso 1: implícito en x
        theta_star = np.zeros_like(theta_int)
        for j in range(1, Ny + 1):
            b = np.zeros(Nx)
            for i in range(1, Nx + 1):
                b[i - 1] = (1 - 2 * ry) * theta_int[i, j] + ry * (theta_int[i, j - 1] + theta_int[i, j + 1])
            theta_star[1:-1, j] = solve_banded((1, 1), Ax_ab, b)

        # Paso 2: implícito en y
        theta_next = np.zeros_like(theta_int)
        for i in range(1, Nx + 1):
            b = np.zeros(Ny)
            for j in range(1, Ny + 1):
                b[j - 1] = (1 - 2 * rx) * theta_star[i, j] + rx * (theta_star[i - 1, j] + theta_star[i + 1, j])
            theta_next[i, 1:-1] = solve_banded((1, 1), Ay_ab, b)

        # Aplicar condiciones de borde
        theta_next[0, :] = theta_next[-1, :] = 0.0
        theta_next[:, 0] = theta_next[:, -1] = 0.0
        theta_int = theta_next.copy()

    # Solución analítica y error
    theta_an = analytical_solution_2d(x, y, T_final, D0, Lx, Ly)
    dxdy = dx * dy
    error_L2 = np.sqrt(np.sum((theta_int - theta_an) ** 2) * dxdy)
    return theta_int, theta_an, error_L2

def grafica_convergencia(M_values, errors_L2):
    """
    Genera una gráfica log-log del error L2 en función del tamaño de malla M
    para estimar el orden de convergencia del método ADI.
    """
    plt.figure(figsize=(7,5))
    plt.loglog(M_values, errors_L2, 'o-', linewidth=2, markersize=6, label='Error L2')
    plt.xlabel('Número de puntos interiores (M)')
    plt.ylabel('Error L2')
    plt.title('Convergencia método ADI 2D')
    plt.grid(True, which='both', linestyle='--', linewidth=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig('actividadC_convergencia.png', dpi=200)
    plt.close()
