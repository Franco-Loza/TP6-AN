import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

def resolucion_ecuacion_richards_2D_no_lineal(D_func, Lx, Ly, T_final, Nx, Ny, N, picard_tol=1e-5, picard_maxiter=30):
    """
    ADI + Picard simplificado para 2D no lineal:
    - D_func: función que acepta array de theta y devuelve array de D (mismo shape)
    - Nx, Ny: puntos interiores (malla Nx+2 x Ny+2)
    """
    dx = Lx / (Nx + 1)
    dy = Ly / (Ny + 1)
    dt = T_final / N

    x = np.linspace(0, Lx, Nx + 2)
    y = np.linspace(0, Ly, Ny + 2)
    X, Y = np.meshgrid(x, y, indexing='ij')

    # condición inicial: gota circular en centro (valor 1 dentro radio R)
    theta = np.zeros_like(X)
    center_x = Lx / 2
    center_y = Ly / 2
    R = min(Lx, Ly) / 4
    mask = (X - center_x)**2 + (Y - center_y)**2 <= R**2
    theta[mask] = 1.0
    theta[~mask] = 0.0

    # interior indices 1..Nx, 1..Ny
    for n in range(N):
        theta_old = theta.copy()
        theta_k = theta_old.copy()  # Picard initial guess

        for it in range(picard_maxiter):
            theta_prev = theta_k.copy()

            # Evaluar D en nodos y puntos medios
            D_nodes = D_func(theta_k)
            # D at midpoints in x: shape (Nx+1, Ny+2) for midpoints between i and i+1
            D_x_mid = 0.5 * (D_nodes[1:, :] + D_nodes[:-1, :])   # shape (Nx+1, Ny+2)
            D_y_mid = 0.5 * (D_nodes[:, 1:] + D_nodes[:, :-1])   # shape (Nx+2, Ny+1)

            # construir matrices aproximadas por filas/columnas usando promedio en la fila/columna
            # paso 1: implícito en x (resolviendo por cada j)
            theta_star = theta_k.copy()
            for j in range(1, Ny + 1):
                # para cada fila j, construir tridiagonal de tamaño Nx
                # calcular D_{i+1/2, j} y D_{i-1/2, j} para i=1..Nx
                D_iphalf = D_x_mid[1:, j]   # length Nx
                D_imhalf = D_x_mid[:-1, j]  # length Nx

                a = - dt / (dx*dx) * D_imhalf
                c = - dt / (dx*dx) * D_iphalf
                b = 1.0 - (a + c)

                ab = np.zeros((3, Nx))
                ab[0, 1:] = c[:-1]
                ab[1, :] = b
                ab[2, :-1] = a[1:]

                rhs = theta_k[1:-1, j].copy()
                theta_star[1:-1, j] = solve_banded((1, 1), ab, rhs)

            # paso 2: implícito en y (resolviendo por cada i)
            theta_next = theta_star.copy()
            for i in range(1, Nx + 1):
                D_jphalf = D_y_mid[i, 1:]   # length Ny
                D_jmhalf = D_y_mid[i, :-1]  # length Ny

                a = - dt / (dy*dy) * D_jmhalf
                c = - dt / (dy*dy) * D_jphalf
                b = 1.0 - (a + c)

                ab = np.zeros((3, Ny))
                ab[0, 1:] = c[:-1]
                ab[1, :] = b
                ab[2, :-1] = a[1:]

                rhs = theta_star[i, 1:-1].copy()
                theta_next[i, 1:-1] = solve_banded((1, 1), ab, rhs)

            theta_k = theta_next.copy()
            # Bordes Dirichlet homogéneas
            theta_k[0,:] = 0.0
            theta_k[-1,:] = 0.0
            theta_k[:,0] = 0.0
            theta_k[:,-1] = 0.0

            relerr = np.linalg.norm(theta_k - theta_prev) / (np.linalg.norm(theta_prev) + 1e-12)
            if relerr < picard_tol:
                break

        theta = theta_k.copy()

    return X, Y, theta

def extract_radial_profile(X, Y, theta, Lx, Ly, num_points=100):
    center_x = Lx / 2.0
    center_y = Ly / 2.0
    R = np.sqrt((X - center_x)**2 + (Y - center_y)**2)
    r_max = min(Lx/2.0, Ly/2.0)
    bins = np.linspace(0, r_max, num_points)
    r_centers = 0.5*(bins[:-1] + bins[1:])
    theta_radial = np.zeros(num_points-1)
    for k in range(num_points-1):
        mask = (R >= bins[k]) & (R < bins[k+1])
        if np.any(mask):
            theta_radial[k] = np.mean(theta[mask])
        elif k>0:
            theta_radial[k] = theta_radial[k-1]
    return r_centers, theta_radial

def grafico_comparacion_circular_eliptica(X_d, Y_d, theta_d, X_e, Y_e, theta_e, Lx, Ly_e, T_final):
    """Genera comparación entre gota circular y elíptica"""
    plt.figure(figsize=(14, 6))

    # Gota circular
    plt.subplot(1, 2, 1)
    plt.contourf(X_d, Y_d, theta_d, levels=20, cmap='viridis')
    plt.colorbar(label='Saturación $\\theta$')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.title(f'Actividad D: Gota Circular - t={T_final:.3f}s')
    plt.axis('equal')

    # Gota elíptica
    plt.subplot(1, 2, 2)
    plt.contourf(X_e, Y_e, theta_e, levels=20, cmap='viridis')
    plt.colorbar(label='Saturación $\\theta$')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.title('Actividad E: Gota Elíptica 2:1')
    plt.axis('equal')

    plt.tight_layout()
    plt.savefig('actividadDE_gotas.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("DE: Gráfico de gotas circular vs elíptica generado")
