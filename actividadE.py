import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

def solver_richards_2d_eliptica(D_func, Lx, Ly, T_final, Nx, Ny, N):
    """
    Resuelve Richards 2D para una gota ELÍPTICA (Relación 2:1).
    """
    dx, dy = Lx / (Nx + 1), Ly / (Ny + 1)
    dt = T_final / N
    x, y = np.linspace(0, Lx, Nx + 2), np.linspace(0, Ly, Ny + 2)
    X, Y = np.meshgrid(x, y, indexing='ij')

    # --- CONDICIÓN INICIAL ELÍPTICA ---
    cx, cy = Lx / 2.0, Ly / 2.0
    R_base = min(Lx, Ly) / 6.0
    a, b = R_base * 2.0, R_base

    theta = np.ones_like(X) * 1e-4
    mask = ((X - cx) ** 2 / a ** 2) + ((Y - cy) ** 2 / b ** 2) <= 1.0
    theta[mask] = 0.90

    for n in range(N):
        theta_k = theta.copy()
        for _ in range(15):
            theta_prev = theta_k.copy()
            D_vals = D_func(theta_k)

            # Promedios de D
            D_x = 0.5 * (D_vals[:-1, :] + D_vals[1:, :])
            D_y = 0.5 * (D_vals[:, :-1] + D_vals[:, 1:])

            rx, ry = dt / dx ** 2, dt / dy ** 2

            # --- Paso 1: Implícito X ---
            theta_half = theta_k.copy()
            for j in range(1, Ny + 1):
                term_y = ry * (D_y[1:-1, j] * (theta_k[1:-1, j + 1] - theta_k[1:-1, j]) -
                               D_y[1:-1, j - 1] * (theta_k[1:-1, j] - theta_k[1:-1, j - 1]))

                D_w = D_x[:-1, j]  # i-1/2
                D_e = D_x[1:, j]  # i+1/2

                main = 1.0 + rx * (D_w + D_e)
                upper = -rx * D_e[:-1]
                lower = -rx * D_w[1:]

                rhs = theta[1:-1, j] + term_y

                ab = np.zeros((3, Nx))
                ab[0, 1:] = upper
                ab[1, :] = main
                ab[2, :-1] = lower
                theta_half[1:-1, j] = solve_banded((1, 1), ab, rhs)

            # --- Paso 2: Implícito Y ---
            theta_next = theta_half.copy()
            for i in range(1, Nx + 1):
                term_x = rx * (D_x[i, 1:-1] * (theta_half[i + 1, 1:-1] - theta_half[i, 1:-1]) -
                               D_x[i - 1, 1:-1] * (theta_half[i, 1:-1] - theta_half[i - 1, 1:-1]))

                D_s = D_y[i, :-1]  # j-1/2
                D_n = D_y[i, 1:]  # j+1/2

                main = 1.0 + ry * (D_s + D_n)
                upper = -ry * D_n[:-1]
                lower = -ry * D_s[1:]

                rhs = theta[i, 1:-1] + term_x

                ab = np.zeros((3, Ny))
                ab[0, 1:] = upper
                ab[1, :] = main
                ab[2, :-1] = lower
                theta_next[i, 1:-1] = solve_banded((1, 1), ab, rhs)

            theta_next[[0, -1], :] = 1e-4
            theta_next[:, [0, -1]] = 1e-4
            theta_k = theta_next
            if np.linalg.norm(theta_k - theta_prev) < 1e-4: break
        theta = theta_k

    return X, Y, theta

def graficar_resultados_E(X, Y, theta, Lx, Ly, T_final):
    """Grafica la gota elíptica."""
    plt.figure(figsize=(7, 6))
    c = plt.contourf(X, Y, theta, 20, cmap='Oranges')
    plt.colorbar(c, label='Humedad')
    plt.title(f"Actividad E: Gota Elíptica 2:1 (t={T_final}s)")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.axis('equal')

    th = np.linspace(0, 2 * np.pi, 100)
    cx, cy = Lx / 2, Ly / 2
    R_base = min(Lx, Ly) / 6.0
    plt.plot(cx + 2 * R_base * np.cos(th), cy + R_base * np.sin(th), 'k--', alpha=0.5, label='Borde Inicial')
    plt.legend()

    plt.tight_layout()
    plt.savefig('actividadE_eliptica.png')
    print("E: Gráfico 'actividadE_eliptica.png' generado.")