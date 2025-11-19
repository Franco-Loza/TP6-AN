import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

def solver_richards_2d_circular(D_func, Lx, Ly, T_final, Nx, Ny, N):
    """
    Resuelve Richards 2D para una gota CIRCULAR.
    """
    dx, dy = Lx / (Nx + 1), Ly / (Ny + 1)
    dt = T_final / N
    x, y = np.linspace(0, Lx, Nx + 2), np.linspace(0, Ly, Ny + 2)
    X, Y = np.meshgrid(x, y, indexing='ij')

    # --- CONDICIÓN INICIAL CIRCULAR ---
    cx, cy = Lx / 2.0, Ly / 2.0
    R_gota = min(Lx, Ly) / 5.0

    # Inicialización
    theta = np.ones_like(X) * 1e-4
    mask = (X - cx) ** 2 + (Y - cy) ** 2 <= R_gota ** 2
    theta[mask] = 0.90

    # --- BUCLE TEMPORAL (ADI) ---
    for n in range(N):
        theta_k = theta.copy()
        for _ in range(15):
            theta_prev = theta_k.copy()
            D_vals = D_func(theta_k)

            # Promedios de D
            D_x = 0.5 * (D_vals[:-1, :] + D_vals[1:, :])
            D_y = 0.5 * (D_vals[:, :-1] + D_vals[:, 1:])

            rx, ry = dt / dx ** 2, dt / dy ** 2

            # --- Paso 1: Implícito en X ---
            theta_half = theta_k.copy()
            for j in range(1, Ny + 1):
                term_y = ry * (D_y[1:-1, j] * (theta_k[1:-1, j + 1] - theta_k[1:-1, j]) -
                               D_y[1:-1, j - 1] * (theta_k[1:-1, j] - theta_k[1:-1, j - 1]))

                # Construcción de diagonales X
                D_w = D_x[:-1, j]
                D_e = D_x[1:, j]

                main = 1.0 + rx * (D_w + D_e)
                upper = -rx * D_e[:-1]
                lower = -rx * D_w[1:]

                rhs = theta[1:-1, j] + term_y

                ab = np.zeros((3, Nx))
                ab[0, 1:] = upper
                ab[1, :] = main
                ab[2, :-1] = lower
                theta_half[1:-1, j] = solve_banded((1, 1), ab, rhs)

            # --- Paso 2: Implícito en Y ---
            theta_next = theta_half.copy()
            for i in range(1, Nx + 1):
                term_x = rx * (D_x[i, 1:-1] * (theta_half[i + 1, 1:-1] - theta_half[i, 1:-1]) -
                               D_x[i - 1, 1:-1] * (theta_half[i, 1:-1] - theta_half[i - 1, 1:-1]))

                # Construcción de diagonales Y
                D_s = D_y[i, :-1]
                D_n = D_y[i, 1:]

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

    return X, Y, theta, R_gota


def solver_1d_radial(D_func, R_max, T_final, Nr, N, R_gota):
    """Solver 1D en coordenadas cilíndricas para validar."""
    dr, dt = R_max / Nr, T_final / N
    r = np.linspace(0, R_max, Nr + 1)

    theta = np.ones(Nr + 1) * 1e-4
    theta[r <= R_gota] = 0.90

    for n in range(N):
        theta_k = theta.copy()
        for _ in range(10):  # Picard
            D_mid = 0.5 * (D_func(theta_k)[:-1] + D_func(theta_k)[1:])
            main, upper, lower, rhs = np.zeros(Nr + 1), np.zeros(Nr), np.zeros(Nr), np.zeros(Nr + 1)

            # Singularidad r=0
            alpha_0 = 4.0 * dt / dr ** 2 * D_mid[0]
            main[0], upper[0], rhs[0] = 1 + alpha_0, -alpha_0, theta[0]

            # Nodos internos
            for i in range(1, Nr):
                c_p = (dt / (r[i] * dr ** 2)) * (r[i] + dr / 2) * D_mid[i]
                c_m = (dt / (r[i] * dr ** 2)) * (r[i] - dr / 2) * D_mid[i - 1]
                main[i], upper[i], lower[i - 1], rhs[i] = 1 + c_p + c_m, -c_p, -c_m, theta[i]

            main[Nr], rhs[Nr] = 1.0, 1e-4
            ab = np.vstack([np.insert(upper, 0, 0), main, np.append(lower, 0)])
            theta_new = solve_banded((1, 1), ab, rhs)
            if np.linalg.norm(theta_new - theta_k) < 1e-5: break
            theta_k = theta_new
        theta = theta_k
    return r, theta


def graficar_resultados_D(X, Y, theta_2d, r_1d, theta_1d, Lx, Ly, T_final):
    """Genera el gráfico de validación recibiendo los datos calculados."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Mapa 2D
    c = ax1.contourf(X, Y, theta_2d, 20, cmap='viridis')
    plt.colorbar(c, ax=ax1, label='Humedad')
    ax1.set_title(f"Gota Circular 2D (t={T_final}s)")
    ax1.set_aspect('equal')

    center = Lx / 2
    R_dist = np.sqrt((X - center) ** 2 + (Y - center) ** 2).flatten()
    th_flat = theta_2d.flatten()
    mask = R_dist < Lx / 1.8

    ax2.plot(R_dist[mask], th_flat[mask], 'o', color='lightblue', alpha=0.3, label='Nodos 2D')
    ax2.plot(r_1d, theta_1d, 'r-', lw=2, label='Solución 1D Radial')
    ax2.set_title("Validación: Perfil Radial")
    ax2.set_xlabel("Radio (m)")
    ax2.set_ylabel("Saturación")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('actividadD_validacion.png')
    print("D: Gráfico 'actividadD_validacion.png' generado.")