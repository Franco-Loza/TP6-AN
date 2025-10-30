import numpy as np
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt
from models_soil_models import diffusivity_brooks_corey
from boltzmann_edo import solve_boltzmann_edo


def resolucion_ecuacion_richards_1D_no_lineal(D_func, L, T_final, M, N, theta_initial=None, picard_tol=1e-6,
                                              picard_maxiter=20):
    """
    Resuelve la ecuación de Richards 1D no lineal usando Euler Implícito + Picard.
    """
    # Discretización
    dx = L / (M + 1)
    dt = T_final / N

    print(f"B: Discretización - dx={dx:.4f}, dt={dt:.6f}, M={M}, N={N}")

    # Dominio espacial
    x = np.linspace(0, L, M + 2)

    # Condición inicial
    if theta_initial is None:
        theta = np.zeros(M + 2)
        center_idx = len(theta) // 2
        theta[max(0, center_idx - 3):min(M + 2, center_idx + 3)] = 0.8
    else:
        theta = theta_initial.copy()

    # Bucle temporal
    for n in range(N):
        if n % 100 == 0:
            print(f"B: Paso temporal {n}/{N}")

        theta_old = theta.copy()
        theta_k = theta_old.copy()

        for picard_iter in range(picard_maxiter):
            theta_prev = theta_k.copy()

            # Evaluar D en los puntos medios
            D_nodes = D_func(theta_k)
            D_half = np.zeros(M + 1)
            for i in range(M + 1):
                D_half[i] = 0.5 * (D_nodes[i] + D_nodes[i + 1])

            # Construir matriz del sistema
            main_diag = np.ones(M)
            lower_diag = np.zeros(M - 1)
            upper_diag = np.zeros(M - 1)
            rhs = np.zeros(M)

            for i in range(M):
                idx = i + 1
                a = -dt / (dx ** 2) * D_half[i]
                c = -dt / (dx ** 2) * D_half[i + 1]
                b = 1.0 - (a + c)

                main_diag[i] = b
                if i > 0:
                    lower_diag[i - 1] = a
                if i < M - 1:
                    upper_diag[i] = c

                rhs[i] = theta_old[idx]

                if i == 0:
                    rhs[i] -= a * theta_k[0]
                if i == M - 1:
                    rhs[i] -= c * theta_k[-1]

            # Resolver sistema tridiagonal
            A_banded = np.vstack([np.insert(upper_diag, 0, 0),
                                  main_diag,
                                  np.append(lower_diag, 0)])

            try:
                theta_interior = solve_banded((1, 1), A_banded, rhs)
                theta_k[1:M + 1] = theta_interior
            except Exception as e:
                print(f"Error resolviendo sistema: {e}")
                break

            diff = np.linalg.norm(theta_k - theta_prev)
            if diff < picard_tol:
                break

        theta = theta_k.copy()

    computational_cost = f"O(N * M * picard_iter) ≈ O({N} * {M})"
    return x, theta, computational_cost


def validacion_con_boltzmann(D_func, L, T_final, theta_initial):
    """
    Valida la solución numérica con la transformación de Boltzmann CORREGIDA.
    """
    try:
        print("B: Iniciando validación con Boltzmann...")

        # Encontrar el valor máximo de saturación para usar como theta(0)
        theta_max = np.max(theta_initial)
        print(f"B: Valor máximo de saturación: {theta_max:.3f}")

        # Para Boltzmann: θ(0) = theta_max, θ(∞) = 0
        theta_a = theta_max
        theta_b = 0.0

        # Resolver EDO de Boltzmann
        eta, theta_eta = solve_boltzmann_edo(theta_a, theta_b)

        # Transformación CORREGIDA: x = η * √(4*D0*t)
        # Usamos D0 promedio para la transformación
        D0_promedio = D_func(theta_max / 2)
        x_boltzmann = eta * np.sqrt(4 * D0_promedio * T_final)
        theta_boltzmann = theta_eta

        print(f"B: Boltzmann completado - {len(x_boltzmann)} puntos, D0_promedio={D0_promedio:.2e}")
        return x_boltzmann, theta_boltzmann, eta

    except Exception as e:
        print(f"Error en validación Boltzmann: {e}")
        print("  Continuando sin validación Boltzmann...")
        return None, None, None


def graficar_comparacion_B(x_numeric, theta_numeric, x_boltzmann, theta_boltzmann, T_final, computational_cost):
    """Genera gráfico comparando solución numérica y Boltzmann."""
    plt.figure(figsize=(12, 8))

    # Solución numérica
    plt.plot(x_numeric, theta_numeric, 'b-', linewidth=2, label='Solución Numérica (Richards 1D)')

    # Solución Boltzmann (si está disponible)
    if x_boltzmann is not None and theta_boltzmann is not None:
        # Solo graficar donde x_boltzmann está dentro del dominio
        mask = x_boltzmann <= np.max(x_numeric)
        plt.plot(x_boltzmann[mask], theta_boltzmann[mask], 'ro', markersize=4,
                 label='Solución Boltzmann', alpha=0.7)

    plt.xlabel('Posición x (m)')
    plt.ylabel('Saturación $\\theta$')
    plt.title(f'Actividad B: Richards 1D No Lineal - t={T_final:.3f}s')
    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.figtext(0.02, 0.02, f'Costo: {computational_cost}',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    plt.tight_layout()
    plt.savefig('actividadB_comparacion.png', dpi=300, bbox_inches='tight')
    plt.close()


def analizar_frente_humedad(x, theta, threshold=0.01):
    """Analiza la posición del frente de humedad."""
    mask = theta > threshold
    if np.any(mask):
        frente_pos = x[mask][-1]
        theta_max = np.max(theta)
        mask_10 = theta > 0.1 * theta_max
        mask_90 = theta > 0.9 * theta_max
        if np.any(mask_10) and np.any(mask_90):
            ancho_frente = x[mask_10][-1] - x[mask_90][0] if len(x[mask_10]) > 0 and len(x[mask_90]) > 0 else 0.0
        else:
            ancho_frente = 0.0
    else:
        frente_pos = 0.0
        ancho_frente = 0.0

    return frente_pos, ancho_frente