import numpy as np
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt
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
    
    NOTA IMPORTANTE: Para difusividad fuertemente no lineal (Brooks-Corey con n=4.795),
    la transformación de Boltzmann proporciona solo una aproximación cualitativa.
    La solución exacta requeriría métodos numéricos más sofisticados.
    """
    try:
        print("B: Iniciando validación con Boltzmann...")

        # Encontrar el valor máximo de saturación para usar como theta(0)
        theta_max = np.max(theta_initial)
        print(f"B: Valor máximo de saturación: {theta_max:.3f}")

        # Para Boltzmann: θ(0) = theta_max, θ(∞) = 0
        theta_a = theta_max
        theta_b = 0.0

        # Para difusividad fuertemente no lineal, usar un η_max moderado
        # que capture la región de transición principal
        D0_promedio = D_func(theta_max / 2)
        
        # Elegir η_max para capturar ~80% del dominio espacial
        # x ≈ 0.8*L => η ≈ 0.8*L / √(4*D0*t)
        eta_max_practico = 0.8 * L / np.sqrt(4 * D0_promedio * T_final)
        # Limitar a un rango manejable computacionalmente
        eta_max = min(500.0, max(50.0, eta_max_practico))
        
        print(f"B: Usando eta_max = {eta_max:.1f} (D0={D0_promedio:.2e} m^2/s)")

        # Resolver EDO de Boltzmann con parámetros robustos
        eta, theta_eta = solve_boltzmann_edo(theta_a, theta_b, eta_max=eta_max, n_points=300)

        # Transformación: x = η * √(4*D0*t)
        x_boltzmann = eta * np.sqrt(4 * D0_promedio * T_final)
        theta_boltzmann = theta_eta

        print(f"B: Boltzmann completado - {len(x_boltzmann)} puntos")
        print(f"B: Rango de x_boltzmann: [{x_boltzmann[0]:.4f}, {x_boltzmann[-1]:.4f}] m")
        print(f"B: Rango de theta_boltzmann: [{np.min(theta_boltzmann):.4f}, {np.max(theta_boltzmann):.4f}]")
        print(f"B: NOTA: Solucion de Boltzmann es aproximada para D(theta) fuertemente no lineal")
        
        return x_boltzmann, theta_boltzmann, eta

    except Exception as e:
        print(f"Error en validación Boltzmann: {e}")
        import traceback
        traceback.print_exc()
        print("  Continuando sin validación Boltzmann...")
        return None, None, None


def graficar_comparacion_B(x_numeric, theta_numeric, x_boltzmann, theta_boltzmann, T_final, computational_cost):
    """Genera gráfico comparando solución numérica y Boltzmann con subplot adicional."""
    fig = plt.figure(figsize=(14, 10))
    
    # Subplot 1: Comparación completa
    ax1 = plt.subplot(2, 1, 1)
    ax1.plot(x_numeric, theta_numeric, 'b-', linewidth=2, label='Solución Numérica (Richards 1D)')
    
    if x_boltzmann is not None and theta_boltzmann is not None:
        mask = x_boltzmann <= np.max(x_numeric)
        ax1.plot(x_boltzmann[mask], theta_boltzmann[mask], 'ro', markersize=3,
                 label=f'Solución Boltzmann (cubre {x_boltzmann[-1]:.3f}m)', alpha=0.7)
        
        # Agregar región sombreada donde Boltzmann es válido
        ax1.axvspan(0, x_boltzmann[-1], alpha=0.1, color='red', 
                    label=f'Región Boltzmann ({x_boltzmann[-1]/np.max(x_numeric)*100:.1f}% del dominio)')
    
    ax1.set_xlabel('Posición x (m)', fontsize=11)
    ax1.set_ylabel('Saturación θ', fontsize=11)
    ax1.set_title(f'Actividad B: Richards 1D No Lineal - t={T_final:.3f}s (Vista Completa)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9)
    
    # Subplot 2: Zoom en la región de Boltzmann
    ax2 = plt.subplot(2, 1, 2)
    
    if x_boltzmann is not None and theta_boltzmann is not None and len(x_boltzmann) > 0:
        # Zoom en la región donde Boltzmann tiene datos
        x_zoom_max = min(x_boltzmann[-1] * 1.5, 0.15)  # 50% más allá o 15cm max
        mask_numeric = x_numeric <= x_zoom_max
        mask_boltz = x_boltzmann <= x_zoom_max
        
        ax2.plot(x_numeric[mask_numeric], theta_numeric[mask_numeric], 'b-', 
                linewidth=2, label='Solución Numérica (Richards 1D)')
        ax2.plot(x_boltzmann[mask_boltz], theta_boltzmann[mask_boltz], 'ro-', 
                markersize=4, linewidth=1, label='Solución Boltzmann', alpha=0.7)
        
        ax2.set_xlabel('Posición x (m)', fontsize=11)
        ax2.set_ylabel('Saturación θ', fontsize=11)
        ax2.set_title(f'Zoom: Región de Validez de Boltzmann (0 - {x_zoom_max:.3f}m)', fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=9)
        
        # Nota explicativa
        nota = (f"NOTA: La solución de Boltzmann solo cubre {x_boltzmann[-1]:.3f}m "
                f"({x_boltzmann[-1]/np.max(x_numeric)*100:.1f}% del dominio)\n"
                f"debido a la fuerte no linealidad de D(θ) con n=4.795.\n"
                f"Para cubrir L={np.max(x_numeric):.1f}m se requeriría η ≈ {np.max(x_numeric)/x_boltzmann[-1]*500:.0f} (numéricamente impracticable).")
        ax2.text(0.02, 0.05, nota, transform=ax2.transAxes, fontsize=8,
                verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    else:
        ax2.text(0.5, 0.5, 'Solución de Boltzmann no disponible', 
                ha='center', va='center', fontsize=12, transform=ax2.transAxes)
    
    # Costo computacional
    fig.text(0.02, 0.98, f'Costo: {computational_cost}',
            fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle="round,pad=0.4", facecolor="lightblue", alpha=0.8))
    
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig('actividadB_comparacion.png', dpi=300, bbox_inches='tight')
    plt.close()