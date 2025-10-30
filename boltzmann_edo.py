import numpy as np
from scipy.integrate import solve_ivp
from models_soil_models import diffusivity_brooks_corey, dD_dtheta_brooks_corey


def boltzmann_ode(eta, y):
    """Ecuación de Boltzmann para difusión no lineal."""
    theta, dtheta = y

    # Asegurar que son escalares
    theta = float(theta)
    dtheta = float(dtheta)

    # Limitar valores físicamente razonables
    theta = np.clip(theta, 1e-4, 0.99)
    dtheta = np.clip(dtheta, -10.0, -1e-6)  # Solo valores negativos pequeños

    D = diffusivity_brooks_corey(theta)

    # Evitar división por cero
    if D < 1e-12:
        d2theta = -eta * dtheta / 2.0
    else:
        dD = dD_dtheta_brooks_corey(theta)
        d2theta = -((eta / 2) * dtheta + dD * (dtheta ** 2)) / D

    return [dtheta, d2theta]


def solve_boltzmann_edo(theta_a, theta_b, eta_max=5.0, n_points=200):
    """
    Método de Boltzmann MEJORADO con manejo robusto.
    """
    print(f"Boltzmann: Resolviendo EDO - θ(0)={theta_a:.3f}, θ(∞)={theta_b}")

    # Usar diferentes pendientes iniciales para encontrar la correcta
    pendientes_prueba = [-0.001, -0.01, -0.1, -0.5, -1.0]
    mejor_sol = None
    mejor_error = float('inf')

    for s0 in pendientes_prueba:
        try:
            eta_eval = np.linspace(0, eta_max, n_points)
            sol = solve_ivp(
                boltzmann_ode,
                [0, eta_max],
                [theta_a, s0],
                t_eval=eta_eval,
                method='RK45',
                rtol=1e-4,
                atol=1e-6,
                max_step=0.1
            )

            if sol.success:
                theta_final = sol.y[0][-1]
                error = abs(theta_final - theta_b)

                if error < mejor_error:
                    mejor_error = error
                    mejor_sol = (sol.t, sol.y[0])
                    print(f"  s0={s0:.3f}, θ(η_max)={theta_final:.3f}, error={error:.3f}")

        except Exception as e:
            print(f"  Error con s0={s0:.3f}: {e}")
            continue

    if mejor_sol is not None:
        eta, theta_eta = mejor_sol
        print(f"  Mejor solución: s0={pendientes_prueba[pendientes_prueba.index(s0)]:.3f}, error={mejor_error:.3f}")
        return eta, theta_eta
    else:
        # Fallback: solución exponencial
        print("  Usando solución exponencial de fallback")
        eta_fallback = np.linspace(0, eta_max, n_points)
        theta_fallback = theta_a * np.exp(-eta_fallback)
        return eta_fallback, theta_fallback