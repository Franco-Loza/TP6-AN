import numpy as np

# Parámetros Brooks-Corey (Tabla 1)
THETA_R = 2.378e-5    # saturación residual
THETA_S = 1.0         # saturación máxima (asumida 1)
D_SAT = 3.983e-6      # Ks/alpha aproximado [m^2/s]
N_BC = 4.795          # exponente Brooks-Corey

def _ensure_array(theta):
    """Interna: convierte a numpy array, guardando si la entrada era escalar."""
    was_scalar = False
    if np.isscalar(theta):
        theta = np.array([theta], dtype=np.float64)
        was_scalar = True
    else:
        theta = np.asarray(theta, dtype=np.float64)
    return theta, was_scalar

def diffusivity_brooks_corey(theta):
    """
    D(theta) según Brooks-Corey (forma empírica dada en la consigna).
    Acepta theta escalar o array. Devuelve escalar si la entrada fue escalar,
    o array si la entrada fue array.
    """
    theta_arr, was_scalar = _ensure_array(theta)

    delta_theta = THETA_S - THETA_R
    # Evitamos valores por debajo de theta_r
    theta_eff = np.maximum(theta_arr - THETA_R, 1e-12)
    Se = theta_eff / delta_theta
    Se = np.minimum(Se, 1.0)

    D_theta = D_SAT * (Se ** N_BC)

    if was_scalar:
        return float(D_theta[0])
    return D_theta

def dD_dtheta_brooks_corey(theta):
    """
    Derivada dD/dtheta analítica para Brooks-Corey.
    Acepta scalar o array; devuelve del mismo tipo (scalar o array).
    """
    theta_arr, was_scalar = _ensure_array(theta)

    delta_theta = THETA_S - THETA_R
    theta_eff = np.maximum(theta_arr - THETA_R, 1e-12)
    Se = theta_eff / delta_theta
    Se = np.minimum(Se, 1.0)

    # Para Se muy pequeñas, la derivada tiende a 0
    deriv = np.zeros_like(theta_arr, dtype=np.float64)
    mask = Se > 1e-12
    if np.any(mask):
        deriv[mask] = D_SAT * N_BC * (Se[mask] ** (N_BC - 1)) / delta_theta

    if was_scalar:
        return float(deriv[0])
    return deriv
