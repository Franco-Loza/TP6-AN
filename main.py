# main.py - VERSIÓN CORREGIDA
import numpy as np
import matplotlib.pyplot as plt

from actividadA import resolucion_ecuacion_difusion_linea_1D, graficar_imagen_A
from actividadB import (resolucion_ecuacion_richards_1D_no_lineal, validacion_con_boltzmann,
                        graficar_comparacion_B, analizar_frente_humedad)
from actividadC import resolucion_ecuacion_difusion_2D, grafica_convergencia
from actividadDE import resolucion_ecuacion_richards_2D_no_lineal, extract_radial_profile, \
    grafico_comparacion_circular_eliptica
from models_soil_models import diffusivity_brooks_corey

if __name__ == "__main__":
    print("\n--- Ejecutando Actividad A: Difusión 1D Lineal (DF vs Analítica) ---")
    D0 = 0.01
    L = 1.0
    T_final = 0.5
    M = 49
    N = 1000
    x, theta_num, theta_an, error_L2, r = resolucion_ecuacion_difusion_linea_1D(D0, L, T_final, M, N)
    print(f"A: r={r:.4e}, error L2={error_L2:.4e}")
    graficar_imagen_A(x, theta_num, theta_an, T_final, r, error_L2, D0)

    print("\n--- Ejecutando Actividad B: Richards 1D No Lineal y validación con Boltzmann ---")
    # Parámetros para Richards 1D no lineal
    L_b = 0.5
    T_final_b = 0.1
    M_b = 100
    N_b = 200

    # Condición inicial - pulso gaussiano suave
    x_initial = np.linspace(0, L_b, M_b + 2)
    theta_initial = np.zeros(M_b + 2)
    center = L_b / 2
    sigma = 0.05
    theta_initial = 0.8 * np.exp(-((x_initial - center) ** 2) / (2 * sigma ** 2))

    print("B: Iniciando cálculo Richards 1D...")
    x_b, theta_b, cost_b = resolucion_ecuacion_richards_1D_no_lineal(
        diffusivity_brooks_corey, L_b, T_final_b, M_b, N_b, theta_initial
    )

    print("B: Richards 1D completado, iniciando Boltzmann...")
    # Validación CORREGIDA - solo pasar theta_initial
    x_boltzmann, theta_boltzmann, eta = validacion_con_boltzmann(
        diffusivity_brooks_corey, L_b, T_final_b, theta_initial
    )

    # Análisis del frente de humedad
    frente_pos, ancho_frente = analizar_frente_humedad(x_b, theta_b)
    print(f"B: Posición del frente = {frente_pos:.4f} m, Ancho del frente = {ancho_frente:.4f} m")
    print(f"B: Costo computacional = {cost_b}")

    # Graficar comparación
    graficar_comparacion_B(x_b, theta_b, x_boltzmann, theta_boltzmann, T_final_b, cost_b)
    print("B: Gráfico generado - actividadB_comparacion.png")

    print("\n--- Ejecutando Actividad C: Difusión 2D Lineal (ADI y Convergencia) ---")
    D0 = 0.01
    Lx = 1.0
    Ly = 1.0
    T_final = 0.1
    N = 100
    M_values = [19, 39, 79]
    errors_L2 = []
    for Mv in M_values:
        theta_num, theta_an, error = resolucion_ecuacion_difusion_2D(D0, Lx, Ly, T_final, Mv, Mv, N)
        errors_L2.append(error)
        print(f"C: M={Mv}, L2 error={error:.6e}")
    grafica_convergencia(M_values, errors_L2)

    print("\n--- Ejecutando Actividades D y E: Richards 2D No Lineal (Circular y Elíptica) ---")
    D_func = diffusivity_brooks_corey
    Lx = 1.0
    Ly = 1.0
    T_final = 0.05
    Nx = 50
    Ny = 50
    N = 200

    # Actividad D: Gota circular
    print("D: Resolviendo gota circular 2D...")
    X_d, Y_d, theta_d = resolucion_ecuacion_richards_2D_no_lineal(D_func, Lx, Ly, T_final, Nx, Ny, N)

    # Actividad E: Gota elíptica
    print("E: Resolviendo gota elíptica 2D...")
    Ly_e = 0.5 * Ly
    X_e, Y_e, theta_e = resolucion_ecuacion_richards_2D_no_lineal(D_func, Lx, Ly_e, T_final, Nx, Ny // 2, N)

    grafico_comparacion_circular_eliptica(X_d, Y_d, theta_d, X_e, Y_e, theta_e, Lx, Ly_e, T_final)

    # Validación radial - CORREGIDO (sin radio_gota)
    print("D: Realizando validación radial...")
    r_d, theta_radial_d = extract_radial_profile(X_d, Y_d, theta_d, Lx, Ly)
    L_1d = Lx / 2.0

    # Usar la función de la Actividad B para validación 1D - CONDICIÓN INICIAL SIMILAR
    x1d_initial = np.linspace(0, L_1d, M_b + 2)
    theta1d_initial = np.zeros(M_b + 2)

    # Crear condición inicial similar a la gota circular 2D
    radio_1d = 0.2  # Radio para la simulación 1D
    mask_1d = x1d_initial < radio_1d
    theta1d_initial[mask_1d] = 1.0  # Mismo valor máximo que en 2D

    x1d, theta1d, _ = resolucion_ecuacion_richards_1D_no_lineal(
        diffusivity_brooks_corey, L_1d, T_final, M_b, N_b, theta1d_initial
    )

    plt.figure(figsize=(10, 6))
    plt.plot(r_d, theta_radial_d, 'o', markersize=4, label='Perfil radial 2D (Circular)')
    plt.plot(x1d, theta1d, '-', linewidth=2, label='Perfil 1D (Validación)')
    plt.xlabel('Radio/Distancia (m)')
    plt.ylabel('Saturación $\\theta$')
    plt.title('Actividad D: Validación Radial - Comparación 2D vs 1D')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('actividadD_validacion_radial.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("D: Validación radial completada y gráfico guardado")

    print("\n✅ Ejecución principal finalizada exitosamente!")