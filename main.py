import numpy as np
from actividadA import resolucion_ecuacion_difusion_linea_1D, graficar_imagen_A
from actividadB import (resolucion_ecuacion_richards_1D_no_lineal, validacion_con_boltzmann, graficar_comparacion_B)
from actividadC import resolucion_ecuacion_difusion_2D, grafica_convergencia
from actividadD import solver_richards_2d_circular, solver_1d_radial, graficar_resultados_D
from actividadE import solver_richards_2d_eliptica, graficar_resultados_E
from actividadF import analisis_costo_computacional, graficar_escalamiento_computacional
from models_soil_models import diffusivity_brooks_corey

def separador(titulo):
    print("\n" + "=" * 60)
    print(f" {titulo}")
    print("=" * 60)

if __name__ == "__main__":
    # =========================================================
    # EJERCICIO A
    # =========================================================
    separador("EJERCICIO A: Difusión 1D Lineal")

    D0_a, L_a, T_a = 0.01, 1.0, 0.5
    M_a, N_a = 49, 1000

    x_a, th_num_a, th_an_a, err_a, r_a = resolucion_ecuacion_difusion_linea_1D(D0_a, L_a, T_a, M_a, N_a)
    print(f"A: Error L2 = {err_a:.4e}")
    graficar_imagen_A(x_a, th_num_a, th_an_a, T_a, r_a, err_a, D0_a)

    # =========================================================
    # EJERCICIO B
    # =========================================================
    separador("EJERCICIO B: Richards 1D No Lineal")

    L_b, T_b = 0.5, 0.1
    M_b, N_b = 100, 200

    # CI Gaussiana
    x_init = np.linspace(0, L_b, M_b + 2)
    th_init = 0.8 * np.exp(-((x_init - L_b / 2) ** 2) / (2 * 0.05 ** 2))

    print("B: Calculando numérica...")
    x_b, th_b, cost_b = resolucion_ecuacion_richards_1D_no_lineal(diffusivity_brooks_corey, L_b, T_b, M_b, N_b, th_init)

    print("B: Calculando Boltzmann...")
    x_boltz, th_boltz, eta = validacion_con_boltzmann(diffusivity_brooks_corey, L_b, T_b, th_init)

    graficar_comparacion_B(x_b, th_b, x_boltz, th_boltz, T_b, cost_b)

    # =========================================================
    # EJERCICIO C
    # =========================================================
    separador("EJERCICIO C: Difusión 2D Lineal")

    M_vals_c = [19, 39, 79]
    errs_c = []
    for m in M_vals_c:
        _, _, err = resolucion_ecuacion_difusion_2D(0.01, 1.0, 1.0, 0.1, m, m, 100)
        errs_c.append(err)
        print(f"C: Malla {m}x{m} -> Error {err:.2e}")
    grafica_convergencia(M_vals_c, errs_c)

    # =========================================================
    # EJERCICIO D
    # =========================================================
    separador("EJERCICIO D: Richards 2D Circular + Validación Radial")

    Lx, Ly, T_de = 0.4, 0.4, 10.0
    Nx, Ny, N_de = 50, 50, 100

    print("D: Resolviendo Richards 2D Circular...")
    X_d, Y_d, th_d, R_g = solver_richards_2d_circular(diffusivity_brooks_corey, Lx, Ly, T_de, Nx, Ny, N_de)

    print("D: Resolviendo Richards 1D Radial para validar...")
    r_d, th_radial_d = solver_1d_radial(diffusivity_brooks_corey, Lx / 1.5, T_de, 100, 100, R_g)

    graficar_resultados_D(X_d, Y_d, th_d, r_d, th_radial_d, Lx, Ly, T_de)

    # =========================================================
    # EJERCICIO E
    # =========================================================
    separador("EJERCICIO E: Richards 2D Elíptica")

    print("E: Resolviendo Richards 2D Elíptica...")
    X_e, Y_e, th_e = solver_richards_2d_eliptica(diffusivity_brooks_corey, Lx, Ly, T_de, 60, 60, N_de)

    graficar_resultados_E(X_e, Y_e, th_e, Lx, Ly, T_de)

    # =========================================================
    # EJERCICIO F
    # =========================================================
    separador("EJERCICIO F: Análisis de Costo y Discretización")

    res_f = analisis_costo_computacional()
    graficar_escalamiento_computacional(res_f)

    print("\n" + "=" * 85)
    print("RESUMEN FINAL: DISCRETIZACIÓN Y COSTO COMPUTACIONAL (Consigna F)")
    print("=" * 85)
    print(f"{'Actividad':<20} | {'Malla (Espacial)':<15} | {'Pasos (N)':<10} | {'Costo / Tiempo'}")
    print("-" * 85)

    # --- 1. Datos de Actividad A (Tomamos el caso más costoso del benchmark) ---
    idx_a = -1  # Último elemento
    m_a = res_f['A']['M'][idx_a]
    n_a = res_f['A']['N'][idx_a]
    t_a = res_f['A']['tiempos'][idx_a]
    print(f"{'A (1D Lineal)':<20} | {m_a:<15} | {n_a:<10} | {t_a:.4f} s (Benchmark)")

    # --- 2. Datos de Actividad B (Usamos las variables definidas arriba: M_b, N_b) ---
    # Costo teórico: O(N * M * Picard_iter). Asumimos ~5 iteraciones promedio.
    ops_b = N_b * M_b * 5
    print(f"{'B (1D No Lineal)':<20} | {M_b:<15} | {N_b:<10} | ~{ops_b:.1e} ops (Teórico)")

    # --- 3. Datos de Actividad C (Tomamos el caso más costoso del benchmark) ---
    idx_c = -1
    m_c = res_f['C']['M'][idx_c]
    # En actividadF.py, N está fijo en 100 para el benchmark de C
    n_c = 100
    t_c = res_f['C']['tiempos'][idx_c]
    print(f"{'C (2D Lineal ADI)':<20} | {m_c}x{m_c:<12} | {n_c:<10} | {t_c:.4f} s (Benchmark)")

    # --- 4. Datos de Actividad D (Usamos variables definidas: Nx, Ny, N_de) ---
    # Costo teórico: O(N * Nx * Ny * Picard). Asumimos ~10 iteraciones.
    ops_d = N_de * Nx * Ny * 10
    print(f"{'D (2D Circular)':<20} | {Nx}x{Ny:<12} | {N_de:<10} | ~{ops_d:.1e} ops (Teórico)")

    # --- 5. Datos de Actividad E (Usamos variables definidas) ---
    # En E usamos Nx=60, Ny=60 (definido en la llamada a la función)
    Nx_e, Ny_e = 60, 60
    ops_e = N_de * Nx_e * Ny_e * 10
    print(f"{'E (2D Elíptica)':<20} | {Nx_e}x{Ny_e:<12} | {N_de:<10} | ~{ops_e:.1e} ops (Teórico)")

    print("-" * 85)
    print("NOTA: 'ops' = Operaciones de punto flotante estimadas.")
    print("      Las actividades A y C reportan tiempo real de ejecución medido.")
    print("      Las actividades B, D y E reportan complejidad teórica (O) basada en la malla.")
    print("=" * 85)

    print("\nFIN DEL PROGRAMA.")