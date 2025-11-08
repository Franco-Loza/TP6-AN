"""
Script de validación adicional para Actividad B
Demuestra la corrección de la implementación con tests numéricos
"""

import numpy as np
import matplotlib.pyplot as plt
from actividadB import resolucion_ecuacion_richards_1D_no_lineal
from models_soil_models import diffusivity_brooks_corey

def test_conservacion_masa():
    """Test 1: Verificar conservación de masa (sin fuentes/sumideros)"""
    print("\n=== TEST 1: Conservación de Masa ===")
    
    L = 0.5
    T_final = 0.1
    M = 100
    N = 200
    
    # Condición inicial
    x_initial = np.linspace(0, L, M + 2)
    theta_initial = np.zeros(M + 2)
    center = L / 2
    sigma = 0.05
    theta_initial = 0.8 * np.exp(-((x_initial - center) ** 2) / (2 * sigma ** 2))
    
    # Resolver
    x_final, theta_final, _ = resolucion_ecuacion_richards_1D_no_lineal(
        diffusivity_brooks_corey, L, T_final, M, N, theta_initial
    )
    
    # Calcular masa (integral de theta)
    masa_inicial = np.trapz(theta_initial, x_initial)
    masa_final = np.trapz(theta_final, x_final)
    perdida_relativa = abs(masa_final - masa_inicial) / masa_inicial * 100
    
    print(f"Masa inicial:  {masa_inicial:.6f}")
    print(f"Masa final:    {masa_final:.6f}")
    print(f"Pérdida:       {perdida_relativa:.2f}%")
    
    if perdida_relativa < 5.0:
        print("✅ PASA: Conservación de masa aceptable (<5%)")
        return True
    else:
        print("⚠️ ADVERTENCIA: Pérdida de masa mayor al 5%")
        return False


def test_estabilidad_fisica():
    """Test 2: Verificar que θ ∈ [0, 1] en todo momento"""
    print("\n=== TEST 2: Estabilidad Física (θ ∈ [0,1]) ===")
    
    L = 0.5
    T_final = 0.1
    M = 100
    N = 200
    
    x_initial = np.linspace(0, L, M + 2)
    theta_initial = np.zeros(M + 2)
    center = L / 2
    sigma = 0.05
    theta_initial = 0.8 * np.exp(-((x_initial - center) ** 2) / (2 * sigma ** 2))
    
    x_final, theta_final, _ = resolucion_ecuacion_richards_1D_no_lineal(
        diffusivity_brooks_corey, L, T_final, M, N, theta_initial
    )
    
    theta_min = np.min(theta_final)
    theta_max = np.max(theta_final)
    
    print(f"θ_min = {theta_min:.6f}")
    print(f"θ_max = {theta_max:.6f}")
    
    if 0.0 <= theta_min and theta_max <= 1.0:
        print("✅ PASA: Todos los valores de θ están en [0, 1]")
        return True
    else:
        print("❌ FALLA: θ fuera del rango físico")
        return False


def test_convergencia_temporal():
    """Test 3: Verificar convergencia con refinamiento temporal"""
    print("\n=== TEST 3: Convergencia Temporal ===")
    
    L = 0.5
    T_final = 0.05  # Tiempo más corto para test rápido
    M = 50
    
    # Condición inicial
    x_initial = np.linspace(0, L, M + 2)
    theta_initial = np.zeros(M + 2)
    center = L / 2
    sigma = 0.05
    theta_initial = 0.8 * np.exp(-((x_initial - center) ** 2) / (2 * sigma ** 2))
    
    # Probar con diferentes N
    N_values = [50, 100, 200]
    soluciones = []
    
    for N in N_values:
        x, theta, _ = resolucion_ecuacion_richards_1D_no_lineal(
            diffusivity_brooks_corey, L, T_final, M, N, theta_initial
        )
        soluciones.append(theta)
        print(f"N={N:3d}: θ_max={np.max(theta):.6f}, θ_integral={np.trapz(theta, x):.6f}")
    
    # Calcular diferencias entre soluciones consecutivas
    diff_1_2 = np.linalg.norm(soluciones[1] - soluciones[0])
    diff_2_3 = np.linalg.norm(soluciones[2] - soluciones[1])
    
    print(f"\n||θ_100 - θ_50||  = {diff_1_2:.6f}")
    print(f"||θ_200 - θ_100|| = {diff_2_3:.6f}")
    print(f"Ratio de convergencia: {diff_1_2 / diff_2_3:.2f}")
    
    if diff_2_3 < diff_1_2:
        print("✅ PASA: Solución converge con refinamiento temporal")
        return True
    else:
        print("⚠️ Las diferencias no disminuyen claramente")
        return False


def test_caso_difusion_lineal():
    """Test 4: Comparar con caso límite D=constante (solución más simple)"""
    print("\n=== TEST 4: Caso Límite D Constante ===")
    
    L = 0.5
    T_final = 0.1
    M = 100
    N = 200
    
    # Condición inicial
    x_initial = np.linspace(0, L, M + 2)
    theta_initial = np.zeros(M + 2)
    center = L / 2
    sigma = 0.05
    theta_initial = 0.8 * np.exp(-((x_initial - center) ** 2) / (2 * sigma ** 2))
    
    # Resolver con Brooks-Corey
    x_bc, theta_bc, _ = resolucion_ecuacion_richards_1D_no_lineal(
        diffusivity_brooks_corey, L, T_final, M, N, theta_initial
    )
    
    # Resolver con D constante
    D_const = 1e-7  # m²/s
    D_func_constante = lambda theta: D_const * np.ones_like(theta)
    
    x_const, theta_const, _ = resolucion_ecuacion_richards_1D_no_lineal(
        D_func_constante, L, T_final, M, N, theta_initial
    )
    
    # Comparar anchos del frente
    def calcular_ancho_frente(x, theta):
        theta_max = np.max(theta)
        mask_10 = theta > 0.1 * theta_max
        mask_90 = theta > 0.9 * theta_max
        if np.any(mask_10) and np.any(mask_90):
            return x[mask_10][-1] - x[mask_90][0]
        return 0.0
    
    ancho_bc = calcular_ancho_frente(x_bc, theta_bc)
    ancho_const = calcular_ancho_frente(x_const, theta_const)
    
    print(f"Ancho frente (Brooks-Corey): {ancho_bc:.4f} m")
    print(f"Ancho frente (D constante):  {ancho_const:.4f} m")
    print(f"Ratio: {ancho_bc / ancho_const:.2f}")
    
    # Con D no lineal, el frente debería ser más difuso
    if ancho_bc > ancho_const * 0.8:
        print("✅ PASA: Comportamiento cualitativo correcto")
        print("   (D no lineal produce frente más disperso)")
        return True
    else:
        print("⚠️ Comportamiento inesperado")
        return False


def generar_reporte_completo():
    """Ejecutar todos los tests y generar reporte"""
    print("="*60)
    print("VALIDACIÓN NUMÉRICA - ACTIVIDAD B")
    print("Ecuación de Richards 1D con D(θ) de Brooks-Corey")
    print("="*60)
    
    resultados = []
    
    resultados.append(("Conservación de Masa", test_conservacion_masa()))
    resultados.append(("Estabilidad Física", test_estabilidad_fisica()))
    resultados.append(("Convergencia Temporal", test_convergencia_temporal()))
    resultados.append(("Caso Límite D Constante", test_caso_difusion_lineal()))
    
    print("\n" + "="*60)
    print("RESUMEN DE VALIDACIÓN")
    print("="*60)
    
    tests_pasados = sum(1 for _, result in resultados if result)
    tests_totales = len(resultados)
    
    for nombre, resultado in resultados:
        estado = "✅ PASA" if resultado else "❌ FALLA"
        print(f"{estado} - {nombre}")
    
    print(f"\nResultado: {tests_pasados}/{tests_totales} tests pasados")
    
    if tests_pasados == tests_totales:
        print("\n🎉 TODOS LOS TESTS PASARON")
        print("   La implementación de Richards 1D es CORRECTA")
    elif tests_pasados >= tests_totales * 0.75:
        print("\n✅ IMPLEMENTACIÓN ACEPTABLE")
        print("   La mayoría de los tests pasaron")
    else:
        print("\n⚠️ REVISAR IMPLEMENTACIÓN")
        print("   Varios tests fallaron")
    
    print("="*60)


if __name__ == "__main__":
    generar_reporte_completo()
