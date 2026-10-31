#Codigo para encontrar raíces de polinomios de grado 6 con interfaz grafica.
#laboratorio VI.1 de Calculo de Raices de Polinomios.
#Integrantes: Francisco Vecchione - Walter Ariel Paz Estrada
#URL de GitHub: https://github.com/Vecchio51/Lab.-VI.1-C-lculo-de-ra-ces-derivadas-sistemas-de-ecuaciones-lineales-e-integrales-con-GUI

import tkinter as tk
from tkinter import messagebox
import math

#--------------------------------------------------------------------
#Funciones auxiliares para evaluar polinomios, derivadas y tolerancia:
#--------------------------------------------------------------------

###Evalúa el polinomio dado un x y una lista de coeficientes [a0, a1, ..., a6].

def evaluar_polinomio(x, coeficientes):
    resultado = 0          # Inicializa el resultado
    for i, coef in enumerate(coeficientes):
        resultado += coef * (x ** i) # Evalúa el término correspondiente
    return resultado

###Calcula la derivada analítica para el método de Newton-Raphson.

def derivada_analitica(x, coeficientes):
    resultado = 0          # Inicializa el resultado
    for i in range(1, len(coeficientes)): # a0 desaparece al derivar
        resultado += coeficientes[i] * i * (x ** (i - 1)) # Evalúa el término correspondiente
    return resultado

###Calcula la tolerancia es a partir de las cifras significativas.

def calcular_es(cifras_sig):
    return 0.5 * (10 ** (2 - cifras_sig)) # Evalúa la tolerancia según las cifras significativas

#----------------------------------------------
#metodos de bisección, Newton-Raphson y secante:
#----------------------------------------------

def biseccion(xl, xu, coef, tol, max_iter):
    if evaluar_polinomio(xl, coef) * evaluar_polinomio(xu, coef) > 0:
        raise ValueError("El intervalo no presenta cambio de signo.") # Verifica que el intervalo tenga un cambio de signo
    
    xr_ant = 0 # Inicializa el valor anterior de xr para calcular el error aproximado
    for i in range(max_iter):
        xr = (xl + xu) / 2 # Calcula el punto medio del intervalo
        fxr = evaluar_polinomio(xr, coef) # Evalúa el polinomio en el punto medio
        
        ea = abs((xr - xr_ant) / xr) * 100 if xr != 0 else 0 # Calcula el error aproximado
        if i > 0 and ea < tol: # Verifica si el error es menor que la tolerancia
            return xr, i+1, ea, fxr # Devuelve la raíz, el número de iteraciones, el error aproximado y el valor del polinomio en la raíz
            
        test = evaluar_polinomio(xl, coef) * fxr # evalúa el producto de los valores de la función en los extremos del intervalo y el punto medio
        if test < 0:
            xu = xr # Actualiza el límite superior del intervalo
        elif test > 0: 
            xl = xr # Actualiza el límite inferior del intervalo
        else:
            ea = 0 # Si el producto es cero, significa que se encontró una raíz exacta
            break
        xr_ant = xr # Actualiza el valor anterior de xr para la siguiente iteración
    raise ValueError("Se alcanzó el máximo de iteraciones sin converger.")
def newton_raphson(x0, coef, tol, max_iter):
    xr_ant = x0
    for i in range(max_iter):
        fx = evaluar_polinomio(xr_ant, coef)
        dfx = derivada_analitica(xr_ant, coef)
        if dfx == 0:
            raise ValueError("Derivada nula. No se puede continuar.")
        
        xr = xr_ant - (fx / dfx)
        ea = abs((xr - xr_ant) / xr) * 100 if xr != 0 else 0
        fxr = evaluar_polinomio(xr, coef)
        
        if ea < tol:
            return xr, i + 1, ea, fxr
        xr_ant = xr
        
    raise ValueError("Se alcanzó el máximo de iteraciones sin converger.")

def secante(x_menos_1, x0, coef, tol, max_iter):
    for i in range(max_iter):
        fx0 = evaluar_polinomio(x0, coef)
        fx_menos_1 = evaluar_polinomio(x_menos_1, coef)
        denominador = fx0 - fx_menos_1
        
        if denominador == 0:
            raise ValueError("División por cero. No se puede continuar.")
            
        xr = x0 - fx0 * (x0 - x_menos_1) / denominador
        ea = abs((xr - x0) / xr) * 100 if xr != 0 else 0
        fxr = evaluar_polinomio(xr, coef)
        
        if ea < tol:
            return xr, i + 1, ea, fxr
            
        x_menos_1 = x0
        x0 = xr
        
    raise ValueError("Se alcanzó el máximo de iteraciones sin converger.")

def simpson_13_compuesta(a, b, n, coef):
    if n % 2 != 0:
        raise ValueError("La cantidad de segmentos (n) debe ser par.")
        
    h = (b - a) / n
    suma = evaluar_polinomio(a, coef) + evaluar_polinomio(b, coef)
    
    for i in range(1, n):
        x_i = a + i * h
        if i % 2 == 0:
            suma += 2 * evaluar_polinomio(x_i, coef)
        else:
            suma += 4 * evaluar_polinomio(x_i, coef)
            
    return (h / 3) * suma
#-------------------------------------------------------------------------
#Funciones de derivacion y metodo de Simpson 1/3 compuesta para integrales:
#-------------------------------------------------------------------------

def derivada_centrada(x, h, coef):
    fx_mas_h = evaluar_polinomio(x + h, coef)
    fx_menos_h = evaluar_polinomio(x - h, coef)
    return (fx_mas_h - fx_menos_h) / (2 * h)

def simpson_13_compuesta(a, b, n, coef):
    if n % 2 != 0:
        raise ValueError("La cantidad de segmentos (n) debe ser par.")
    #ACA SIGUE EL CODIGO PARA EL METODO DE SIMPSON 1/3 COMPUESTA, QUE CALCULA LA INTEGRAL DEFINIDA DE UN POLINOMIO


#-----------------------------------------------------------
#Funcion de Sistemas de Ecuaciones (Gauss) con y sin pivoteo:
#-----------------------------------------------------------

def gauss(A_input, b_input, pivoteo=False):
    # Copias locales para no mutar la entrada original
    A = [fila[:] for fila in A_input]
    b = b_input[:]
    n = len(b)
    
    # Eliminación hacia adelante
    for i in range(n):
        if pivoteo:
            # Búsqueda del mayor elemento en la columna para pivotar
            max_row = i
            for k in range(i + 1, n):
                if abs(A[k][i]) > abs(A[max_row][i]):
                    max_row = k
            # Intercambio de filas si es necesario
            if max_row != i:
                A[i], A[max_row] = A[max_row], A[i]
                b[i], b[max_row] = b[max_row], b[i]

        # Detección de pivote nulo
        if A[i][i] == 0:
            raise ValueError(f"Pivote nulo detectado en la fila {i+1}. Seleccione 'Con pivoteo parcial' o modifique el sistema.")

        # Hacer ceros debajo del pivote
        for j in range(i + 1, n):
            factor = A[j][i] / A[i][i]
            for k in range(i, n):
                A[j][k] -= factor * A[i][k]
            b[j] -= factor * b[i]

    # Sustitución hacia atrás
    x = [0] * n
    for i in range(n - 1, -1, -1):
        suma = sum(A[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (b[i] - suma) / A[i][i]

    return x

#--------------------------------------------------------------------
#Interfaz gráfica con Tkinter para el laboratorio de cálculo numérico:
#--------------------------------------------------------------------

class LaboratorioNumericoApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Laboratorio VI.1 - Cálculo Numérico")
#ACA  SE DEFINE EL CODIGO PARA ARMAR LA INTERFAZ GRAFICA, CON LOS CAMPOS DE ENTRADA PARA LOS COEFICIENTES DEL POLINOMIO, EL METODO A UTILIZAR, Y LOS BOTONES PARA EJECUTAR EL CALCULO DE RAICES.