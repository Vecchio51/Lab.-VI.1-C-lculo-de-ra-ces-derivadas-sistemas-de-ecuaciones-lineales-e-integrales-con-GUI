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
    
    xr_ant = 0
    for i in range(max_iter):
        xr = (xl + xu) / 2
        fxr = evaluar_polinomio(xr, coef) # Evalúa el polinomio en el punto medio
        
        ea = abs((xr - xr_ant) / xr) * 100 if xr != 0 else 0
        if i > 0 and ea < tol:
            return xr, i+1, ea, fxr
            
        test = evaluar_polinomio(xl, coef) * fxr # evalúa el producto de los valores de la función en los extremos del intervalo y el punto medio
        if test < 0:
            xu = xr
        elif test > 0: 
            xl = xr 
        else:
            ea = 0
            break
        xr_ant = xr
    raise ValueError("Se alcanzó el máximo de iteraciones sin converger.")

#-------------------------------------------------------------------------
#Funciones de derivacion y metodo de Simpson 1/3 compuesta para integrales:
#-------------------------------------------------------------------------

def derivada_centrada(x, h, coef):
    fx_mas_h = evaluar_polinomio(x + h, coef)
    fx_menos_h = evaluar_polinomio(x - h, coef)
    return (fx_mas_h - fx_menos_h) / (2 * h)

def simpson_13_compuesta(a, b, n, coef):
    #ACA SIGUE EL CODIGO PARA EL METODO DE SIMPSON 1/3 COMPUESTA, QUE CALCULA LA INTEGRAL DEFINIDA DE UN POLINOMIO


#-----------------------------------------------------------
#Funcion de Sistemas de Ecuaciones (Gauss) con y sin pivoteo:
#-----------------------------------------------------------

def gauss(A_input, b_input, pivoteo=False):
    # Copias locales para no mutar la entrada original
    A = [fila[:] for fila in A_input]
    b = b_input[:]
    n = len(b)
    #ACA SIGUE EL CODIGO PARA EL METODO DE GAUSS, QUE RESUELVE SISTEMAS DE ECUACIONES LINEALES, CON Y SIN PIVOTEO PARCIAL.
    #Si se requiere pivoteo parcial, se realiza antes de la eliminación
    #Hay que implementar la lógica de pivoteo parcial, eliminación hacia adelante y sustitución hacia atrás para obtener la solución del sistema de ecuaciones.

#--------------------------------------------------------------------
#Interfaz gráfica con Tkinter para el laboratorio de cálculo numérico:
#--------------------------------------------------------------------

class LaboratorioNumericoApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Laboratorio VI.1 - Cálculo Numérico")
#ACA  SE DEFINE EL CODIGO PARA ARMAR LA INTERFAZ GRAFICA, CON LOS CAMPOS DE ENTRADA PARA LOS COEFICIENTES DEL POLINOMIO, EL METODO A UTILIZAR, Y LOS BOTONES PARA EJECUTAR EL CALCULO DE RAICES.