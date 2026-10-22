#Codigo para el método de Newton-Raphson para encontrar raíces de polinomios de grado 6.
#laboratorio VI.1 de Calculo de Raices de Polinomios.
#Integrantes: Francisco Vecchione - Walter Ariel Paz Estrada
#URL de GitHub: https://github.com/Vecchio51/Lab.-VI.1-C-lculo-de-ra-ces-derivadas-sistemas-de-ecuaciones-lineales-e-integrales-con-GUI

import tkinter as tk
from tkinter import messagebox
import math

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

