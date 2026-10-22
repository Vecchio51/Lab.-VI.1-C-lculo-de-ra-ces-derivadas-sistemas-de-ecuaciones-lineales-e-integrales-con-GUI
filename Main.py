import tkinter as tk
from tkinter import messagebox
import math

###Evalúa el polinomio dado un x y una lista de coeficientes [a0, a1, ..., a6].

def evaluar_polinomio(x, coeficientes):
    resultado = 0
    for i, coef in enumerate(coeficientes):
        resultado += coef * (x ** i)
    return resultado

###Calcula la derivada analítica para el método de Newton-Raphson.

def derivada_analitica(x, coeficientes):
    resultado = 0
    for i in range(1, len(coeficientes)): # a0 desaparece al derivar
        resultado += coeficientes[i] * i * (x ** (i - 1))
    return resultado

###Calcula la tolerancia es a partir de las cifras significativas.

def calcular_es(cifras_sig):
    return 0.5 * (10 ** (2 - cifras_sig))

