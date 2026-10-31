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
        
        # --- Frame 1: Coeficientes del polinomio ---
        self.frame_coef = tk.LabelFrame(root, text="Coeficientes del polinomio", padx=5, pady=5)
        self.frame_coef.pack(fill="x", padx=10, pady=5)
        
        self.coef_entries = []
        valores_defecto_pol = [0, 1, 0, -8, 0, 82, -25] 
        etiquetas_coef = ["a6", "a5", "a4", "a3", "a2", "a1", "a0"]
        
        for i in range(7):
            lbl = tk.Label(self.frame_coef, text=etiquetas_coef[i])
            lbl.grid(row=1, column=i)
            ent = tk.Entry(self.frame_coef, width=5)
            ent.insert(0, str(valores_defecto_pol[i]))
            ent.grid(row=0, column=i, padx=2)
            self.coef_entries.append(ent)

        # --- Frame 2: Cálculo de una raíz ---
        self.frame_raiz = tk.LabelFrame(root, text="Cálculo de una raíz", padx=5, pady=5)
        self.frame_raiz.pack(fill="x", padx=10, pady=5)
        
        self.metodo_raiz_var = tk.StringVar(value="Bisección")
        tk.Radiobutton(self.frame_raiz, text="Bisección", variable=self.metodo_raiz_var, value="Bisección").grid(row=0, column=0)
        tk.Radiobutton(self.frame_raiz, text="Newton-Raphson", variable=self.metodo_raiz_var, value="Newton").grid(row=0, column=1)
        tk.Radiobutton(self.frame_raiz, text="Secante", variable=self.metodo_raiz_var, value="Secante").grid(row=0, column=2)
        
        tk.Label(self.frame_raiz, text="xl / x0 / x-1:").grid(row=1, column=0, sticky="e")
        self.ent_xi = tk.Entry(self.frame_raiz, width=8)
        self.ent_xi.insert(0, "0")
        self.ent_xi.grid(row=1, column=1, sticky="w")
        
        tk.Label(self.frame_raiz, text="xu / x0:").grid(row=1, column=2, sticky="e")
        self.ent_xu = tk.Entry(self.frame_raiz, width=8)
        self.ent_xu.insert(0, "1")
        self.ent_xu.grid(row=1, column=3, sticky="w")
        
        tk.Label(self.frame_raiz, text="Tolerancia es (%):").grid(row=2, column=0, sticky="e")
        self.ent_tol = tk.Entry(self.frame_raiz, width=8)
        self.ent_tol.insert(0, "0.5")
        self.ent_tol.grid(row=2, column=1, sticky="w")
        
        tk.Label(self.frame_raiz, text="Máx iter:").grid(row=2, column=2, sticky="e")
        self.ent_iter = tk.Entry(self.frame_raiz, width=8)
        self.ent_iter.insert(0, "100")
        self.ent_iter.grid(row=2, column=3, sticky="w")
        
        self.btn_calc_raiz = tk.Button(self.frame_raiz, text="Calcular", command=self.calcular_raiz)
        self.btn_calc_raiz.grid(row=1, column=4, rowspan=2, padx=10)
        
        self.lbl_res_raiz = tk.Label(self.frame_raiz, text="Raíz:\nIteraciones:\nError:\nf(raíz):", justify="left")
        self.lbl_res_raiz.grid(row=0, column=5, rowspan=3, sticky="w")

        # --- Frame 3: Derivada ---
        self.frame_deriv = tk.LabelFrame(root, text="Cálculo de la derivada primera (diferencia centrada)", padx=5, pady=5)
        self.frame_deriv.pack(fill="x", padx=10, pady=5)
        
        tk.Label(self.frame_deriv, text="x0:").pack(side="left")
        self.ent_deriv_x0 = tk.Entry(self.frame_deriv, width=5)
        self.ent_deriv_x0.insert(0, "1")
        self.ent_deriv_x0.pack(side="left", padx=5)
        
        tk.Label(self.frame_deriv, text="dx (h):").pack(side="left")
        self.ent_deriv_h = tk.Entry(self.frame_deriv, width=5)
        self.ent_deriv_h.insert(0, "0.1")
        self.ent_deriv_h.pack(side="left", padx=5)
        
        self.btn_calc_deriv = tk.Button(self.frame_deriv, text="Calcular", command=self.calcular_derivada)
        self.btn_calc_deriv.pack(side="left", padx=10)
        
        self.lbl_res_deriv = tk.Label(self.frame_deriv, text="Derivada f'(x0): ")
        self.lbl_res_deriv.pack(side="left")

        # --- Frame 4: Sistema de Ecuaciones ---
        self.frame_sis = tk.LabelFrame(root, text="Resolución de un sistema de 4 ecuaciones", padx=5, pady=5)
        self.frame_sis.pack(fill="x", padx=10, pady=5)
        self.matriz_entries = []
        self.b_entries = []
        
        valores_A = [[2,2,1,-1], [2,-2,6,-4], [2,2,1,1], [1,5,-4,-2]]
        valores_b = [5, 0, 13, -9]
        
        for i in range(4):
            fila_entries = []
            for j in range(4):
                ent = tk.Entry(self.frame_sis, width=5)
                ent.insert(0, str(valores_A[i][j]))
                ent.grid(row=i+1, column=j)
                fila_entries.append(ent)
            self.matriz_entries.append(fila_entries)
            
            ent_b = tk.Entry(self.frame_sis, width=5)
            ent_b.insert(0, str(valores_b[i]))
            ent_b.grid(row=i+1, column=4, padx=10)
            self.b_entries.append(ent_b)
            
        self.pivoteo_var = tk.BooleanVar(value=False)
        tk.Radiobutton(self.frame_sis, text="Sin pivoteo", variable=self.pivoteo_var, value=False).grid(row=1, column=5)
        tk.Radiobutton(self.frame_sis, text="Con pivoteo parcial", variable=self.pivoteo_var, value=True).grid(row=2, column=5)
        
        self.btn_calc_sis = tk.Button(self.frame_sis, text="Calcular", command=self.calcular_sistema)
        self.btn_calc_sis.grid(row=3, column=5)
        
        self.lbl_res_sis = tk.Label(self.frame_sis, text="x1=\nx2=\nx3=\nx4=", justify="left")
        self.lbl_res_sis.grid(row=1, column=6, rowspan=4, padx=10)

        # --- Frame 5: Integral ---
        self.frame_int = tk.LabelFrame(root, text="Cálculo de la integral definida (Simpson 1/3 compuesta)", padx=5, pady=5)
        self.frame_int.pack(fill="x", padx=10, pady=5)
        
        tk.Label(self.frame_int, text="Límite a:").pack(side="left")
        self.ent_int_a = tk.Entry(self.frame_int, width=5)
        self.ent_int_a.insert(0, "0")
        self.ent_int_a.pack(side="left", padx=5)
        
        tk.Label(self.frame_int, text="Límite b:").pack(side="left")
        self.ent_int_b = tk.Entry(self.frame_int, width=5)
        self.ent_int_b.insert(0, "2")
        self.ent_int_b.pack(side="left", padx=5)
        
        tk.Label(self.frame_int, text="Segmentos (n par):").pack(side="left")
        self.ent_int_n = tk.Entry(self.frame_int, width=5)
        self.ent_int_n.insert(0, "2")
        self.ent_int_n.pack(side="left", padx=5)
        
        self.btn_calc_int = tk.Button(self.frame_int, text="Calcular", command=self.calcular_integral)
        self.btn_calc_int.pack(side="left", padx=10)
        
        self.lbl_res_int = tk.Label(self.frame_int, text="Integral definida: ")
        self.lbl_res_int.pack(side="left")

    # --- Métodos de ejecución ---
    def obtener_coeficientes(self):
        return [float(e.get()) for e in self.coef_entries][::-1] # Invierte para que el índice coincida con el grado

    def calcular_raiz(self):
        try:
            coef = self.obtener_coeficientes()
            metodo = self.metodo_raiz_var.get()
            tol = float(self.ent_tol.get())
            max_iter = int(self.ent_iter.get())
            
            if metodo == "Bisección":
                xl = float(self.ent_xi.get())
                xu = float(self.ent_xu.get())
                xr, iteraciones, ea, fxr = biseccion(xl, xu, coef, tol, max_iter)
            elif metodo == "Newton":
                x0 = float(self.ent_xi.get())
                xr, iteraciones, ea, fxr = newton_raphson(x0, coef, tol, max_iter)
            elif metodo == "Secante":
                x_menos_1 = float(self.ent_xi.get())
                x0 = float(self.ent_xu.get())
                xr, iteraciones, ea, fxr = secante(x_menos_1, x0, coef, tol, max_iter)
                
            self.lbl_res_raiz.config(text=f"Raíz: {xr:.6f}\nIteraciones: {iteraciones}\nError: {ea:.4f}%\nf(raíz): {fxr:.6e}")
        except ValueError as e:
            messagebox.showwarning("Aviso", str(e))
        except Exception as e:
            messagebox.showerror("Error", f"Verifique los datos ingresados.\nDetalle: {e}")

    def calcular_derivada(self):
        try:
            coef = self.obtener_coeficientes()
            x0 = float(self.ent_deriv_x0.get())
            h = float(self.ent_deriv_h.get())
            res = derivada_centrada(x0, h, coef)
            self.lbl_res_deriv.config(text=f"Derivada f'(x0): {res:.6f}")
        except Exception as e:
            messagebox.showerror("Error", "Verifique los datos de la derivada.")

    def calcular_sistema(self):
        try:
            A = [[float(ent.get()) for ent in fila] for fila in self.matriz_entries]
            b = [float(ent.get()) for ent in self.b_entries]
            pivoteo = self.pivoteo_var.get()
            
            x = gauss(A, b, pivoteo)
            texto_res = "\n".join([f"x{i+1} = {val:.4f}" for i, val in enumerate(x)])
            self.lbl_res_sis.config(text=texto_res)
        except ValueError as e:
            messagebox.showwarning("Aviso de Sistema", str(e))
        except Exception as e:
            messagebox.showerror("Error", "Verifique los coeficientes del sistema.")

    def calcular_integral(self):
        try:
            coef = self.obtener_coeficientes()
            a = float(self.ent_int_a.get())
            b = float(self.ent_int_b.get())
            n = int(self.ent_int_n.get())
            
            res = simpson_13_compuesta(a, b, n, coef)
            self.lbl_res_int.config(text=f"Integral definida: {res:.6f}")
        except ValueError as e:
            messagebox.showwarning("Aviso", str(e))
        except Exception as e:
            messagebox.showerror("Error", "Verifique los datos de la integral.")

if __name__ == "__main__":
    root = tk.Tk()
    app = LaboratorioNumericoApp(root)
    root.mainloop()