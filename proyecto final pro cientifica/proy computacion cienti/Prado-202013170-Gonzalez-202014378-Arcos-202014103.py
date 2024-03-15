##

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg)
from scipy.integrate import odeint
from tkinter import filedialog
matplotlib.use("TkAgg")

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from PIL import ImageTk,Image
import struct as st

fventana= "#fff6ff"
# define el tamaño de la interfaz y el color de fondo
window=tk.Tk()
window.geometry('900x600')      # Tamaño de la ventana
window.config(bg=fventana)
window.title('Modelo matematico para proyecciones sobre la tuberculosis')

def ventana_cerrar():
    box = tk.messagebox.askokcancel("Cerrar Aplicación", "está apunto de salir de la aplicación, presione aceptar para continuar", icon="warning") #pregunta al presionar
    if box == True:
        window.destroy() #se cierra la ventana si da aceptar
    else:
        tk.messagebox.showinfo('Retornando a la aplicación', 'Gracias por seguir con nosotros') #retorna a la interfaz
# boton de cerrar
Cerrar= tk.Button(window, text='X', background='#ED3D1A', width=5, height=1,
                   command=ventana_cerrar,relief="ridge").place(x=1, y=1)

# Frames para tiempo, parametros, metodos y botones
frame_importar = tk.Frame(master=window) #crea el frame en nuestra ventana
frame_importar.place(x=100, y=10) #define la ubicacion del frame dentro de la interfaz
frame_importar.config(bg=fventana, width=600, height=100, relief="flat", bd=8)  #define otras caracteristicas del frame

frame_metodos = tk.Frame(master=window)
frame_metodos.place(x=100,y=430)
frame_metodos.config(bg=fventana, width=300, height=250, relief="flat", bd=8)

parametrosc = "#f2dee3"
frame_parametros = tk.Frame(master=window)
frame_parametros.place(x=600, y=10)
frame_parametros.config(bg=parametrosc, width=300, height=320, relief=tk.GROOVE, bd=8)

frame_circulos = tk.Frame(master=window)
frame_circulos.place(x=70,y=370)
frame_circulos.config(bg=fventana, width=400, height=80, relief="flat", bd=8)

frame_tiempo = tk.Frame(master=window)
frame_tiempo.place(x=605, y=400)
frame_tiempo.config(bg=fventana, width=300, height=100, relief="flat", bd=5)

frame_graficas = tk.Frame(master=window)
frame_graficas.place(x=70,y=100)
frame_graficas.config(bg="white", width=500, height=270, relief=tk.GROOVE, bd=8)

# Titulos para cada frame

titulo_metodos = tk.Label(master=frame_metodos,bg=fventana,
                   font=('Arial', 12,'bold italic'), text=f"Métodos de solución",height=-8).grid(row=0, column=1,pady=5)

titulo_parametros = tk.Label(master=frame_parametros, bg=parametrosc,
                   font=('Arial', 12,'bold italic'), text=f"Parámetros",height=2).grid(row=1,pady=0)

titulo_tiempo = tk.Label(master=frame_tiempo, bg=fventana, font=('Arial', 12,'bold italic'), text=f"Tiempo de simulación",height=2).grid(row=0, padx=25)


# Importar y exportar
letra = "white"

def browseFiles():
    filename = filedialog.askopenfilename(initialdir="/",
                                          title="Select a File",
                                          filetypes=(("Text files",
                                                      "*.txt*"),
                                                     ("all files",
                                                      "*.bin*")))

    # Change label contents
    exp=open('TiempoDeSimulación.bin','rb')
    leer=exp.read()
    desempacar=st.unpack('d'*int(len(leer)/8),leer)
    exp.close()
    T=np.array(desempacar)

    # Change label contents
    exp_y = open(filename, 'rb')
    l = exp_y.read()
    desempacar_y = st.unpack('d' * int(len(l)/8), l)
    S1EulerFor = desempacar_y
    exp_y.close()
    fig = plt.Figure(figsize=(4.8, 2.5), dpi=100)
    grafica = fig.add_subplot(111)

    grafica.plot(T,S1EulerFor,'r')
    Plot = FigureCanvasTkAgg(fig, master=frame_graficas)

    Plot.get_tk_widget().place(x=0, y=0)
    Plot.draw()


# botones para exportar e importar con tk button
colorei ="#F0C138"

importar= tk.Button(master=frame_importar, text='Importar', background=colorei
                    , foreground=letra, width=15, height=2,
                    command = browseFiles,font=('Arial', 12, 'bold italic'),relief="ridge").place(x=250, y=20)



# entry para los valores de cda parametro
fondo = '#a167c9'

# 1
A_label = tk.Label(master=frame_parametros, height=1, width=7, text="ʌ", foreground=letra,
                   background=fondo).grid(row=2,pady=3) #crea el label para nombre del parametro
alphau = tk.DoubleVar()  #define el tipo de parametro que entrara por el entry
A = tk.Entry(master=frame_parametros, textvariable=alphau, width=15).grid(row=2, column=1,pady=4, padx=20)#crea el entry

# 2
B_label = tk.Label(master=frame_parametros, height=1, width=7, text="β", background=fondo, foreground=letra).grid(row=3)
betau = tk.DoubleVar()
Beta = tk.Entry(master=frame_parametros, textvariable=betau, width =15).grid(row=3, column=1 ,pady=4)

# 3
delta_label =tk.Label(master=frame_parametros, height=1, width=7, text="δ", background=fondo,
                      foreground=letra).grid(row=4)
deltau = tk.DoubleVar()
Delta = tk.Entry(master=frame_parametros, textvariable=deltau, width = 15).grid(row=4, column=1,pady=4)

# 4
p_label = tk.Label(master=frame_parametros, height=1, width=7, text="p", background=fondo,
                   foreground=letra).grid(row=5)
pu = tk.DoubleVar()
P = tk.Entry(master=frame_parametros, textvariable=pu, width = 15).grid(row=5, column=1,pady=4)

# 5
u_label = tk.Label(master=frame_parametros, height=1, width=7, text="μ", background=fondo,
                   foreground=letra).grid(row=6)
uu = tk.DoubleVar()
U = tk.Entry(master=frame_parametros, textvariable=uu, width = 15).grid(row=6, column=1,pady=4)

# 6
k_label = tk.Label(master=frame_parametros, height=1, width=7, text="k", background=fondo, foreground=letra).grid(row=7)
kappau = tk.DoubleVar()
kappa= tk.Entry(master=frame_parametros, textvariable= kappau, width = 15).grid(row=7, column=1,pady=4)

# 7
r1_label = tk.Label(master=frame_parametros, height=1, width=7, text="r1", background=fondo,
                    foreground=letra).grid(row=8)
r1u = tk.DoubleVar()
R1 = tk.Entry(master=frame_parametros, textvariable=r1u, width = 15).grid(row=8, column=1,pady=4)

# 8
r2_label = tk.Label(master=frame_parametros, height=1, width=7, text="r2", background=fondo,
                    foreground=letra).grid(row=9)
r2u = tk.DoubleVar()
R2 = tk.Entry(master=frame_parametros, textvariable=r2u, width = 15).grid(row=9, column=1,pady=4)

# 9
phi_label = tk.Label(master=frame_parametros, height=1, width=7, text="φ", background=fondo,
                     foreground=letra).grid(row=10)
phiu = tk.DoubleVar()
Phi = tk.Entry(master=frame_parametros, textvariable=phiu, width = 15).grid(row=10, column=1,pady=4)


# 10
y_label = tk.Label(master=frame_parametros, height=1, width=7, text="γ", background=fondo,
                   foreground=letra).grid(row=11)
yu = tk.DoubleVar()
Y = tk.Entry(master=frame_parametros, textvariable=yu, width = 15).grid(row=11, column=1, pady=4)


# 11
d1_label = tk.Label(master=frame_parametros, height=1, width=7, text="d1", background=fondo,
                    foreground=letra).grid(row=12)
d1u = tk.DoubleVar()
D1= tk.Entry(master=frame_parametros, textvariable=d1u, width =15).grid(row=12, column=1,pady=4)

# 12
d2_label = tk.Label(master=frame_parametros, height=1, width=7, text="d2 ", background=fondo,
                    foreground=letra).grid(row=13)
d2u = tk.DoubleVar()
D2 = tk.Entry(master=frame_parametros, textvariable=d2u, width = 15).grid(row=13, column=1,pady=6)

# entry y frame para el tiempo
color_tiempo="#F0C138"
fondo_tiempo= fventana
tin = tk.DoubleVar()
tii = tk.Entry(master=frame_tiempo, textvariable=tin, width= 15,highlightbackground = color_tiempo,
              highlightthickness=2,highlightcolor= color_tiempo).grid(row=2)
li = tk.Label(master=frame_tiempo, width= 15,text="tiempo inicial", bg=fondo_tiempo, ).grid(row=1)
tfi = tk.DoubleVar()
tff = tk.Entry(master=frame_tiempo, textvariable=tfi, width = 15, highlightbackground = color_tiempo,
              highlightthickness=2,highlightcolor= color_tiempo).grid(row=4)
lfi = tk.Label(master=frame_tiempo, width= 15, bg=fondo_tiempo, text="tiempo final",).grid(row=3)
th = tk.DoubleVar()
hh = tk.Entry(master=frame_tiempo, textvariable=th, width= 15, highlightbackground = color_tiempo,
             highlightthickness=2, highlightcolor= color_tiempo).grid(row=6, pady=1)
lh = tk.Label(master=frame_tiempo, width= 15, text="años", bg=fondo_tiempo,).grid(row=5)



#valores ficticios por ahora--->QUITAR PARA LA ENTREGA FINAL
alphau.set(2)
betau.set(0.025)
deltau.set(1)
pu.set(0.3)
uu.set(0.0101)
kappau.set(0.005)
r1u.set(0)
r2u.set(0.8182)
phiu.set(0.02)
yu.set(0.01)
d1u.set(0.0227)
d2u.set(0.2)
tin.set(0.0)
tfi.set(20)
th.set(0.5)

alpha = alphau.get()
beta = betau.get()
delta = deltau.get()
u = uu.get()
p = pu.get()
r2 = r2u.get()
r1 = r1u.get()
y = yu.get()
d1 = d1u.get()
d2 = d2u.get()
phi = phiu.get()
kappa = kappau.get()
h = th.get()
ti = tin.get()
tf = tfi.get()
T = np.arange(ti, tf + h, h)
# Sistema de ecuaciones


S10 = 198 # se define el valor inicial de los suceptibles
E20 = 1  # se define el valor inicial de los expuestos
I30 = 0  # se define el valor inicial de los infectados
L40 = 0  # se define el valor inicial de los recuperados

def ecuacion_dS(S, I, L, alpha, beta, delta, u):
    return alpha - (beta * S * (I + delta * L)) - u * S


def ecuacion_dE(S, E, I, L, beta, delta, p, u, kappa, r1, r2):
    return (beta * (1 - p) * S * (I + delta * L)) + (r2 * I) - (u + kappa * (1 - r1)) * E


def ecuacion_dI(S, E, I, L, beta, delta, p, u, kappa, r1, r2, phi, y, d1):
    return beta * p * S * (I + delta * L) + (kappa * (1 - r1) * E) + (y * L) - (u + d1 + phi * (1 - r2) + r2) * I


def ecuacion_dL(I, L, u, r2, phi, y, d2):
    return (phi * (1 - r2) * I) - (u + d2 + y) * L



#funcion para solve ivp
def FSystem(Y,t):
    S = Y[0]
    E = Y[1]
    I = Y[2]
    L = Y[3]
    dsdt = alpha-beta*S*(I+delta*L)-u*S
    dedt = beta*(1-p)*S*(I+delta*L)+r2* I-(u+kappa*(1-r1))*E
    didt = beta*p*S*(I+delta*L)+kappa*(1-r1)*E+y*L-(u+d1+phi*(1-r2)+r2)*I
    dldt = phi*(1-r2)*I-(u+d2+y)*L
    z = [dsdt, dedt,didt,dldt]
    return z

#funcion para eular backward
def FEulerBackRoot(yt2, St1, Et1, It1, Lt1, alpha, beta, delta, p, u,kappa, r1, r2, phi, y, d1, d2):
    return [St1 + h * ecuacion_dS(yt2[0], yt2[2], yt2[3], alpha, beta, delta, u) - yt2[0],
            Et1 + h * ecuacion_dE(yt2[0], yt2[1], yt2[2], yt2[3], beta, delta, p, u, kappa, r1, r2) - yt2[1],
            It1 + h * ecuacion_dI(yt2[0], yt2[1], yt2[2], yt2[3], beta, delta, p, u, kappa, r1, r2, phi, y,d1) - yt2[2],
            Lt1 + h * ecuacion_dL(yt2[2], yt2[3], u, r2, phi, y, d2) - yt2[3]]

#funcion para euler modificado
def FEulerModRoot(yt2, St1, Et1, It1, Lt1, alpha, beta, delta, p, u, kappa, r1, r2, phi, y, d1, d2):
    return [St1 + (h / 2.0) * (ecuacion_dS(St1, It1, Lt1, alpha, beta, delta, u) +
                               ecuacion_dS(yt2[0], yt2[2], yt2[3], alpha, beta, delta, u)) - yt2[0],
            Et1 + (h / 2.0) * (ecuacion_dE(St1, Et1, It1, Lt1, beta, delta, p, u, kappa, r1, r2) +
                               ecuacion_dE(yt2[0], yt2[1], yt2[2], yt2[3], beta, delta, p, u, kappa, r1, r2)) - yt2[1],
            It1 + (h / 2.0) * (ecuacion_dI(St1, Et1, It1, Lt1, beta, delta, p, u, kappa, r1, r2, phi, y, d1) +
                               ecuacion_dI(yt2[0], yt2[1], yt2[2], yt2[3], beta, delta, p, u, kappa, r1, r2, phi, y,
                                           d1)) - yt2[2],
            Lt1 + (h / 2.0) * (ecuacion_dL(It1, Lt1, u, r2, phi, y, d2) +
                               ecuacion_dL(yt2[2], yt2[3], u, r2, phi, y, d2)) - yt2[3]]




# Crear vectores para guardar valores y guardar valores iniciales
# se crea el arreglo para guardar los valores calculados en cada iteracion de cada metodo para la variable de suceptibles
S1EulerFor = np.zeros(len(T))
S1EulerBack = np.zeros(len(T))
S1EulerMod = np.zeros(len(T))
S1RK2 = np.zeros(len(T))
S1RK4 = np.zeros(len(T))
# se crea el arreglo para guardar los valores calculados en cada iteracion de cada metodo para la variable de expuestos
E2EulerFor = np.zeros(len(T))
E2EulerBack = np.zeros(len(T))
E2EulerMod = np.zeros(len(T))
E2RK2 = np.zeros(len(T))
E2RK4 = np.zeros(len(T))
# se crea el arreglo para guardar los valores calculados en cada iteracion de cada metodo para la variable de infectados
I3EulerFor = np.zeros(len(T))
I3EulerBack = np.zeros(len(T))
I3EulerMod = np.zeros(len(T))
I3RK2 = np.zeros(len(T))
I3RK4 = np.zeros(len(T))
# se crea el arreglo para guardar los valores calculados en cada iteracion de cada metodo para la variable de recuperados/muertos
L4EulerFor = np.zeros(len(T))
L4EulerBack = np.zeros(len(T))
L4EulerMod = np.zeros(len(T))
L4RK2 = np.zeros(len(T))
L4RK4 = np.zeros(len(T))

# La primera posición del arreglo corresponde a los valores iniciales de los suceptibles
S1EulerFor[0] = S10
S1EulerBack[0] = S10
S1EulerMod[0] = S10
S1RK2[0] = S10
S1RK4[0] = S10
# La primera posición del arreglo corresponde a los valores iniciales de los expuestos
E2EulerFor[0] = E20
E2EulerBack[0] = E20
E2EulerMod[0] = E20
E2RK2[0] = E20
E2RK4[0] = E20
# La primera posición del arreglo corresponde a los valores iniciales de los infectados
I3EulerFor[0] = I30
I3EulerBack[0] = I30
I3EulerMod[0] = I30
I3RK2[0] = I30
I3RK4[0] = I30
# La primera posición del arreglo corresponde a los valores iniciales de los recuperados/muertos
L4EulerFor[0] = L40
L4EulerBack[0] = L40
L4EulerMod[0] = L40
L4RK2[0] = L40
L4RK4[0] = L40

#comienzo de la solucion

# solucion solve ivp
ICS = [S10, E20, I30, L40]
SolRK41 = odeint(FSystem, ICS, T)


for iter in range(1, len(T)):
    # Euler hacia adelante
    S1EulerFor[iter] = S1EulerFor[iter - 1] + h * ecuacion_dS(S1EulerFor[iter - 1], I3EulerFor[iter - 1],
                                                              L4EulerFor[iter - 1], alpha, beta, delta, u)
    E2EulerFor[iter] = E2EulerFor[iter - 1] + h * ecuacion_dE(S1EulerFor[iter - 1], E2EulerFor[iter - 1],
                                                              I3EulerFor[iter - 1], L4EulerFor[iter - 1],
                                                              beta, delta, p, u, kappa, r1, r2)
    I3EulerFor[iter] = I3EulerFor[iter - 1] + h * ecuacion_dI(S1EulerFor[iter - 1], E2EulerFor[iter - 1],
                                                              I3EulerFor[iter - 1], L4EulerFor[iter - 1],beta, delta,
                                                              p, u, kappa, r1, r2, phi, y, d1)
    L4EulerFor[iter] = L4EulerFor[iter - 1] + h * ecuacion_dL(I3EulerFor[iter - 1], L4EulerFor[iter - 1],
                                             u, r2, phi, y, d2)
    # Euler hacia atras
    # Se obtuvieron los coeficientes del euler
    SolBack = opt.fsolve(FEulerBackRoot,  # Función
                         np.array([S1EulerBack[iter - 1],  # x0
                                   E2EulerBack[iter - 1],
                                   I3EulerBack[iter - 1],
                                   L4EulerBack[iter - 1]]),
                         (S1EulerBack[iter - 1],  # Parámetros de función
                          E2EulerBack[iter - 1],
                          I3EulerBack[iter - 1],
                          L4EulerBack[iter - 1],
                          alpha,
                          beta,
                          delta,
                          p,
                          u,
                          kappa,
                          r1,
                          r2,
                          phi,
                          y,
                          d1,
                          d2),
                         xtol = 10 ** -15)  # Tolerancia

    S1EulerBack[iter] = SolBack[0] # Susceptibles
    E2EulerBack[iter] = SolBack[1] # Expuestos
    I3EulerBack[iter] = SolBack[2] # Infectados
    L4EulerBack[iter] = SolBack[3] # Recuperados/muertos
    # Euler modificado
    solMod = opt.fsolve(FEulerModRoot,
                        np.array([S1EulerMod[iter - 1],
                                  E2EulerMod[iter - 1],
                                  I3EulerMod[iter - 1],
                                  L4EulerMod[iter - 1]]),
                        (S1EulerMod[iter - 1],
                         E2EulerMod[iter - 1],
                         I3EulerMod[iter - 1],
                         L4EulerMod[iter - 1],
                         alpha,
                         beta,
                         delta,
                         p,
                         u, kappa, r1, r2, phi, y, d1, d2),
                        xtol=10 ** -5)
    S1EulerMod[iter] = solMod[0]
    E2EulerMod[iter] = solMod[1]
    I3EulerMod[iter] = solMod[2]
    L4EulerMod[iter] = solMod[3]
    # rk2
    k11 = ecuacion_dS(S1RK2[iter - 1], I3RK2[iter - 1], L4RK2[iter - 1], alpha, beta, delta, u)
    k21 = ecuacion_dE(S1RK2[iter - 1], E2RK2[iter - 1], I3RK2[iter - 1], L4RK2[iter - 1], beta, delta, p, u, kappa,
                      r1, r2)
    k31 = ecuacion_dI(S1RK2[iter - 1], E2RK2[iter - 1], I3RK2[iter - 1], L4RK2[iter - 1], beta, delta, p, u, kappa, r1,
                      r2, phi, y, d1)
    k41 = ecuacion_dL(I3RK2[iter - 1], L4RK2[iter - 1], u, r2, phi, y, d2)

    k12 = ecuacion_dS(S1RK2[iter - 1] + k11 * h, I3RK2[iter - 1] + k31 * h, L4RK2[iter - 1] + k41 * h, alpha, beta,
                      delta, u)
    k22 = ecuacion_dE(S1RK2[iter - 1] + k11 * h, E2RK2[iter - 1] + k21 * h, I3RK2[iter - 1] + k31 * h,
                      L4RK2[iter - 1] + k41 * h, beta, delta, p, u, kappa, r1, r2)
    k32 = ecuacion_dI(S1RK2[iter - 1] + k11 * h, E2RK2[iter - 1] + k21 * h, I3RK2[iter - 1] + k31 * h,
                      L4RK2[iter - 1] + k41 * h, beta, delta, p, u, kappa, r1, r2, phi, y, d1)
    k42 = ecuacion_dL(I3RK2[iter - 1] + k31 * h, L4RK2[iter - 1] + k41 * h, u, r2, phi, y, d2)

    S1RK2[iter] = S1RK2[iter - 1] + ((h / 2.0) * (k11 + k12))
    E2RK2[iter] = E2RK2[iter - 1] + ((h / 2.0) * (k21 + k22))
    I3RK2[iter] = I3RK2[iter - 1] + ((h / 2.0) * (k31 + k32))
    L4RK2[iter] = L4RK2[iter - 1] + ((h / 2.0) * (k41 + k42))

    # RK4

    k11 = ecuacion_dS(S1RK4[iter - 1], I3RK4[iter - 1], L4RK4[iter - 1], alpha, beta, delta, u)
    k21 = ecuacion_dE(S1RK4[iter - 1], E2RK4[iter - 1], I3RK4[iter - 1], L4RK4[iter - 1], beta, delta, p, u, kappa, r1,
                      r2)
    k31 = ecuacion_dI(S1RK4[iter - 1], E2RK4[iter - 1], I3RK4[iter - 1], L4RK4[iter - 1], beta, delta, p, u, kappa, r1,
                      r2, phi, y, d1)
    k41 = ecuacion_dL(I3RK4[iter - 1], L4RK4[iter - 1], u, r2, phi, y, d2)

    k12 = ecuacion_dS(S1RK4[iter - 1] + (0.5 * k11 * h), I3RK4[iter - 1] + (0.5 * k31 * h),
                      L4RK4[iter - 1] + (0.5 * k41 * h), alpha, beta, delta, u)
    k22 = ecuacion_dE(S1RK4[iter - 1] + (0.5 * k11 * h), E2RK4[iter - 1] + (0.5 * k21 * h),
                      I3RK4[iter - 1] + (0.5 * k31 * h), L4RK4[iter - 1] + (0.5 * k41 * h), beta, delta, p, u, kappa,
                      r1, r2)
    k32 = ecuacion_dI(S1RK4[iter - 1] + (0.5 * k11 * h), E2RK4[iter - 1] + (0.5 * k21 * h),
                      I3RK4[iter - 1] + (0.5 * k31 * h), L4RK4[iter - 1] + (0.5 * k41 * h), beta, delta, p, u, kappa,
                      r1, r2, phi, y, d1)
    k42 = ecuacion_dL(I3RK4[iter - 1] + (0.5 * k31 * h), L4RK4[iter - 1] + (0.5 * k41 * h), u, r2, phi, y, d2)

    k13 = ecuacion_dS(S1RK4[iter - 1] + (0.5 * k12 * h), I3RK4[iter - 1] + (0.5 * k32 * h),
                      L4RK4[iter - 1] + (0.5 * k42 * h), alpha, beta, delta, u)
    k23 = ecuacion_dE(S1RK4[iter - 1] + (0.5 * k12 * h), E2RK4[iter - 1] + (0.5 * k22 * h),
                      I3RK4[iter - 1] + (0.5 * k32 * h), L4RK4[iter - 1] + (0.5 * k42 * h), beta, delta, p, u, kappa,
                      r1, r2)
    k33 = ecuacion_dI(S1RK4[iter - 1] + (0.5 * k12 * h), E2RK4[iter - 1] + (0.5 * k22 * h),
                      I3RK4[iter - 1] + (0.5 * k32 * h), L4RK4[iter - 1] + (0.5 * k42 * h), beta, delta, p, u, kappa,
                      r1, r2, phi, y, d1)
    k43 = ecuacion_dL(I3RK4[iter - 1] + (0.5 * k32 * h), L4RK4[iter - 1] + (0.5 * k42 * h), u, r2, phi, y, d2)

    k14 = ecuacion_dS(S1RK4[iter - 1] + (0.5 * k13 * h), I3RK4[iter - 1] + (0.5 * k33 * h),
                      L4RK4[iter - 1] + (0.5 * k43 * h), alpha, beta, delta, u)
    k24 = ecuacion_dE(S1RK4[iter - 1] + (0.5 * k13 * h), E2RK4[iter - 1] + (0.5 * k23 * h),
                      I3RK4[iter - 1] + (0.5 * k33 * h), L4RK4[iter - 1] + (0.5 * k43 * h), beta, delta, p, u, kappa,
                      r1, r2)
    k34 = ecuacion_dI(S1RK4[iter - 1] + (0.5 * k13 * h), E2RK4[iter - 1] + (0.5 * k23 * h),
                      I3RK4[iter - 1] + (0.5 * k33 * h), L4RK4[iter - 1] + (0.5 * k43 * h), beta, delta, p, u, kappa,
                      r1, r2, phi, y, d1)
    k44 = ecuacion_dL(I3RK4[iter - 1] + (0.5 * k33 * h), L4RK4[iter - 1] + (0.5 * k43 * h), u, r2, phi, y, d2)

    S1RK4[iter] = S1RK4[iter - 1] + (h / 6) * (k11 + 2 * k12 + 2 * k13 + k14)
    E2RK4[iter] = E2RK4[iter - 1] + (h / 6) * (k21 + 2 * k22 + 2 * k23 + k24)
    I3RK4[iter] = I3RK4[iter - 1] + (h / 6) * (k31 + 2 * k32 + 2 * k33 + k34)
    L4RK4[iter] = L4RK4[iter - 1] + (h / 6) * (k41 + 2 * k42 + 2 * k43 + k44)

def resultados(boton):
    s = ""
    e = ""
    i = ""
    l = ""
    fig = plt.Figure(figsize=(4.8, 2.5), dpi=100)
    grafica = fig.add_subplot(111)
    if boton == 1:
        s = S1EulerFor
        e = E2EulerFor
        i = I3EulerFor
        l = L4EulerFor
    elif boton==2:
        s = S1EulerBack
        e = E2EulerBack
        i = I3EulerBack
        l = L4EulerBack
    elif boton==3:
        s = S1EulerMod
        e = E2EulerMod
        i = I3EulerMod
        l = L4EulerMod
    elif boton == 4:
        s = S1RK2
        e = E2RK2
        i = I3RK2
        l = L4RK2
    elif boton == 5:
        s = S1RK4
        e = E2RK4
        i = I3RK4
        l = L4RK4
    elif boton == 6:
        s = SolRK41[:,0]
        e = SolRK41[:,1]
        i = SolRK41[:,2]
        l = SolRK41[:,3]

    if CheckVarS.get() == 1:
        grafica.plot(T, s, 'r')
    if CheckVarE.get() == 1:
        grafica.plot(T, e, 'b')
    if CheckVarI.get() == 1:
        grafica.plot(T,i, 'm')
    if CheckVarL.get() == 1:
        grafica.plot(T, l, 'g')
    grafica.legend(['s(t)', 'e(t)', 'i(t)', 'l(t)'])

    plt.close()
    Plot = FigureCanvasTkAgg(fig, master=frame_graficas)
    Plot.draw()
    Plot.get_tk_widget().place(x=0, y=0)




#QUÉ MÉTODO USAR
colorf="#8148B6"
fuente= ('Arial', 8, 'bold italic')
Euler_adelante = tk.Button(master=frame_metodos, text='Euler adelante',command= lambda : [resultados(1),exportar(1)], background=colorf,
                          foreground=letra, width=15, height=2,relief="ridge", font=fuente).grid(row=2, column=0,pady=5)
Euler_atras = tk.Button(master=frame_metodos, text='Euler atrás',command=lambda : [resultados(2),exportar(2)], background=colorf,
                       foreground=letra, width=15, height=2,relief="ridge", font=fuente).grid(row=2, column=1,pady=5)
Euler_modificado = tk.Button(master=frame_metodos, text='Euler modificado',command=lambda : [resultados(3),exportar(3)], background=colorf,
                             foreground=letra, width=15, height=2,relief="ridge", font=fuente).grid(row=2, column=2,pady=5)
Runge_Kutta2 = tk.Button(master=frame_metodos, text='Runge-Kutta 2', command=lambda : [resultados(4),exportar(4)], background=colorf,
                        foreground=letra, width=15, height=2,relief="ridge", font=fuente).grid(row=3, column=0,pady=5)
Runge_Kutta4 = tk.Button(master=frame_metodos, text='Runge-Kutta 4',command=lambda : [resultados(5),exportar(5)], background=colorf,
                         foreground=letra, width=15, height=2,relief="ridge", font=fuente).grid(row=3, column=1,pady=5)
solve_ivp = tk.Button(master=frame_metodos, text='solve_ivp',command= lambda: [resultados(6),exportar(6)], background=colorf,
                     foreground=letra, width=15, height=2,relief="ridge", font=fuente).grid(row=3, column=2,pady=5)

#BOTÓN EXPORTAR
#colorei ="#F0C138"
def exportar(boton):
    datos=open('TiempoDeSimulación.bin','wb')
    empacar_t=st.pack('d'*int(len(T)),*T)
    datos.write(empacar_t)
    datos.close()
    datos_s= open('DatosEnS.bin', 'wb')
    datos_e = open('DatosEnE.bin', 'wb')
    datos_i = open('DatosEnI.bin', 'wb')
    datos_l = open('DatosEnL.bin', 'wb')
    s = ""
    e = ""
    i = ""
    l = ""
    if boton == 1:
        s = S1EulerFor
        e = E2EulerFor
        i = I3EulerFor
        l = L4EulerFor
    elif boton == 2:
        s = S1EulerBack
        e = E2EulerBack
        i = I3EulerBack
        l = L4EulerBack
    elif boton == 3:
        s = S1EulerMod
        e = E2EulerMod
        i = I3EulerMod
        l = L4EulerMod
    elif boton == 4:
        s = S1RK2
        e = E2RK2
        i = I3RK2
        l = L4RK2
    elif boton == 5:
        s = S1RK4
        e = E2RK4
        i = I3RK4
        l = L4RK4
    elif boton == 6:
        s = SolRK41[:, 0]
        e = SolRK41[:, 1]
        i = SolRK41[:, 2]
        l = SolRK41[:, 3]

    empacar_s = st.pack('d' * int(len(s)), *s)
    datos_s.write(empacar_s)
    datos_s.close()
    empacar_e = st.pack('d' * int(len(e)), *e)
    datos_e.write(empacar_e)
    datos_e.close()
    empacar_i = st.pack('d' * int(len(i)), *i)
    datos_i.write(empacar_i)
    datos_i.close()
    empacar_l = st.pack('d' * int(len(l)), *l)
    datos_l.write(empacar_l)
    datos_l.close()





Exportar = tk.Button(master=frame_importar, text='Exportar', background=colorei,foreground=letra , width=15, height=2,
                   font=('Arial', 12, 'bold italic'),relief="ridge",command=exportar).place(x=2, y=20)


# BOTONES REDONDOS

CheckVarS=tk.IntVar()
CheckVarE=tk.IntVar()
CheckVarI=tk.IntVar()
CheckVarL=tk.IntVar()

#crea el boton y lo ubica en el frame
botonS =tk.PhotoImage(file = 'BotonS.png')
botonl =tk.PhotoImage(file = 'BotonL.png')
botoni =tk.PhotoImage(file = 'Botoni.png')
botonE =tk.PhotoImage(file = 'BotonE.png')
S=tk.Checkbutton(master=frame_circulos,image=botonS, command=resultados,relief=tk.GROOVE,height=28, bd=4, onvalue=1,
                 offvalue=0,variable=CheckVarS, width=28, bg ="#E0CEFA" ).place(x=80,y=1)
E=tk.Checkbutton(master=frame_circulos,image=botonE,command=resultados,relief=tk.GROOVE,height=28,bd=4, onvalue=1,
                 offvalue=0,variable=CheckVarE, width=28,  bg ="#E0CEFA").place(x=160,y=1)
I=tk.Checkbutton(master=frame_circulos,image=botoni,command=resultados,relief=tk.GROOVE,height=28,bd=4, onvalue=1,
                 offvalue=0,variable=CheckVarI, width=28,  bg ="#E0CEFA").place(x=240,y=1)
L=tk.Checkbutton(master=frame_circulos,image=botonl,command=resultados,relief=tk.GROOVE,height=28, bd=4, onvalue=1,
                 offvalue=0,variable=CheckVarL, width=28,  bg ="#E0CEFA").place(x=320,y=1)


window.mainloop()


