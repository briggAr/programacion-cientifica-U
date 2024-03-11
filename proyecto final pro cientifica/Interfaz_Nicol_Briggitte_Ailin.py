##
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg)

matplotlib.use("TkAgg")

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from PIL import ImageTk,Image

#define el tamaño de la interfaz y el color de fondo
window=tk.Tk()
window.geometry('1000x600')      # Tamaño de la ventana
window.config(bg='white')

def importar():
    return None
#boton de cerrar
Cerrar= tk.Button(window, text='X', background='#ED3D1A', width=5, height=1,
                   command=window.destroy,relief="ridge").place(x=1, y=1)

#Frames para tiempo, parametros, metodos y botones
frame_importar = tk.Frame(master=window) #crea el frame en nuestra ventana
frame_importar.place(x=80, y=10) #define la ubicacion del frame dentro de la interfaz
frame_importar.config(bg="white", width=600, height=100, relief="flat", bd=8)  #define otras caracteristicas del frame

frame_metodos = tk.Frame(master=window)
frame_metodos.place(x=620, y=400)
frame_metodos.config(bg="white", width=300, height=250, relief="flat", bd=8)


frame_parametros = tk.Frame(master=window)
frame_parametros.place(x=650, y=10)
frame_parametros.config(bg="#77cfe7", width=300, height=320, relief=tk.GROOVE, bd=8)

frame_circulos=tk.Frame(master=window)
frame_circulos.place(x=70,y=350)
frame_circulos.config(bg="white", width=400, height=80, relief="flat", bd=8)

frame_tiempo=tk.Frame(master=window)
frame_tiempo.place(x=80,y=480)
frame_tiempo.config(bg="white", width=300, height=100, relief="flat", bd=8)


#Importar y exportar
letra = "white"
#botones para exportar e importar con tk button
Exportar = tk.Button(master=frame_importar, text='Exportar', background="#2159AD",foreground=letra , width=15, height=2,
                   command=importar,font=('Arial', 12, 'bold italic'),relief="ridge").place(x=2, y=20)
Importar= tk.Button(master=frame_importar, text='Importar', background="#2159AD", foreground=letra, width=15, height=2,
                   command=importar,font=('Arial', 12, 'bold italic'),relief="ridge").place(x=250, y=20)


#funcion para grafica para cada metodo
def grafica():
    fig = plt.Figure(figsize=(4, 2), dpi=100)
    t = np.arange(0,10, 0.01)
    fig.add_subplot(111).plot(t, fun(t))     # subplot(filas, columnas, item)
    fig.suptitle(opcion.get())

    plt.close()
    plt.style.use('seaborn-darkgrid')
    Plot = FigureCanvasTkAgg(fig, master=window)
    Plot.draw()

#botones para metodos
colorf="#fad000"

Euler_adelante= tk.Button(master=frame_metodos, text='Euler adelante',command=grafica, background=colorf, foreground=letra, width=20, height=2,relief="ridge").grid(row=2, column=0,pady=5, padx=10)
Euler_atras= tk.Button(master=frame_metodos, text='Euler atrás',command=grafica, background=colorf, foreground=letra, width=20, height=2,relief="ridge").grid(row=3, column=0,pady=5)
Runge_Kutta2= tk.Button(master=frame_metodos, text='Runge-Kutta 2', command=grafica, background=colorf, foreground=letra, width=20, height=2,relief="ridge").grid(row=2, column=1,pady=5)
Runge_Kutta4= tk.Button(master=frame_metodos, text='Runge-Kutta 4',command=grafica, background=colorf, foreground=letra, width=20, height=2,relief="ridge").grid(row=3, column=1,pady=5)




#entry para los valores de cda parametro
fondo = '#009dc8'

# 1
A_label = tk.Label(master=frame_parametros, height=1, width=7, text="ʌ", foreground=letra, background=fondo).grid(row=2,pady=3, padx= 20) #crea el label para nombre del parametro
tasa_A = tk.StringVar()  #define el tipo de parametro que entrara por el entry
A = tk.Entry(master=frame_parametros, textvariable=tasa_A, width=20).grid(row=2, column=1,pady=4, padx=20) #crea el entry

#2
B_label = tk.Label(master=frame_parametros, height=1, width=7, text="β", background=fondo, foreground=letra).grid(row=3)
beta= tk.StringVar()
Beta = tk.Entry(master=frame_parametros, textvariable=beta, width =20).grid(row=3, column=1 ,pady=4)

#3
delta_label=tk.Label(master=frame_parametros, height=1, width=7, text="δ", background=fondo, foreground=letra).grid(row=4)
tasa_delta= tk.StringVar()
Delta= tk.Entry(master=frame_parametros, textvariable=tasa_delta, width = 20).grid(row=4, column=1,pady=4)

#4
p_label=tk.Label(master=frame_parametros, height=1, width=7, text="p", background=fondo, foreground=letra).grid(row=5)
proporcion= tk.StringVar()
P= tk.Entry(master=frame_parametros, textvariable=proporcion, width = 20).grid(row=5, column=1,pady=4)

#5
u_label=tk.Label(master=frame_parametros, height=1, width=7, text="μ", background=fondo, foreground=letra).grid(row=6)
tasa_u= tk.StringVar()
U= tk.Entry(master=frame_parametros, textvariable=tasa_u, width = 20).grid(row=6, column=1,pady=4)

#6
k_label=tk.Label(master=frame_parametros, height=1, width=7, text="k", background=fondo, foreground=letra).grid(row=7)
tasa_k= tk.StringVar()
K= tk.Entry(master=frame_parametros, textvariable=tasa_k, width = 20).grid(row=7, column=1,pady=4)

#7
r1_label=tk.Label(master=frame_parametros, height=1, width=7, text="r1", background=fondo, foreground=letra).grid(row=8)
tasa_r1= tk.StringVar()
r1= tk.Entry(master=frame_parametros, textvariable=tasa_r1, width = 20).grid(row=8, column=1,pady=4)

#8
r2_label=tk.Label(master=frame_parametros, height=1, width=7, text="r2", background=fondo, foreground=letra).grid(row=9)
tasa_r2= tk.StringVar()
r2= tk.Entry(master=frame_parametros, textvariable=tasa_r2, width = 20).grid(row=9, column=1,pady=4)

#9
phi_label=tk.Label(master=frame_parametros, height=1, width=7, text="φ", background=fondo, foreground=letra).grid(row=10)
tasa_phi= tk.StringVar()
Phi= tk.Entry(master=frame_parametros, textvariable=tasa_phi, width = 20).grid(row=10, column=1,pady=4)


#10
y_label=tk.Label(master=frame_parametros, height=1, width=7, text="γ", background=fondo, foreground=letra).grid(row=11)
tasa_y= tk.StringVar()
y= tk.Entry(master=frame_parametros, textvariable=tasa_y, width = 20).grid(row=11, column=1,pady=4)

#11
d1_label=tk.Label(master=frame_parametros, height=1, width=7, text="d1", background=fondo, foreground=letra).grid(row=12)
tasa_d1 = tk.StringVar()
D1= tk.Entry(master=frame_parametros, textvariable=tasa_d1, width = 20).grid(row=12, column=1,pady=4)

#12
d2_label=tk.Label(master=frame_parametros, height=1, width=7, text="d2 ", background=fondo, foreground=letra).grid(row=13)
tasa_d2= tk.StringVar()
D2= tk.Entry(master=frame_parametros, textvariable=tasa_d2, width = 20).grid(row=13, column=1,pady=6)

#entry y frame para el tiempo en dia, mes y años
tdia= tk.IntVar()
dia = tk.Entry(master=frame_tiempo, textvariable=tdia, width= 15).grid(row=0, column=0,padx=15)
tmes= tk.IntVar()
mes = tk.Entry(master=frame_tiempo, textvariable=tmes, width= 15).grid(row=0, column=1,padx=15)
tanio= tk.IntVar()
anio = tk.Entry(master=frame_tiempo, textvariable=tanio, width= 15).grid(row=0, column=2,padx=15)
lanio = tk.Label(master=frame_tiempo, width= 4, text="años", bg="white").grid(row=0, column=3)


# BOTONES REDONDOS

boton= tk.PhotoImage(file='Boton.PNG')
s_t=tk.Button(master=frame_circulos, image=boton, command=None,relief="flat",height=28, width=28).place(x=80,y=1) #crea el boton y lo ubica en el frame
E_t=tk.Button(master=frame_circulos, image=boton, command=None,relief="flat",height=28, width=28).place(x=160,y=1)
I_t=tk.Button(master=frame_circulos, image=boton, command=None,relief="flat",height=28, width=28).place(x=240,y=1)
L_t=tk.Button(master=frame_circulos, image=boton, command=None,relief="flat",height=28, width=28).place(x=320,y=1)

#crea el label que dara nombre a cada boton circular
label_S=tk.Label(master=frame_circulos, text="S(t)",bg="white",
                   font=('Arial', 12, 'bold italic')).place(x=80,y=35)
label_E=tk.Label(master=frame_circulos, text="E(t)",bg="white",
                   font=('Arial', 12, 'bold italic')).place(x=160,y=35)
label_I=tk.Label(master=frame_circulos, text="I(t)",bg="white",
                   font=('Arial', 12, 'bold italic')).place(x=240,y=35)
label_L=tk.Label(master=frame_circulos, text="L(t)",bg="white",
                   font=('Arial', 12, 'bold italic')).place(x=320,y=35)

#imagen de sumplemento para la grafica
imgen=ImageTk.PhotoImage(Image.open("Ima.PNG")) #importa la imagen
lab=tk.Label(image=imgen) #crea el label
lab.place(x=100,y=100) #ubica la imagen

#Titulos para cada frame

titulo_metodos = tk.Label(master=frame_metodos,bg="white",
                   font=('Arial', 12,'bold italic'), text=f"Métodos de solución",height=-8).grid(row=0, column=0,pady=5)

titulo_parametros = tk.Label(master=frame_parametros, bg="#77cfe7",
                   font=('Arial', 12,'bold italic'), text=f"Parámetros",height=2).grid(row=1,pady=3)
titulo_tiempo=tk.Label(window, bg="white",font=('Arial', 12,'bold italic'), text=f"Tiempo de simulación",height=2).place(x=200,y=440)
#titulo_tiempos= tk.Label(master=frame_tiempo, bg="white",font=('Arial', 12, 'bold italic'), text=f"Tiempo de simulacion",height=8).place(x=1, y=1)
window.mainloop()
