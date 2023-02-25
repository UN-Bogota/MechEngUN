#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 12:41:53 2023

@author: bojack
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
from scipy.optimize import fsolve
from sympy import *


# Coeficientes estequeometricos de los componentes de entrada

n_entrada = np.array([0.7478,
                      0.1099,
                      0.0512,
                      0.0223,
                      0.0076,
                      0.0027,
                      0.0022,
                      0.0506,
                      0.0057])

# C H O N  -----> Elementos reactivos 
comp_entrada = np.array([[1, 4, 0, 0], 
                        [2, 6, 0, 0],
                        [3, 8, 0, 0],
                        [4, 10, 0, 0],
                        [5, 12, 0, 0],
                        [6, 14, 0, 0],
                        [7, 16, 0, 0],
                        [1, 0, 2, 0],
                        [0, 0, 0, 2]])


# Coeficientes estequeometricos de la reacción de salida real (se deja vacio)
n_salida_real = np.array([])
# Coeficientes estequeometricos de la reacción de salida ideal (se deja vacio)
n_salida_est = np.array([])

# C H O N  -----> Elementos que componen los productos (reac. estequeometrica)
comp_salida_est = np.array([[1, 0, 2, 0],
                            [0, 2, 1, 0],
                            [0, 0, 0, 2]])

# C H O N  -----> Elementos que componen los productos reales
comp_salida_real = np.array([[1, 0, 2, 0],
                             [1, 0, 1, 0],
                             [0, 2, 1, 0],
                             [0, 2, 0, 0],
                             [0, 1, 1, 0],
                             [0, 0, 2, 0],
                             [0, 0, 0, 2],
                             [0, 0, 1, 1]])

    
def getTotalFuelElements(n_entrada, comp_entrada):
    Carbon = 0
    Hydrogen = 0
    Oxygen = 0
    Nitrogen = 0
    
    for i in range(len(comp_entrada)):
        Carbon += n_entrada[i] * comp_entrada[i][0]
        Hydrogen += n_entrada[i] * comp_entrada[i][1]
        Oxygen += n_entrada[i] * comp_entrada[i][2]
        Nitrogen += n_entrada[i] * comp_entrada[i][3]
        
    return Carbon, Hydrogen, Oxygen, Nitrogen

def getA_est(Oxygen, Carbon, Hydrogen):
    
    A_est = ((Carbon*2 + Hydrogen/2) - Oxygen)/(2*0.21) # Balance de Oxigenos
    
    return A_est

def getA_real(A_est):
    
    A_real = phi*A_est
    
    return A_real

Carbon, Hydrogen, Oxygen, Nitrogen = getTotalFuelElements(n_entrada, comp_entrada)
phi = 1.2
A_est = getA_est(Oxygen, Carbon, Hydrogen)
A_real = getA_real(A_est)

N_total = Nitrogen + A_real*0.79*2; # Nitrógeno total

n_salida_est = np.append(n_salida_est, np.array([Carbon, Hydrogen/2, N_total, A_est]))




C_total = Carbon # Carbono total
O_real = A_real*2*0.21 + Oxygen #Oxígeno total
H_total = Hydrogen # Hidrogeno total 


#------------------------------------------------------------

presion = np.arange(50, 300, 50)  # [kPa]

T_produc = np.arange(1000, 3000, 500)

#Cantidad total de O:



#print(O_real)

# B CO2 + c CO + + D H2O + E H2 + F OH + G O2 + h NO + I N2

#Balance de O : Ot = 2B + c + D + F

#Ecuaciones de equilibrio

def k1(T): # H2O--------- OH + (1/2) H2
    n= -11.280*(T == 1000) - 6.344 * (T == 1500) -3.776*(T==2000)- 2.27*(T == 2500)
    return 10**n

def k2(T): # CO2 + H2 ------- CO + H2O
    n=-0.159*(T==1000)+0.4035*(T==1500)+0.656*(T==2000)+0.784*(T==2500)
    return 10**n

def k3(T): # (1/2) * O2 + (1/2) * N2 -------- NO
    n=-4.062*(T==1000)-2.501*(T==1500)-1.699*(T==2000)-1.227*(T==2500)
    return 10**n

def k4(T): # CO2 -------- CO + (1/2) * O2
    n=-10.221*(T==1000)-5.36*(T==1500)-2.884*(T==2000)-1.44*(T==2500)
    return 10**n

# def equations(p):
#     x, y, z = p
#     eqn1 = x+z-1
#     eqn2 = x+2*y+2*z-2
#     eqn3 = x * y**(1/2) / (z) * (1 / (x+y+z)**(1/2)) - 0.0363
#     return (eqn1, eqn2,eqn3)

# x, y ,z =  fsolve(equations, (0.5, 1, 1))

P=50 #kPa
T=1000 #K 

def equations(vars):
    """
    Se usan las 4 ecuaciones de balance de masa y 4 de equilibrio para 
    encontrar cada coeficiente, los coeficientes corresponden a:
                
        B: CO2
        C: CO
        D: H2O
        E: H2
        F: OH
        G: O2
        H: N2
        I: NO
        J: C3H8
        K: H
        L: O
        M: NO2
        falta el N
    """
    
    B, C, D, E, F, G, H, I, J, K, L, M = vars
    n = B + C + D + E + F + G + H + I + J + K + L + M 
    
    f = ((P/101.325)/n)**(1/2)
    
    # Todas se igualan a 0
    
    C_balance = B + C + 3*J - C_total
    H_balance = 2*D + 2*E + F + 8*J + K - H_total
    O_balance = 2*B + C + D + F + 2*G + I + L + 2*M - O_real 
    N_balance = 2*H + I + M - N_total
    
    eqn5 = F**2 * E * f - (k1(T) * D)**2 # # H2O ----- OH+1/2H2
    #eqn5 = ((F * E**(1/2)) / D) * f-k1(T)
    eqn6 = c * D - k2(T) * (B * E) # CO2 + H2 ----- CO + H2O
    eqn7 = h**2 - k3(T)**2 * (G * I) # 1/2 * O2 + 1/2 * N2---------NO
    #eqn7 = h / (G**(1/2) * I **(1/2)) - k3(T)
    eqn8 = c**2 * G * f - (k4(T) * B)**2 # CO2 ------ CO + 1/2 * O2
    #eqn8 = c * G**(1/2) / B * f - k4(T)
    
    
    # eqn5=F1**2*e1-(k1(T))**2*(n/f)*D1**2 # #H2O-OH+1/2H2
    # eqn6=c1*D1-k2(T)*B1*e1*n #CO2+H2-CO+H2O
    #eqn7=h1**2-(k3(T))**2*G1*I1*n**2 
    #eqn8=c1**2*G1-(k4(T))**2*(n/f)*B1**2 
    
    return [C_balance, H_balance, O_balance, N_balance, 
            eqn5, eqn6, eqn7, eqn8]

CI = [1,1,2,1,0.1,0.1,1,8]
B, c, D, E, F, G, h, I =  fsolve(equations, CI)


# import warnings
# ivar = True
# while ivar: 

#     try: 
#         warnings.simplefilter("error", category=RuntimeWarning)
        
#     except : 
#         warnings.resetwarnings()
        
#         for i in range(len(CI)):
#             if CI[i] >= 10:
#                 print('entro')
#                 ivar = False
#                 break
#             CI[i] *= 1.1
#         print(CI)
        
#     else: 
#         ivar = False
            
print(B, c, D, E, F, G, h, I)








