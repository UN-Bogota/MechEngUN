#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 14:34:22 2023

@author: bojack
"""

import numpy as np
from kineticsEquationsObject import *
from resources import *
from scikits.odes import ode
#from scipy.integrate import ode
import matplotlib.pyplot as plt



############ Mecanismo simplificado CH4 (WD) ##################

# Mecanismo: 
    # CH4 + 3/2 O2 ----> CO + 2H2O
    # CO + 0.5 O2 -------------> CO2
    # CO2 ----> CO + 0.5O2


# Matriz de coeficientes estequeometericos

                #CH4 O2 CO H2O CO2
nup = np.array([[1, 3/2, 0, 0, 0],
                [0, 0.5, 1, 0, 0],
                [0, 0, 0, 0, 1]
               ])

                #CH4 O2 CO H2O CO2
nupp = np.array([[0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0]
               ])



T = 1500 #K

# la presion siempre en pascales

P = 101325 #Pa
AR = 1/0.42
phi = 1
# Ru = 1,987 # [cal/mol.K]
R = 8.31446261815324 # [J/mol.K]

lab = ['CH4',' O2', 'CO', 'H2O', 'CO2']

##### Construir las ecuaciones #################

ecuaciones = kinetics(nup, nupp)
ecuaciones.setT(T)


# Condiciones iniciales

y0=np.array([1, 2, 0, 0, 0])

#y0=y0/sum(y0)


v_molar_CH4 = R*T/P # [m³/mol CH4]

v_molar_CH4 = v_molar_CH4*1000 # [L/mol CH4]
y0 = y0*(1/v_molar_CH4)

t0=0
tf=1
dt=1e-6


solution = ode('cvode', ecuaciones.getDiffEq, old_api=False).solve(
    np.arange(t0,tf,dt), y0)


V= P/(R*T)
#sol=solution.values.y[:,]
sol=np.divide(solution.values.y[:,0:5],np.sum(solution.values.y,axis=1).reshape(len(solution.values.t),1))
plt.figure(1)
plt.plot(solution.values.t, solution.values.y*V)
plt.legend(lab,loc='upper right')
plt.title('$\phi$ ='+str(phi)+', T = '+str(T)+' [K], P = '+str(P/101325)+' [atm]')
plt.ylabel('Concentración [kmol/$m^3$]')
plt.xlabel('tiempo [s]')

plt.figure(2)
plt.plot(solution.values.t[:4], sol[:4])
plt.legend(lab,loc='upper right')
plt.title('$\phi$ ='+str(phi)+', T = '+str(T)+' [K], P = '+str(P/101.325)+' [atm]')
plt.ylabel('Fracción molar')
plt.xlabel('tiempo [s]')
#plt.plot(solution.values.t, solution.values.y[:,], label='Kinetishe')

