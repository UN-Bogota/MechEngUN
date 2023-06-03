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
nup = np.array([[0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0]
               ])



T = 1500 #K
P = 101.325 #Pa
AR = 1/0.42
phi = 1
R = 8.3144 #J/(mol*K)

lab = ['CH4',' O2', 'CO', 'H2O', 'CO2']

##### Construir las ecuaciones #################

ecuaciones = kinetics(nup, nupp)
ecuaciones.setT(T)

y0=np.array([0,1,0,0.21*AR/phi,0,0,0.79*AR/phi,0, 0])
y0=y0/sum(y0)
y0=y0*((P/101.325)/(T/1000))
t0=0
tf=1e-4
dt=1e-10


solution = ode('cvode', ecuaciones.getDiffEq, old_api=False).solve(
    np.arange(t0,tf,dt), y0)


V= P/(R*T)
#sol=solution.values.y[:,]
sol=np.divide(solution.values.y[:,0:9],np.sum(solution.values.y,axis=1).reshape(len(solution.values.t),1))
plt.figure(1)
plt.plot(solution.values.t, solution.values.y*V)
plt.legend(lab,loc='upper right')
plt.title('$\phi$ ='+str(phi)+', T = '+str(T)+' [K], P = '+str(P/101325)+' [atm]')
plt.ylabel('Concentración [kmol/$m^3$]')
plt.xlabel('tiempo [s]')

plt.figure(2)
plt.plot(solution.values.t, sol)
plt.legend(lab,loc='upper right')
plt.title('$\phi$ ='+str(phi)+', T = '+str(T)+' [K], P = '+str(P/101.325)+' [atm]')
plt.ylabel('Fracción molar')
plt.xlabel('tiempo [s]')
#plt.plot(solution.values.t, solution.values.y[:,], label='Kinetishe')

