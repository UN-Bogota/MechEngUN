#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 14:34:22 2023

@author: bojack
"""

import numpy as np
from equation import *
from resources import *
from scikits.odes import ode
#from scipy.integrate import ode
import matplotlib.pyplot as plt

                #H H2 O O2 OH H2O N2 HO2 H2O2 M5 M6 M7 M8 M9 M15
nup = np.array([[1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],
                [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]              
               ])
                 #H H2 O O2 OH H2O N2 HO2 H2O2 M5 M6 M7 M8 M9 M15 
nupp = np.array([[0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                 [0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                 [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                 [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0],
                 [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],                                        
                 [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                 [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0]                                       
                ])

T=1500 #K
P=101.325 #Pa
AR=1/0.42
phi=1
R=8.3144 #J/(mol*K)
lab=['H2O2',' H2', 'O', 'O2', 'OH' ,'H2O', 'N2', 'HO2', 'H']
ecuaciones = kinetics(nup, nupp)
ef = get_ef()
ecuaciones.setEfficiency(ef)
ecuaciones.setT(T)
y0=np.array([0,1,0,0.21*AR/phi,0,0,0.79*AR/phi,0, 0])
y0=y0/sum(y0)
y0=y0*((P/101.325)/(T/1000))
t0=0
tf=1e-4
dt=1e-10
#Parametrizar con phi y temperatura
# r = ode(ecuaciones.getDiffEq).set_integrator('vode', method='bdf')
# r.set_initial_value(y0, t0)
# t=[]
# sol=[]
# while r.successful() and r.t < 5*dt:
#     t.append(r.t+dt)
#     sol.append(r.integrate(r.t+dt))
#     print(sol)
#print(k_values(1000))
#print(ecuaciones.getDiffEq(0,y0,np.zeros(9)))
solution = ode('cvode', ecuaciones.getDiffEq, old_api=False).solve(np.arange(t0,tf,dt), y0)
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

