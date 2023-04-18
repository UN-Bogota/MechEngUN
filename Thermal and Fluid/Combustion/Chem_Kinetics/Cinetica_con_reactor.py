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
#from scipy.integrate import solve_ivp
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


T=1000 #K
P=101325 #Pa
AR=1/0.42
phi=1
R=8.3144 #J/(mol*K)
lab=['H2O2',' H2', 'O', 'O2', 'OH' ,'H2O', 'N2', 'HO2', 'H']

y0=np.array([0,1,0,0.21*AR/phi,0,0,0.79*AR/phi,0, 0])
y0=y0/sum(y0)
y0=y0*P/(R*T)
y0=np.append(y0,[T,P])

t0=0
tf=5e-6
dt=1e-12
ecuaciones = kinetics(nup, nupp)
ef = get_ef()
ecuaciones.setEfficiency(ef)
ecuaciones.setReactorType('constant-volume')
#print(ecuaciones.getDiffEq(0,y0,np.zeros(10)))

solution = ode('cvode', ecuaciones.getDiffEq, old_api=False).solve(np.arange(t0,tf,dt), y0)
V= solution.values.y[:,10]/(R*solution.values.y[:,9])
sol=np.divide(solution.values.y[:,0:9],V.reshape(len(V),1))
#plt.subplot(3,1,1)

plt.figure(1)
plt.plot(solution.values.t, solution.values.y[:,0:9])
plt.legend(lab,loc='upper right')
plt.title('Cinética en reactor para $\phi$ ='+str(phi)+', T = '+str(T)+' [K], P = '+str(P/101325)+' [atm]')
plt.ylabel('concentración [kmol/$m^3$]')
plt.xlabel('tiempo [s]')
#plt.xlim([2.82e-8,2.86e-8])
#plt.subplot(3,1,2)

plt.figure(2)
plt.plot(solution.values.t, sol)
plt.legend(lab,loc='upper right')
plt.title('Cinética en reactor para $\phi$ ='+str(phi)+', T = '+str(T)+' [K], P = '+str(P/101325)+' [atm]')
plt.ylabel('Fracción molar')
plt.xlabel('tiempo [s]')

plt.figure(3)
plt.plot(solution.values.t, solution.values.y[:,9])
plt.title('Temperatura en reactor para $\phi$ ='+str(phi)+', T = '+str(T)+' [K], P = '+str(P/101325)+' [atm]')
plt.ylabel('Temperatura [K]')
plt.xlabel('tiempo [s]')
#plt.subplot(3,1,3)
plt.figure(4)
plt.plot(solution.values.t, solution.values.y[:,10]/101325)
plt.title('Presión en reactor para $\phi$ ='+str(phi)+', T = '+str(T)+' [K], P = '+str(P/101325)+' [atm]')
plt.ylabel('Presión [atm]')
plt.xlabel('tiempo [s]')
#plt.plot(solution.values.t, solution.values.y[:,], label='Kinetishe')
