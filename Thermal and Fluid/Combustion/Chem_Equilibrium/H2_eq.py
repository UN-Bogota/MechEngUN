#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 11:23:50 2023

@author: HIDES
"""

import random
import numpy as np
import matplotlib.pyplot as plt
from reactionObject import Reaction
import sys
import time

sys.setrecursionlimit(2000)

inicio = time.time()
def genInitCon(solveWithEnergy = False):
    
    #random.seed(10, CO H2O, H2, OH, O2, N2, NO, C3H8, H, O, N)
    n = random.random()*10
    
    if solveWithEnergy:
        CI = [0.4*n, 0.01*n, 0.45*n, 0.1*n, 0.02*n, 0.1*n, 0.1*n, 0.1*n,
              0.8*n, 2000]
    
    else:
        #CI = [2.870e-05*n, 0.1, 0.00062*n, 0.2*n, 4.42e-4*n, 0.86, 0.000001*n, 0.00001*n, 1.8809]
        if H2Reac.phi < 1.02:
            CI =[0.002717*n, 0.1, 0.00367782*n, 0.2327979*n, 0.000285*n, 0.9828002, 0.031395*n, 5.306e-12*n, 2.849]
            #CI = [9.84e-04*n, 0.00156, 0.0003276*n, 0.28049*n, 0.00387*n, 0.9961, 0.00056*n, 2.26035517e-05*n, 2.93]
        else:
            #CI = [2.870e-05*n, 0.1, 0.00062*n, 0.2*n, 4.42e-4*n, 0.86, 0.000001*n, 0.00001*n, 1.8809]
            CI = [2.6053e-06*n, 0.2308, 3.716e-05*n, 1.67e-07*n, 6.3883e-4*n, 0.76915, 2.51152e-06*n, 1.5388e-05*n, 1.4468]
            #CI = [0.001147*n, 0.19967223, 1.92052e-06*n, 3.3413e-04*n, 0.0004743*n, 0.7995, 1.608e-05, 1.76e-05, 1.50]
    return CI

def molarProdFrac(products):
    
    molarProdFrac = []
    N_tot = sum(products)
    
    for i in products:
        molarProdFrac.append(i/N_tot)
    return molarProdFrac

def findSolution(reaction, epsilon = 1E-9):
    """
     Function thats verify both the positive of all concentrations
     and the preservation of mass conservation.

    Parameters
    ----------
    reaction : Object 
        DESCRIPTION.
    epsilon : value of convergence
        DESCRIPTION. The default is 1E-6.

    Returns
    -------
    TYPE
        The concentrations at equilibrium state

    """
    
    
    deltaMass = 1
    noNegativity = False
    
    sol = None
    
    while (not noNegativity) and  abs(deltaMass) >= epsilon :

        CI = genInitCon(reaction.solveWithEnergy)
        
        sol_i = reaction.solveH2System(CI)
        deltaMass = reaction.getReactmass() - reaction.getProdMass()
        
        count = True
        
        for i in sol_i:
            if i <= 0:
                count = False
                break
            
        if count:
            sol = sol_i
            noNegativity = True
            
    if sol != None:    
        return sol
    else:
        return findSolution(reaction)

# Creo el objeto
H2Reac = Reaction()

# Se definen las especies del combustible

hydrogen = [0, 2, 0, 0]
H2_name = 'H2'

temperatura = 2000
presion = 101.325*10

productNames = ['H', 'H2', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2', 'N2']

products_real = np.array([
                    [0, 1, 0, 0], # H
                    [0, 2, 0, 0], # H2
                    [0, 0, 1, 0], # O
                    [0, 0, 2, 0], # O2
                    [0, 1, 1, 0], # OH
                    [0, 2, 1, 0], # H2O 
                    [0, 1, 2, 0], # HO2 
                    [0, 2, 2, 0], # H2O2 
                    [0, 0, 0, 2], # N2
])

products_est = np.array([
                    [0, 2, 1, 0], # H2O
                    [0, 0, 0, 2] # N2
])

#-------------------------------------------------------------

H2Reac.addComp_salida_est(products_est)
H2Reac.addComp_salida_real(products_real, productNames)
H2Reac.addFuelSpecies(hydrogen, H2_name)
H2Reac.addProductTemperature(temperatura)
H2Reac.addProductPressure(presion)
H2Reac.addFirstLaw(False)

#-------------------------------------------------------------

phi = np.arange(0.5, 1.7, 0.01)
phi = np.round(phi, 2)

#------------------------------------------------------------


resul= np.zeros([len(phi),len(productNames)])
j=0
for i in phi:
    print('para el phi: ', i, ' --------------')
    H2Reac.addPhi(i)
    sol = findSolution(H2Reac)
    print(sol)
    resul[j,:]=sol/sum(sol)
    j+=1
    
plt.figure(1)

for k in range(len(productNames)):
    plt.plot(phi,resul[:,k],label=productNames[k])
    
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('Fracción molar en equilibrio a P= ' + str(presion) + ' kPa y T= ' + str(temperatura) + ' K')
plt.ylabel('Fracción molar [mol/mol H2]')
plt.xlabel('Radio de equivalencia en la entrada $\phi$')
plt.grid(color='k', linestyle=':')
plt.show()

fin = time.time()

print(fin-inicio)

titles = ['$phi$', 'H', 'H2', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2', 'N2']

result = np.hstack((phi.reshape(len(phi),1), resul))
result = np.vstack((titles, result))
import pandas as pd

df = pd.DataFrame(result)
filename = 'equlibrium_result_at_T_' + str(temperatura) + '_P_' + str(presion)+'.csv'
df.to_csv(filename, index=False)

