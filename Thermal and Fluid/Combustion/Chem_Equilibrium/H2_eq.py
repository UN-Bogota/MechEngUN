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

sys.setrecursionlimit(3000)

inicio = time.time()
def genInitCon(solveWithEnergy = False):
    
    #['H', 'H2', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2', 'N2']
    n = random.random()*10
    
    if solveWithEnergy:
        CI = [0.4*n, 0.01*n, 0.45*n, 0.1*n, 0.02*n, 0.1*n, 0.1*n, 0.1*n,
              0.8*n, 2000]
    
    else:
        #CI = [2.870e-05*n, 0.1, 0.00062*n, 0.2*n, 4.42e-4*n, 0.86, 0.000001*n, 0.00001*n, 1.8809]
        if H2Reac.phi < 1.01:
            
            n = random.uniform(0.1, 0.7)
            n = round(n, 4)
            #CI = [0.0001, 0.002*n, 0.0001, 0.5*n, 0.006, 1.5*n, 0.0000001, 0.0000001, 3*n] #->1_2000
            #CI = [0.0000001, 0.000003*n, 0.00001, 0.3*n, 0.0001, 1.5*n, 0.0000001, 0.000001, 3*n] #->1_1500
            
            #CI = [0.00013*n, 0.01, 0.000703*n, 0.2121*n, 0.0102*n, 0.9934, 7.74325e-05*n, 5.278e-05*n, 2.687]
        elif H2Reac.phi >= 1.01 and H2Reac.phi < 1.15:
            n = random.uniform(0.4, 1)
            n = round(n, 4)
            CI = [0.01, 0.1*n, 0.00001, 0.000001, 0.0001, 0.8*n, 0.000001, 0.0000001, 1.5*n]
        else:
            n = random.uniform(0.1, 1)
            n = round(n, 4)
            CI = [0.01, 0.3*n, 0.00001, 0.000001, 0.0001, 0.8*n, 0.000001, 0.0000001, 1.3*n]
            
    return CI

def molarProdFrac(products):
    
    molarProdFrac = []
    N_tot = sum(products)
    
    for i in products:
        molarProdFrac.append(i/N_tot)
    return molarProdFrac

def findSolution(reaction, epsilon = 1E-10):
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

temperatura = 1500
presion = 101.325

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

