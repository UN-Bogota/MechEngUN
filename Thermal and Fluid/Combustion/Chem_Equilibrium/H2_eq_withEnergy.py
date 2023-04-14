#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

sys.setrecursionlimit(2000)


def genInitCon(solveWithEnergy = False):
    
    #random.seed(10, CO H2O, H2, OH, O2, N2, NO, C3H8, H, O, N)
    n = random.random()*10
    
    if solveWithEnergy:
        CI = [0.01*n, 0.1, 0.01*n, 0.002, 0.01*n, 0.7, 0.3*n, 0.01*n,
              2, 2000]
    
    else:
        CI = [0.01*n, 0.1, 0.01*n, 0.002, 0.01*n, 0.7, 0.3*n, 0.01*n,
              2]
        #CI = [1.745e-09, 0.20, 1.0195e-16, 4.714e-19, 1.37e-08, 0.576, 2.51726e-05, 0.20644, 1.8]
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


#phi = 1
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

H2Reac.addComp_salida_est(products_est)

H2Reac.addComp_salida_real(products_real)

H2Reac.addFuelSpecies(hydrogen, H2_name)

H2Reac.addProductPressure(presion)
H2Reac.addProductSpecies(productNames)
H2Reac.addFirstLaw(True)

#sol = findSolution(H2Reac)

#print(sol/sum(sol))
phi = np.linspace(0.5, 1.7, 26)

resul= np.zeros([len(phi),9])
adiaTemps =np.zeros([len(phi),1])
j=0
for i in phi:
    print('para el phi: ', i, ' --------------------------------')
    H2Reac.addPhi(i)
    sol = findSolution(H2Reac)
    print(sol)
    input()
    resul[j,:]=sol[0:-1]/sum(sol[0:-1])
    adiaTemps[j,:] = sol[-1]
    j+=1
    
plt.figure(1)
for k in range(9):
    plt.plot(phi,resul[:,k],label=productNames[k])
plt.legend()
plt.title('Fracción molar en equilibrio a P= ' + str(presion) + ' kPa y T= ' + str(temperatura) + ' K')
plt.ylabel('Fracción molar [mol/mol H2]')
plt.xlabel('Radio de equivalencia en la entrada $\phi$')
plt.grid(color='k', linestyle=':')
plt.show()

