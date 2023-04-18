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
        
        if H2Reac.phi < 0.9:
            #CI =[0.002717*n, 0.1, 0.00367782*n, 0.2327979*n, 0.000285*n, 0.9828002, 0.031395*n, 5.306e-12*n, 2.849, 2000]
            CI = [9.84e-04*n, 0.00156, 0.003276*n, 0.28049*n, 0.00387*n, 0.9961, 0.00056*n, 2.26035517e-05*n, 2.93, 1900]
            #CI = [0.0001197*n, 0.00479, 0.000215*n, 0.078*n, 1.72e-04*n, 0.9, 0.019*n, 0.00034*n, 2.5, 1900]
        else:
            #CI = [2.870e-05*n, 0.1, 0.00062*n, 0.2*n, 4.42e-4*n, 0.86, 0.000001*n, 0.00001*n, 1.8809]
            #CI = [2.6053e-06*n, 0.2308, 3.716e-05*n, 1.67e-07*n, 6.3883e-4*n, 0.76915, 2.51152e-06*n, 1.5388e-05*n, 1.4468, 2200]
            #CI = [0.001147*n, 0.19967223, 1.92052e-06*n, 3.3413e-04*n, 0.0004743*n, 0.7995, 1.608e-05, 1.76e-05, 1.50, 2200]
            CI = [0.0057*n, 0.2, 1.724e-04*n, 1.33e-03*n, 0.0015*n, 0.8, 1.67e-05*n, 1.0551e-05*n, 1.5, 2200]    
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

#--------------------------------------------------
# Creo el objeto
H2Reac = Reaction()
# Se definen las especies del combustible
hydrogen = [0, 2, 0, 0]
H2_name = 'H2'

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

#----------------------------------------------
H2Reac.addComp_salida_est(products_est)
H2Reac.addComp_salida_real(products_real, productNames)
H2Reac.addFuelSpecies(hydrogen, H2_name)
H2Reac.addProductPressure(presion)
H2Reac.addFirstLaw(True)

#----------------------------------------------

phi = np.arange(0.5, 1.7, 0.01)
phi = np.round(phi, 2)

#------------------------------------------


resul= np.zeros([len(phi),9])
adiaTemps =np.zeros([len(phi),1])
j=0

for i in phi:
    print('para el phi: ', i, ' --------------------------------')
    H2Reac.addPhi(i)
    sol = findSolution(H2Reac)
    print(sol)
    resul[j,:]=sol[0:-1]/sum(sol[0:-1])
    adiaTemps[j,:] = sol[-1]
    j+=1

    
fig, axs = plt.subplots(2)

for k in range(9):
    axs[0].plot(phi, resul[:,k], label=productNames[k])
    
axs[0].grid(color='k', linestyle=':')
axs[1].grid(color='k', linestyle=':')

axs[1].plot(phi, adiaTemps, color='k')

fig.legend(loc='outside right center')

axs[0].set(xlabel='Radio de equivalencia en la entrada $\phi$', ylabel='Fracción molar')
axs[1].set(xlabel='Radio de equivalencia en la entrada $\phi$', ylabel='$T_{adiabática}$')


# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    

fig.suptitle('Fracción molar en equilibrio a P = ' + str(presion) + ' kPa y temperatura de llama adiabática')

plt.show()

titles = ['$phi$', 'H', 'H2', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2', 'N2', '$T_ad$']

result = np.hstack((phi.reshape(len(phi),1), resul))
result = np.vstack((titles, result))

import pandas as pd

df = pd.DataFrame(result)
filename = 'equlibrium_result_at_P_' + str(presion) + '_T_ad.csv'
df.to_csv(filename, index=False)