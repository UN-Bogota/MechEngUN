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
import pandas as pd

sys.setrecursionlimit(4000)


def genInitCon(solveWithEnergy = False):
    
    n = random.random()*10
    
    ###########################################################################
    #['CO2', 'CO', 'H2O', 'H2', 'OH', 'O2', 'N2', 'H', 'O']#
    ###########################################################################
    
    if solveWithEnergy:
        
        if CH4Reac.phi < 1.01:
            CI = [0.56, 0.001, 0.4, 0.001*n, 0.0001*n, 0.0001*n, 0.67, 7.74325e-05*n, 5.278e-05*n, 0.0001, 1, 0.0001, 0.001, 2200]
        
        else:
            CI = [0.2*n, 0.01*n, 0.3*n, 0.00703*n, 0.0001*n, 0.0001*n, 0.7, 7.74325e-05*n, 5.278e-05*n, 0.0001, 0.0001, 0.001, 2200]
    else:
        
        if CH4Reac.phi < 1.01:
            CI = [0.056, 0.0001*n, 0.2*n, 0.0001*n, 0.001, 0.01*n, 0.7, 0.0001*n, 0.01*n]
            
        else:
            CI = [1.0*n, 1e-06*n, 1.923*n, 1e-06, 0.0, 0.0, 7.38*n, 1e-05, 0.1*n]
            
    return CI

def molarProdFrac(products):
    
    molarProdFrac = []
    N_tot = sum(products)
    
    for i in products:
        molarProdFrac.append(i/N_tot)
    return molarProdFrac

def findSolution(reaction, CI, epsilon = 1E-10):
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
    deltaMass_ant = 1
    noNegativity = False
    
    sol = None
    while abs(deltaMass) >= epsilon or (not noNegativity):
        
        deltaMass_ant = deltaMass
        sol_i = reaction.solveSystem(CI)
        deltaMass = reaction.getReactmass() - reaction.getProdMass()
        
        count = True
        
        for i in sol_i:
            if i <= 0:
                count = False
                break
        if not count:
            n = random.random()*10
            CI = genInitCon()
        if count:
            sol = sol_i
            noNegativity = True
        if deltaMass == deltaMass_ant:
            CI = genInitCon()
    if sol != None:    
        return sol
    else:
        print('se va a putiar')
        return findSolution(reaction, CI)

# Creo el objeto
CH4Reac = Reaction()

# Se definen las especies del combustible
methane = [1, 4, 0, 0]
CH4_name = 'CH4'



# CONDIOCIONES CONOCIDAS

presion = 101.325*2
temperatura = 700 # K


#----------------------------------------------
CH4Reac.addFuelSpecies(methane, CH4_name)
CH4Reac.addProductPressure(presion)
CH4Reac.addProductTemperature(temperatura)
CH4Reac.addFirstLaw(False)


#----------------------------------------------

phi = np.arange(0.05, 1.7, 0.01)
phi = np.round(phi, 2)

#------------------------------------------


resul= np.zeros([len(phi),len(CH4Reac.productNames)])
adiaTemps =np.zeros([len(phi),1])
CI = genInitCon()
j=0

for i in phi:
    try:
        print('para el phi: ', i, ' --------------------------------')
        CH4Reac.addPhi(i)
        CH4Reac.getTotalReacRealElements()
        sol = findSolution(CH4Reac, CI)
        if CH4Reac.solveWithEnergy:
            resul[j,:]=sol[0:-1]/sum(sol[0:-1])
            adiaTemps[j,:] = sol[-1]    
        else:
            resul[j,:]=sol/sum(sol)
        CI = sol
        print(sol)
        print(CH4Reac.get_h_prod(), 'kJ/kmol CH4')
        j+=1
    except KeyboardInterrupt:
        result = resul
        titles = ['$phi$','CO2', 'CO', 'H2O', 'H2', 'OH', 'O2', 'N2', 'H', 'O']
        
        result = np.hstack((phi.reshape(len(phi),1), result))
        result = np.vstack((titles, result))
        df = pd.DataFrame(result)
        if CH4Reac.solveWithEnergy:
            filename = 'equlibrium_result_CH4_at_P_' + str(presion) + '_T_ad.csv'
        else:
            filename = 'equlibrium_result_CH4_at_P_' + str(presion) + '_T_'+ str(temperatura) + '_' + str(phi[0]) + '_phi_'+ str(i-0.01) + '.csv'
            df.to_csv(filename, index=False)
        break

    
if CH4Reac.solveWithEnergy:
    fig, axs = plt.subplots(2)
else:
    fig, axs = plt.subplots(1)
    
    
for k in range(len(CH4Reac.productNames)):
    
    if CH4Reac.solveWithEnergy:
        axs[0].plot(phi, resul[:,k], label=CH4Reac.productNames[k])
    else:
        axs.plot(phi, resul[:,k], label=CH4Reac.productNames[k])
   

if CH4Reac.solveWithEnergy:
    axs[0].grid(color='k', linestyle=':')
    axs[0].set(xlabel='Radio de equivalencia en la entrada $\phi$', ylabel='Fracción molar')
    axs[1].grid(color='k', linestyle=':')
    axs[1].plot(phi, adiaTemps, color='k')
    axs[1].set(xlabel='Radio de equivalencia en la entrada $\phi$', ylabel='$T_{adiabática}$')
    fig.suptitle('Fracción molar en equilibrio a P = ' + str(presion) + ' kPa y temperatura de llama adiabática')
    
else:
    axs.grid(color='k', linestyle=':')
    axs.set(xlabel='Radio de equivalencia en la entrada $\phi$', ylabel='Fracción molar')

    fig.suptitle('Fracción molar en equilibrio a P = ' + str(presion) + ' y temperatura = ' + str(temperatura))

fig.legend(loc='outside right center')

# Hide x labels and tick labels for top plots and y ticks for right plots.
if CH4Reac.solveWithEnergy:
    for ax in axs.flat:
        ax.label_outer()
    
plt.show()

#-------------------------------------------------------------------------------------------------------------------

if CH4Reac.solveWithEnergy:
    titles = ['$phi$', '$T_ad$', 'CO2', 'CO', 'H2O', 'H2', 'OH', 'O2', 'N2', 'H', 'O']
    result = np.hstack((adiaTemps, resul))
else:
    result = resul
    titles = ['$phi$','CO2', 'CO', 'H2O', 'H2', 'OH', 'O2', 'N2', 'H', 'O']

result = np.hstack((phi.reshape(len(phi),1), result))
result = np.vstack((titles, result))
df = pd.DataFrame(result)
if CH4Reac.solveWithEnergy:
    filename = 'equlibrium_result_CH4_at_P_' + str(presion) + '_T_ad.csv'
else:
    filename = 'equlibrium_result_CH4_at_P_' + str(presion) + '_T_'+ str(temperatura) +'.csv'
df.to_csv(filename, index=False)