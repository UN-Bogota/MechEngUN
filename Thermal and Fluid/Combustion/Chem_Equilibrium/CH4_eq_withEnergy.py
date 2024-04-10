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


def genInitCon():
    
    n = random.random()*10
    
    ###########################################################################
    #['CO2', 'CO', 'H2O', 'H2', 'OH', 'O2', 'N2', 'H', 'O']#
    ###########################################################################
    
    if CH4Reac.solveWithEnergy:
        
        if CH4Reac.phi < 0.5:
            CI = [0.9999*n, 0.0001*n, 0.2*n, 0.0001, 0.001, 0.01*n, 7.38*n, 0.0001*n, 0.01*n, 800]
        
        elif CH4Reac.phi < 0.8:
            CI = [0.889*n, 0.0001*n, 0.2*n, 0.0001*n, 0.001, 0.01*n, 7.38*n, 0.0001*n, 0.01*n, 1500]
        
        else:
            CI = [0.3*n, 0.1*n, 1.923*n, 1e-06, 0.0001, 0.0001, 7.38*n, 1e-05, 0.1*n, 2000]
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

def findSolution(reaction, CI, epsilon = 1E-10, args = []):
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
            CI = genInitCon()
        if count:
            sol = sol_i
            noNegativity = True
        if deltaMass == deltaMass_ant:
            CI = genInitCon()
    
        
    if sol != None: 
        if CH4Reac.solveWithEnergy:
            if sol[-1] < args[0] or sol[-1]> 3000:
                print('se va a putiar: T= ',  args[0])
                
                return findSolution(reaction, CI, args=[args[0]-1])
        return sol
    else:
        print('se va a putiar')
        return findSolution(reaction, CI, args=args[0])

# Creo el objeto
CH4Reac = Reaction()

# Se definen las especies del combustible
methane = [1, 4, 0, 0]
CH4_name = 'CH4'


presion = 123.79

#----------------------------------------------
CH4Reac.addFuelSpecies(methane, CH4_name)
CH4Reac.addProductPressure(presion)
CH4Reac.addFirstLaw(True)

#----------------------------------------------

phi = np.arange(0.1, 0.3, 0.01)
phi = np.round(phi, 2)
N2 = False
#------------------------------------------


resul= np.zeros([len(phi),len(CH4Reac.productNames)])
adiaTemps =np.zeros([len(phi),1])
prodEntalphy =np.zeros([len(phi),1])
CH4Reac.addPhi(phi[0])
CI = genInitCon()
j=0

for i in phi:
    try:
        print('para el phi: ', i, ' --------------------------------')
        print('---->', CI)        
        CH4Reac.addPhi(i)
        CH4Reac.getTotalReacRealElements()
        if i == phi[0]:
            sol = findSolution(CH4Reac, CI, args=[500])
        else:
            sol = findSolution(CH4Reac, CI, args=[CI[-1]])
        if CH4Reac.solveWithEnergy:
            resul[j,:]=sol[0:-1]/sum(sol[0:-1])
            adiaTemps[j,:] = sol[-1]    
        else:
            resul[j,:]=sol/sum(sol)
        CI = sol
        print(sol)
        prodEntalphy[j, :] = (CH4Reac.get_cp_prod(sol[-1]))
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
            temperatura = 700
            filename = 'equlibrium_result_CH4_at_P_' + str(presion) + '_T_'+ str(temperatura) + '_' + str(phi[0]) + '_phi_'+ str(i-0.01) + '.csv'
            df.to_csv(filename, index=False)
        break

    
if CH4Reac.solveWithEnergy:
    fig, axs = plt.subplots(2)
else:
    fig, axs = plt.subplots(1)
     
    
for k in range(len(CH4Reac.productNames)):
    
    if CH4Reac.solveWithEnergy:
        if N2 == False and CH4Reac.productNames[k] == 'N2':
            continue
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
plt.savefig('equlibrium_result_CH4_at_P_' + str(presion) + '_T_ad.png', dpi=600, bbox_inches='tight')
plt.show()


#-------------------------------------------------------------------------------------------------------------------

if CH4Reac.solveWithEnergy:
    titles = ['$phi$', '$T_ad$', 'CO2', 'CO', 'H2O', 'H2', 'OH', 'O2', 'N2', 'H', 'O', 'h_prod']
    result = np.hstack((adiaTemps, resul))
else:
    result = resul
    titles = ['$phi$','CO2', 'CO', 'H2O', 'H2', 'OH', 'O2', 'N2', 'H', 'O']

result = np.hstack((phi.reshape(len(phi),1), result))
result = np.hstack((result, prodEntalphy))
result = np.vstack((titles, result))
df = pd.DataFrame(result)

if CH4Reac.solveWithEnergy:
    filename = 'equlibrium_result_CH4_at_P_' + str(presion) + '_T_ad----.csv'
else:
    filename = 'equlibrium_result_CH4_at_P_' + str(presion) + '_T_'+ str(temperatura) +'.csv'
df.to_csv(filename, index=False)
