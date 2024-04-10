#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 18:05:43 2023

@author: bojack
"""

import random

from reactionObject import Reaction

def genInitCon(solveWithEnergy = False):
    
    #random.seed(10, CO H2O, H2, OH, O2, N2, NO, C3H8, H, O, N)
    n = random.random()*10
    # CI = [CO2]
    #CI = [0.4*n, 0.01*n, 0.45*n, 0.1*n, 0.02*n, 0.1*n, 0.7*n, 0.1*n,
    #      0.001*n, 0.01*n, 0.1*n, 0.01*n]
    
    if solveWithEnergy:
        CI = [0.4*n, 0.01*n, 0.45*n, 0.1*n, 0.02*n, 0.1*n, 0.7*n, 0.1*n,
              0.001*n, 0.01*n, 0.1*n, 0.01*n, 2200]
    
    else:
        CI = [0.4*n, 0.01*n, 0.45*n, 0.1*n, 0.02*n, 0.1*n, 0.7*n, 0.1*n,
              0.001*n, 0.01*n, 0.1*n, 0.01*n]
    return CI

def molarProdFrac(products):
    
    molarProdFrac = []
    N_tot = sum(products)
    
    for i in products:
        molarProdFrac.append(i/N_tot)
    return molarProdFrac

def findSolution(reaction, epsilon = 1E-6):
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
        sol_i = reaction.solveSystem(CI)

        deltaMass = reaction.getReactmass() - reaction.getProdMass()
        
        count = True
        
        for i in sol_i[0]:
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
metanoReac = Reaction()

# Se definen las especies del combustible
methane = [1, 4, 0, 0]
methane_name = 'CH4'

hydrogen = [0, 2, 0, 0]

phi = 0.9
metanoReac.addFuelSpecies(methane, methane_name)

metanoReac.addPhi(phi)
metanoReac.addFirstLaw(True)

sol = findSolution(metanoReac)
print(sol)

