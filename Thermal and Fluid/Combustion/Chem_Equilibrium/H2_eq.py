#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 11:23:50 2023

@author: HIDES
"""

import random
import numpy as np

from reactionObject import Reaction

def genInitCon(solveWithEnergy = False):
    
    #random.seed(10, CO H2O, H2, OH, O2, N2, NO, C3H8, H, O, N)
    n = random.random()*10
    
    if solveWithEnergy:
        CI = [0.4*n, 0.01*n, 0.45*n, 0.1*n, 0.02*n, 0.1*n, 0.1*n, 0.1*n,
              0.8*n, 2000]
    
    else:
        CI = [0.4*n, 0.01*n, 0.45*n, 0.1*n, 0.02*n, 0.1*n, 0.1*n, 0.1*n,
              0.8*n]
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
phi = 0.5
temperatura = 1000
presion = 101.325


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

H2Reac.addPhi(phi)
H2Reac.addProductTemperature(temperatura)


H2Reac.addFirstLaw(False)

sol = findSolution(H2Reac)


print(sol)

print(sol/sum(sol))

print(sum(sol/sum(sol)))