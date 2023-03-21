#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 18:05:43 2023

@author: bojack
"""

import random

from reactionObject import Reaction

def genInitCon():
    
    """
    Se usan las 4 ecuaciones de balance de masa y 8 de equilibrio para 
    encontrar cada coeficiente, los coeficientes corresponden a:
                
        B: CO2
        C: CO
        D: H2O
        E: H2
        F: OH
        G: O2
        H: N2
        I: NO
        J: C3H8
        K: H
        L: O
        M: N
        falta el NO2
    """
    #random.seed(10)
    n = random.random()*10
    
    CI = [0.4*n, 0.01*n, 0.45*n, 0.1*n, 0.02*n, 0.1*n, 0.7*n, 0.1*n,
          0.001*n, 0.01*n, 0.1*n, 0.01*n]
    return CI

def findSolution(reaction, epsilon = 1E-6):
    
    deltaMass = 1
    noNegativity = False
    
    sol = None
    
    while (not noNegativity) and  abs(deltaMass) >= epsilon :
       
        CI = genInitCon()
        sol_i = reaction.solveSystem(CI)
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

def molarProdFrac(products):
    
    molarProdFrac = []
    N_tot = sum(products)
    for i in products:
        molarProdFrac.append(i/N_tot)
    return molarProdFrac

# Creo el objeto
metanoReac = Reaction()

# Se definen las especies del combustible
methane = [1, 4, 0, 0]
phi = 0.9
metanoReac.addFuelSpecies(methane)

metanoReac.addPhi(phi)

sol = findSolution(metanoReac)
print(sol)

