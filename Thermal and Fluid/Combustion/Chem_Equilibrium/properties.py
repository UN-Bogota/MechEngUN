#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 00:12:11 2022

@author: bojack
"""

#####   ENTALPIAS DE FORMACION [kJ/kmol] ###

# REACTIVOS ------------
import numpy as np
import ThProperties as Th

def hf_CH4():    
    return 
def hf_C2H6():
    return -84680.0
def hf_C3H8():
    return -103850.0
def hf_C4H10(): 
    return -126150.0
def hf_C5H12(): 
    return -146900.0
def hf_C6H14(): 
    return -167200.0
def hf_C7H16(): 
    return -187000.0

# COMUN --------------
def hf_CO2(): 
    return -393520.0
def hf_N2(): 
    return 0.0
def hf_O2():
    return 0.0

# PRODUCTOS ------------------


def hf_CO(): 
    return -110530.0
def hf_H2O(): 
    return -241820.0
def hf_H2(): 
    return 0.0
def hf_OH(): 
    return 39460.0
def hf_NO(): 
    return 88850.0

species = {
    'CH4': -74850.0,
    'C2H6': -103850.0,
    'C3H8': -103850.0,
    'C4H10': -126150.0,
    'C5H12': -146900.0,
    'C6H14': -167200.0,
    'C7H16': -187000.0,
    'CO2': -393520.0,
    'N2': 0.0,
    'O2': 0.0,
    'CO': -110530.0,
    'H2O': -241820.0,
    'H2': 0.0,
    'OH': 39460.0,
    'NO': 88850.0,
    'N' : 472629.0,
    'O': 249197.0,
    'H': 217977.0,
    'HO2': 2090.0,     #PENDIENTE en ThProoerties
    'H2O2': -136110.0   #PENDIENTE en ThProoerties
    }

def hf_reactivos(reacSpecies):
    hfReac = []
    for i in reacSpecies:
        hfReac.append(species[i])
    return np.array(hfReac)

def hf_productos(prodSpecies):
    hfProd = []
    for i in prodSpecies:
        hfProd.append(species[i])
    return np.array(hfProd)

def deltaH(T, elements):
    
    deltaH_values = []
    
    for i in elements:
        try:
            value = Th.cal_property(T, i, 'h') - Th.cal_property(298.15, i, 'h')
            deltaH_values.append(value)
        except UnboundLocalError:
            print('Se putio con ', i, 'a T= ', T)
            return np.zeros(len(elements))
    return np.array(deltaH_values)

def calcT(cte, n, epsilon = 100):
    delta = 1000
    Li = 290
    Lu = 5000
    T = (Li+Lu)/2
    
    while abs(delta) >= epsilon:
        
        deltaH_pro = np.dot(n, deltaH(T))
        delta = cte - deltaH_pro
        
        if delta <= 0:
            Lu = T
        else:
            Li = T
            
        T = (Li+Lu)/2
    
    return T