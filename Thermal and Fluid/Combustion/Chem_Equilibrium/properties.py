#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 00:12:11 2022

@author: bojack
"""

#####   ENTALPIAS DE FORMACION     ###

# REACTIVOS ------------
import numpy as np
import ThProperties as Th

def hf_CH4():    
    return 
def hf_C2H6():
    return -84680
def hf_C3H8():
    return -103850
def hf_C4H10(): 
    return -126150
def hf_C5H12(): 
    return -146900
def hf_C6H14(): 
    return -167200
def hf_C7H16(): 
    return -187000

# COMUN --------------
def hf_CO2(): 
    return -393520
def hf_N2(): 
    return 0
def hf_O2():
    return 0

# PRODUCTOS ------------------


def hf_CO(): 
    return -110530
def hf_H2O(): 
    return -241820
def hf_H2(): 
    return 0
def hf_OH(): 
    return 39460
def hf_NO(): 
    return 88850
species = {
    'CH4': -74850,
    'C2H6': -103850,
    'C3H8': -103850,
    'C4H10': -126150,
    'C5H12': -146900,
    'C6H14': -167200,
    'C7H16': -187000,
    'CO2': -393520,
    'N2': 0,
    'O2': 0,
    'CO': -110530,
    'H2O': -241820,
    'H2': 0,
    'OH': 39460,
    'NO': 88850
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
        deltaH_values.append(Th.cal_property(T, i, 'h') - Th.cal_property(298.15, i, 'h'))

    return np.array(deltaH_values)



def calcT(cte, n, epsilon = 100):
    delta = 1000
    Li = 300
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
