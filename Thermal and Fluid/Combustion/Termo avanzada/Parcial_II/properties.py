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
    return -74850 
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

def hf_reactivos():
    return np.array([hf_CH4(), hf_C2H6(), hf_C3H8(), hf_C4H10(), hf_C5H12(), 
                     hf_C6H14(), hf_C7H16(), hf_CO2(), hf_N2()])
def hf_productos():
    return np.array([hf_CO2(), hf_CO(), hf_H2O(), hf_H2(), hf_OH(), hf_O2(), 
                     hf_NO(), hf_N2()])

def deltaHp(T):
    deltaH_CO2 = Th.cal_property(T, 'CO2', 'h') - Th.cal_property(298.15, 'CO2', 'h')
    deltaH_CO = Th.cal_property(T, 'CO', 'h') - Th.cal_property(298.15, 'CO', 'h')
    deltaH_H2O = Th.cal_property(T, 'H2O', 'h') - Th.cal_property(298.15, 'H2O', 'h')
    deltaH_H2 = Th.cal_property(T, 'H2', 'h') - Th.cal_property(298.15, 'H2', 'h')
    deltaH_OH = Th.cal_property(T, 'OH', 'h') - Th.cal_property(298.15, 'OH', 'h')
    deltaH_O2 = Th.cal_property(T, 'O2', 'h') - Th.cal_property(298.15, 'O2', 'h')
    deltaH_NO = Th.cal_property(T, 'NO', 'h') - Th.cal_property(298.15, 'NO', 'h')
    deltaH_N2 = Th.cal_property(T, 'N2', 'h') - Th.cal_property(298.15, 'N2', 'h')
    return np.array([deltaH_CO2, deltaH_CO, deltaH_H2O, deltaH_H2, deltaH_OH, 
                     deltaH_O2, deltaH_NO, deltaH_N2])



def calcT(cte, n, epsilon = 100):
    delta = 1000
    Li = 300
    Lu = 5000
    T = (Li+Lu)/2
    
    while abs(delta) >= epsilon:
        
        deltaH_pro = np.dot(n, deltaHp(T))
        delta = cte - deltaH_pro
        
        if delta <= 0:
            Lu = T
        else:
            Li = T
            
        T = (Li+Lu)/2
    
    return T
