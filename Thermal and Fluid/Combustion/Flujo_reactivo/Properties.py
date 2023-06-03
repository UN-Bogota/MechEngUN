#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 18:15:21 2023

@author: HIDES
"""

import numpy as np




def binaryDiffCoef(T, Pressure, especies):
    
    MWa = getMW(especies[0])
    MWb = getMW(especies[1])
    
    sigma_a = getCollisionDiameter(especies[0])
    sigma_b = getCollisionDiameter(especies[1])
    sigma_ab = (sigma_a + sigma_b)/2
    
    MWab = 2 * ((1/MWa) + (1/MWb))**(-1)
    
    A = 1.06036
    B = 0.15610
    C = 0.19300
    D = 0.47635
    E = 1.03587
    F = 1.52996
    G = 1.76474
    H = 3.89411
    
    Kb = 1.380649E-23
    
    epsilon_a = getEpsilon(especies[0])
    epsilon_b = getEpsilon(especies[1])
    
    Tast = Kb*T / (epsilon_a*epsilon_b)**(1/2)
    
    omegaD = A / (Tast**B) + C / (np.exp(D*Tast)) + E / (np.exp(F*Tast)) + G / (np.exp(H*Tast))
    
    Dab = 0.0266 * T**(3/2) / (
        Pressure*MWab**(1/2) * sigma_ab**2 * omegaD)
    
    return Dab
    
    
    

def getMW(specie):
    
    if specie == 'CH4':
        MW = 16.043
    elif specie == 'O2':
        MW = 31.999
    elif specie == 'CO':
        MW = 28.011
    elif specie == 'H2O':
        MW = 18.015
    elif specie == 'CO2':
        MW = 44.01
    elif specie == 'H2':
        MW = 2.016
    elif specie == 'N2':
        MW = 28.013
    
    return MW
        
        
def getCollisionDiameter(specie):
    """
    

    Parameters
    ----------
    specie : str
        Name of specie

    Returns
    -------
    MW : float
        Hard-Sphere Collision diameter

    """
    
    if specie == 'CH4':
        sigma = 3.758
    elif specie == 'O2':
        sigma = 3.467
    elif specie == 'CO':
        sigma = 3.690
    elif specie == 'H2O':
        sigma = 2.641
    elif specie == 'CO2':
        sigma = 3.941
    elif specie == 'H2':
        sigma = 2.827
    elif specie == 'N2':
        sigma = 3.798
    
    return sigma
        
def getEpsilon(specie):
    
    Kb = 1.380649E-23
    
    if specie == 'CH4':
        epsilon = 148.6
    elif specie == 'O2':
        epsilon = 106.7
    elif specie == 'CO':
        epsilon = 91.7
    elif specie == 'H2O':
        epsilon = 809.1
    elif specie == 'CO2':
        epsilon = 195.2
    elif specie == 'H2':
        epsilon = 59.7
    elif specie == 'N2':
        epsilon = 71.4
    
    return epsilon*Kb

def getFmatix(concentraciones, especies, T, Pressure):
    
    len_conc = len(concentraciones)
    L = np.zeros((len_conc, len_conc))
    
    for i in range(len_conc):
        MWi = getMW(especies[i])
        for j in range(len_conc):
            MWj = getMW(especies[j])
            
            L_ij = 0
            
            for k in range(len_conc):
                
                a = especies[i]
                b = especies[k]
                
                Dik = binaryDiffCoef(T, Pressure, [a, b])
                
                L_ij += (concentraciones[k] / (MWi*Dik) ) * (MWj*concentraciones[j]*(1-(i==k))
                                                     - MWi*concentraciones[i]*((i==j)-(j==k)))
            L[i, j] = L_ij
            L_ij = 0
           
    F = np.linalg.inv(L)
    
    return F

def getMultiDiffCoef(concentraciones, especies, T, Pressure):
    
    len_conc = len(concentraciones)
    D = np.zeros((len_conc, len_conc))
    F = getFmatix(concentraciones, especies, T, Pressure)
    
    MWmix = 0
    for especie in range(len_conc):
        
        MWmix += concentraciones[especie]*getMW(especies[especie])
        
    for i in range(len_conc):
        for j in range(len_conc):
            MWj = getMW(especies[j])
            D[i, j] = concentraciones[i]*(MWmix/MWj)*(F[i, j] - F[i, i])
       
    print(D)
    return D

    


################# EJEMPLO #################
concentraciones = [0.15, 0.2, 0.65]
especies = ['H2', 'O2', 'N2']

T = 600

P = 101325

D = getMultiDiffCoef(concentraciones, especies, T, P)    
    
######################################
    
    
    