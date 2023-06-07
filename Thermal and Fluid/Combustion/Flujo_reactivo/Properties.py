#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 18:15:21 2023

@author: HIDES
"""
import numpy as np
import matplotlib.pyplot as plt


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
       
    return D

    


################# EJEMPLO #################
concentraciones = [0.15, 0.2, 0.65]
especies = ['H2', 'O2', 'N2']

T = 600

P = 101325

D = getMultiDiffCoef(concentraciones, especies, T, P)    
    
######################################
# Segunda implementación
# especies=['H2','O2','N2']
# MW=np.array([2,32,28])
# L00_00=np.zeros([len(especies),len(especies)])
# X=np.array([0.15,0.2,0.65])
# T=600 #K
# P= 1 #atm
# MWm=np.dot(X,MW)
# D=np.array([[0,2.5668,2.4095],[2.5668,0,0.6753],[2.4095,0.6753,0]])#/(100**2)
# for i in range(len(especies)):
#     for j in range(len(especies)):
#         x0=0
#         x1=0
#         x2=0
#         x3=0
#         x4=0
#         if not i==j: 
#             #print(i,j)
#             for k in range(len(especies)):
#                 if not i==k:
#                     x0+=X[k]/(MW[i]*D[i,k])*(MW[j]*X[j]*(1-1*(i==k))-MW[i]*X[i]*(1*(i==j)-1*(j==k)))
#                 #x1+=X[j]*X[k]*(1*(i==j)-1*(j==k))*MW[k]*(1.2*Cjk-1)/((MW[j]+MW[k])*D[j,k])
#                 #x2+=MW[i]/MW[j]*X[i]*X[k]/((MW[i]+MW[k])**2*D[i,k])*((1*(i==k)-1*(j==i))*(15/2*MW[j]**2+25/4*MW[k]**2-3*MW[k]*Bik)-4*MW[j]*MW[k]*Aik*((1*(i==k)-1*(j==i)))*(1+5/(3*np.pi())*(c_rot[i]/(kB*C)+c_rot[k]/(kB*50))))
#                 #x3+=MW[j]*Ajk/((MW[j]+MW[k])*D[j,k])(1*(i==k)+1*(j==i))*X[j]*X[k]*c_rot[j]/50
#                 #x4+=X[i]*X[k]/D[i,k]+12/(5*np.pi()*c_int[i]*MW[k]*D[i,k]*50)*X[i]*X[k]*MW[i]*Aik*c_rot[i]
#         L00_00[i,j]=16/25*T/P*x0
#         #L00_10[i,j]=8/5*T/P*x1
#         #L10_10[i,j]=16/25*T/P*x2
#         #L10_01[i,j]=32*T/(5*np.pi()*c_int[j])*x3
#         #L01_01[i,j]=-4*kB*T/(c_int[i]*P)*x4
# Dm=np.zeros([len(especies),len(especies)])
# F = np.linalg.inv(L00_00)
# for i in range(len(especies)):
#     for j in range(len(especies)):
#         Dm[i,j]=X[i]*16/25*T/P*MWm/MW[j]*(F[i,j]-F[i,i])
    
