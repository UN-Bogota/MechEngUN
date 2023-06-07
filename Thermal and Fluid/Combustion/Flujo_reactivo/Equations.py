#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 18:15:21 2023

@author: julian
"""
import numpy as np
from Properties import *
Ru = 8.31446261815324 # [J/mol.K]
L=1
n=3
dx=L/(n-1)
especies=['CH4','O2', 'CO', 'H2O', 'CO2']

# Matiz de pesos moleculares
MW=[]
for i in range(len(especies)):
    MW.append(getMW(especies[i]))
MW=np.array(MW)
P=101325 #Pa
T0= 1500#K
phi = 1
A_est = 2
A_real = A_est/phi

v_molar_CH4 = Ru*T/P # [m³/mol CH4]
v_molar_CH4 = v_molar_CH4*1000 # [L/mol CH4]


x0= [1, A_real, 0, 0, 0]
x0 = x0*MW/MW[0]

Y0 = x0/np.sum(x0)


MWm = 1/(np.dot(Y0, 1/MW))
rho=P/(Ru*T0)*MWm

v=1 #m/s
l=len(especies)
def create_A(Y):
    
    cells = n*(l+1)
    T = Y[l*n:]
    A = np.zeros([cells,cells])
    C = np.zeros([cells,1])
    A[n*l,n*l]=1
    C[n*l]=T0
    #Av=np.zeros([len(especies),len(especies)])
    for m in range(1,n):
        index = [i*n+m  for i in range(l)]
        MWm = 1/(np.dot(Y[index], 1/MW))
        X=MWm/MW*Y[index]
        Dm=getMultiDiffCoef(X,especies,T[m],P)
        for i in range(l):
            #Término convectivo
            A[n*i+m,n*i+m] = rho*v/dx
            A[n*i+m,n*i+m-1] = -rho*v/dx
            C[n*i+m] = np.dot(getW(rho*Y[index]*1/MW, T[m]).T,(1/MW))
            for j in range(l):
                #Término difusivo
                if not m == (n-1):
                    #No hay difusión en la salida
                    A[n*i+m,j*n+m]+=2*rho*MW[i]/(MWm*dx**2)*Dm[i,j]
                    A[n*i+m,j*n+m+1]+=-rho*MW[i]/(MWm*dx**2)*Dm[i,j]
                    A[n*i+m,j*n+m-1]+=-rho*MW[i]/(MWm*dx**2)*Dm[i,j]
        A[l*n+m,m]=1
        C[l*n+m]=P/(rho*Ru)*MWm
    #BC
    for i in range(l):
        A[n*i,n*i]=1
        C[n*i]=x0[i]
    return A,C
    Ys=np.linalg.solve(A, C)

    return Ys


def Sol():
    Y=np.array([])
    for i in range(l):
        Y0=np.ones(n)*x0[i] # pone las condiciones iniciales para las concentraciones
        Y=np.concatenate((Y, Y0)) #crea un vector de solucion
    Y0=np.ones(n)*T0
    Y=np.concatenate((Y, Y0))
    Y=np.array(Y)
    Ys=Y+1 # condicion para que entre al while
    while abs(np.amax(Y-Ys))>1e-6: #
        
        Y = Ys
        #Ys=create_A(Y)
        A,C = create_A(Y)
        
        return A,C




def getW(concentraciones, T):
    print(concentraciones)
    # Recibe las concentraciones en 
    Ru = 1.987207 # [cal/mol.K]
    
    r1 = (5E11) * np.exp(-47800/(Ru*T)) * (concentraciones[0]**(0.7) * concentraciones[1]**(0.8))
    r2 = (2.24E12) * np.exp(-40700/(Ru*T)) * (concentraciones[2] * concentraciones[3])
    r3 = (5E8) * np.exp(-40700/(Ru*T)) * (concentraciones[4])
    
    omega = np.array([[-r1],
                      [-3/2*r1 - 0.5*r2 + 0.5*r3],
                      [r1 - r2 + r3],
                      [2*r1],
                      [r2-r3]
                      ])
    return omega

A,C=Sol()