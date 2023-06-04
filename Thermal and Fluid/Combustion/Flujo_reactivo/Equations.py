#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 18:15:21 2023

@author: julian
"""
import numpy as np
L=1
n=10
dx=L/(n-1)
especies=['','','','','','','']
MW=['CH4',' O2', 'CO', 'H2O', 'CO2']
def create_A(Y):
    l=len(especies)
    cells=n*(l+1)
    A=np.zeros([cells])
    C=np.zeros([cells,1])
    #Av=np.zeros([len(especies),len(especies)])
    for m in range(1,n):
        index=[i*n+m  for i in range(l)]
        MWm=1/(np.dot(Y[index], 1/MW))
        for i in range(l):
            #Término convectivo
            A[n*i+m,n*i+m]=rho*v/dx
            A[n*i+m,n*i+m-1]=-rho*v/dx
            C[n*i+m]=w(Y,T)
            for j in range(l):
                #Término difusivo
                if not m==(n-1):
                    #No hay difusión en la salida
                    A[n*i+m,j*n+m]+=2*rho*MW[i]/(MWm*dx**2)*Dm[i,j]
                    A[n*i+m,j*n+m+1]+=-rho*MW[i]/(MWm*dx**2)*Dm[i,j]
                    A[n*i+m,j*n+m-1]+=-rho*MW[i]/(MWm*dx**2)*Dm[i,j]
        A[l*n+m,m]=1
        C[l*n+m]=P/(rho*R)*MWm
    #BC
    for i in range(l):
        A[n*i,n*i]=1
        C[n*i]=x0[i]
        
    Ys=np.linalg.solve(A, C)
    return Ys

def Sol():
    while abs(np.amax(Y-Ys))>1e-6:
        Y=Ys
        Ys=create_A(Y)
        
#Acoplar con cinética (W)

        

    