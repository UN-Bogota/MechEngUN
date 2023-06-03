#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 18:15:21 2023

@author: julian
"""
import numpy as np
import matplotlib.pyplot as plt
# Multicomponente
especies=['H2','O2','N2']
MW=np.array([2,32,28])
L00_00=np.zeros([len(especies),len(especies)])
X=np.array([0.15,0.2,0.65])
T=600 #K
P= 1 #atm
MWm=np.dot(X,MW)
D=np.array([[0,2.5668,2.4095],[2.5668,0,0.6753],[2.4095,0.6753,0]])#/(100**2)
for i in range(len(especies)):
    for j in range(len(especies)):
        x0=0
        x1=0
        x2=0
        x3=0
        x4=0
        if not i==j: 
            #print(i,j)
            for k in range(len(especies)):
                if not i==k:
                    x0+=X[k]/(MW[i]*D[i,k])*(MW[j]*X[j]*(1-1*(i==k))-MW[i]*X[i]*(1*(i==j)-1*(j==k)))
                #x1+=X[j]*X[k]*(1*(i==j)-1*(j==k))*MW[k]*(1.2*Cjk-1)/((MW[j]+MW[k])*D[j,k])
                #x2+=MW[i]/MW[j]*X[i]*X[k]/((MW[i]+MW[k])**2*D[i,k])*((1*(i==k)-1*(j==i))*(15/2*MW[j]**2+25/4*MW[k]**2-3*MW[k]*Bik)-4*MW[j]*MW[k]*Aik*((1*(i==k)-1*(j==i)))*(1+5/(3*np.pi())*(c_rot[i]/(kB*C)+c_rot[k]/(kB*50))))
                #x3+=MW[j]*Ajk/((MW[j]+MW[k])*D[j,k])(1*(i==k)+1*(j==i))*X[j]*X[k]*c_rot[j]/50
                #x4+=X[i]*X[k]/D[i,k]+12/(5*np.pi()*c_int[i]*MW[k]*D[i,k]*50)*X[i]*X[k]*MW[i]*Aik*c_rot[i]
        L00_00[i,j]=16/25*T/P*x0
        #L00_10[i,j]=8/5*T/P*x1
        #L10_10[i,j]=16/25*T/P*x2
        #L10_01[i,j]=32*T/(5*np.pi()*c_int[j])*x3
        #L01_01[i,j]=-4*kB*T/(c_int[i]*P)*x4
Dm=np.zeros([len(especies),len(especies)])
F = np.linalg.inv(L00_00)
for i in range(len(especies)):
    for j in range(len(especies)):
        Dm[i,j]=X[i]*16/25*T/P*MWm/MW[j]*(F[i,j]-F[i,i])