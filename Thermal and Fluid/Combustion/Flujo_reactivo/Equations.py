#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 18:15:21 2023

@author: julian
"""
import numpy as np
import matplotlib.pyplot as plt
from Properties import *
from ThProperties import *
from scipy.optimize import fsolve
Ru = 8.31446261815324 # [J/mol.K]
L=5
n=25
dx=L/(n-1)
especies=['CH4','O2', 'CO', 'H2O', 'CO2']

# Matiz de pesos moleculares
MW=[]
for i in range(len(especies)):
    MW.append(getMW(especies[i]))
MW=np.array(MW)
P=101.325 #Pa
T0= 900 #K
phi = 1
A_est = 2
A_real = A_est/phi

v_molar_CH4 = Ru*T/P # [m³/mol CH4]
v_molar_CH4 = v_molar_CH4*1000 # [L/mol CH4]


x0= [1, A_real, 1e-5, 1e-5, 1e-5]
x0 = x0*MW/MW[0]

y0 = x0/np.sum(x0)
Y0=y0.copy()
#print(Y0)

MWm = 1/(np.dot(Y0, 1/MW))
rho=P/(Ru*T0)*MWm

v=2 #m/s
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
        MWm = 1/(np.dot(Y[index].T, 1/MW))
        X=MWm/MW*Y[index]
        #print(MWm)
        #print('for m: ', m, ' X: ', X, ' especies: ', especies, ' Temp: ', T[m], ' press: ', P)
        Dm=getMultiDiffCoef(X,especies,T[m],P)
        omega=getW(rho*Y[index]*1/MW, T[m]).T*MW
        #print(rho*Y[index]*1/MW)
        kc=[]
        km=[]
        kp=[]
        hc=[]
        hm=[]
        dYdx=[]
        for i in range(l):
            kc.append(cal_property(T[m], especies[i], 'k'))
            km.append(cal_property(T[m-1], especies[i], 'k'))
            if m<(n-1):
                kp.append(cal_property(T[m+1], especies[i], 'k'))
            else:
                kp=kc
            #print(especies[i])
            hc.append(cal_property(T[m], especies[i], 'h'))
            hm.append(cal_property(T[m-1], especies[i], 'h'))
            dYdx.append((Y[index[i]]-Y[index[i]-1])/dx)
            
            #Término convectivo
            A[n*i+m,n*i+m] = rho*v/dx
            A[n*i+m,n*i+m-1] = -rho*v/dx
            C[n*i+m] = omega[0][i]
            
            for j in range(l):
                
                #Término difusivo
                if not m == (n-1):
                    #No hay difusión en la salida
                    A[n*i+m,j*n+m]+=2*rho*MW[i]/(MWm*dx**2)*Dm[i,j]
                    A[n*i+m,j*n+m+1]+=-rho*MW[i]/(MWm*dx**2)*Dm[i,j]
                    A[n*i+m,j*n+m-1]+=-rho*MW[i]/(MWm*dx**2)*Dm[i,j]
                    
        vi=(Dm@dYdx)/X
        mpp=Y[index]*(v+vi)*rho
        hc=np.array(hc)/MW*1000
        hm=np.array(hm)/MW*1000
        dh=(hc-hm)/dx
        kc=np.array(kc)
        km=np.array(km)
        kp=np.array(kp)
        hc=np.array(hc)
        hm=np.array(hm)
        kmp=(kc+kp)/2
        kmp=np.dot(X,kmp)
        kmm=(kc+km)/2
        kmm=np.dot(X,kmm)
        A[l*n+m,l*n+m]=-(kmp+kmm)/(dx**2)
        A[l*n+m,l*n+m-1]=kmm/(dx**2)
        if not m == (n-1):
            A[l*n+m,l*n+m+1]=kmp/(dx**2)
        C[l*n+m]=-np.dot(mpp,dh)+np.dot(hc,omega[0])
        
    #BC
    for i in range(l):
        A[n*i,n*i]=1
        C[n*i]=y0[i]
    #return A,C
    #Ys=np.linalg.solve(A, C)
    #print(Ys)
    return A,C#Ys.T[0]


def Sol(x):


#    Ys=Y+1
    
    # while abs(np.amax(Y-Ys))>1e-6:
    #     Y=Ys
    #Ys=create_A(Y)
    A,C=create_A(x)
    x=x.reshape(len(x),1)
    a=(np.matmul(A,x)-C).flatten()
    #print(a)
    return a 
    #     print(Ys)
    #return Ys
#Acoplar con cinética (W)


def getW(concentraciones, T):
    #print(concentraciones)
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
    #print(omega)
    return omega
#Sol()
Y=np.array([])
for i in range(l):
    Y0=np.ones(n)*y0[i] # pone las condiciones iniciales para las concentraciones
    Y=np.concatenate((Y, Y0)) #crea un vector de solucion
Y0=np.linspace(1,1.2,n)*T0
Y=np.concatenate((Y, Y0))
#print(Sol(Y))
Ys = fsolve(Sol, Y)
x=np.linspace(0, L,n)
print(Ys)
#Y=np.array(Y)
plt.figure(1)
plt.plot(x,Ys[:n])
plt.title('$CH_4$')
plt.figure(2)
plt.plot(x,Ys[n:2*n])
plt.title('$O_2$')
plt.figure(3)
plt.plot(x,Ys[2*n:3*n])
plt.title('$CO$')
plt.figure(4)
plt.plot(x,Ys[3*n:4*n])
plt.title('$H_2O$')
plt.figure(5)
plt.plot(x,Ys[4*n:5*n])
plt.title('$CO_2$')
plt.figure(6)
plt.plot(x,Ys[5*n:6*n])
plt.title('Temperatura')
