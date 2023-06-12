# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 19:33:15 2023

@author: jolej
"""

from Properties import *
from scipy.optimize import fsolve
Ru = 8.31446261815324 # [J/mol.K]
L=1
n=20
dx=L/(n-1)
especies=['CH4','O2', 'CO', 'H2O', 'CO2']

# Matiz de pesos moleculares
MW=[]
for i in range(len(especies)):
    MW.append(getMW(especies[i]))
MW=np.array(MW)
P=101.325 #Pa
T0= 1000 #K
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

v=1 #m/s
l=len(especies)
def create_A(Y):
    
    cells = n*(l+1)
    T = Y[l*n:]
    A = np.zeros([cells,cells])
    C = np.zeros([cells,1])
    A[n*l,n*l]=1
    C[n*l]=T0
    #print(Y)
    #Av=np.zeros([len(especies),len(especies)])
    for m in range(1,n):
        index = [i*n+m  for i in range(l)]
        MWm = 1/(np.dot(Y[index].T, 1/MW))
        X=MWm/MW*Y[index]
        #print(MWm)
        #print('for m: ', m, ' X: ', X, ' especies: ', especies, ' Temp: ', T[m], ' press: ', P)
        Dm=getMultiDiffCoef(X,especies,T[m],P)
        #print(Y[index],T)
        omega=getW(rho*Y[index]*1/MW, T[m]).T*MW
        #print(omega[0])
        for i in range(l):
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
        A[l*n+m,l*n+m]=1
        C[l*n+m]=P/(rho*Ru)*MWm
        #print(MWm)
    #BC
    for i in range(l):
        A[n*i,n*i]=1
        C[n*i]=y0[i]
    #return A,C
    #Ys=np.linalg.solve(A, C)
    #print(Ys)
    return A,C#Ys.T[0]


def Sol(x):
    #x=x.reshape(len(x),1)
    Y=np.array([])
    for i in range(l):
        Y0=np.ones(n)*y0[i] # pone las condiciones iniciales para las concentraciones
        Y=np.concatenate((Y, Y0)) #crea un vector de solucion
    Y0=np.ones(n)*T0
    Y=np.concatenate((Y, Y0))
    #print(Sol(Y))
    Ys=Y+1
    
    while abs(np.amax(Y-Ys))>1e-6:
        Y=Ys
        print(np.amax(Y-Ys))
        A,C=create_A(Y)
        Ys=np.linalg.solve(A, C)
        Ys=np.where(Ys <= 0, 2e-6, Ys)
        Ys=Ys.flatten()
        #print(Ys)
        
        #  print(Ys)
    #A,C=create_A(x)
    
    #a=(np.matmul(A,x)-C).flatten()
    #print(a)
    return Ys,A,C
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
Ys,A,C=Sol(1)
# Y=np.array([])
# for i in range(l):
#     Y0=np.ones(n)*y0[i] # pone las condiciones iniciales para las concentraciones
#     Y=np.concatenate((Y, Y0)) #crea un vector de solucion
# Y0=np.ones(n)*T0
# Y=np.concatenate((Y, Y0))
#print(Sol(Y))
#Ys = fsolve(Sol, Y)
print(Ys)