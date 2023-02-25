# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
from scipy.optimize import fsolve

import ThProperties as Th
from matplotlib import pyplot as plt
import properties as prop
from scipy.interpolate import lagrange

n_entrada = np.array([0.7478,
             0.1099,
             0.0512,
             0.0223,
             0.0076,
             0.0027,
             0.0022,
             0.0506,
             0.0057])
# C H O N
comp_entrada = np.array([[1, 4, 0, 0], 
                        [2, 6, 0, 0],
                        [3, 8, 0, 0],
                        [4, 10, 0, 0],
                        [5, 12, 0, 0],
                        [6, 14, 0, 0],
                        [7, 16, 0, 0],
                        [1, 0, 2, 0],
                        [0, 0, 0, 2]])

n_salida_real = np.array([])
n_salida_est = np.array([])

comp_salida_est = np.array([[1, 0, 2, 0],
                            [0, 2, 1, 0],
                            [0, 0, 0, 2],
                            [0, 0, -0.21*2, -0.79*2]])

comp_salida_real = np.array([[1, 0, 2, 0],
                             [1, 0, 1, 0],
                             [0, 2, 1, 0],
                             [0, 2, 0, 0],
                             [0, 1, 1, 0],
                             [0, 0, 2, 0],
                             [0, 0, 0, 2],
                             [0, 0, 1, 1]])
phi = 1.2
com = n_entrada@comp_entrada
com_st = np.linalg.solve(comp_salida_est.T, com)
A_st = com_st[3]
A_r = A_st/phi
C_coef = com[0]
H_coef = com[1]
O_coef = com[2]+A_r*0.42
N_coef = com[3]+A_r*0.79*2

B = 0
C = 0
O = 0
N = 0

for i in range(len(comp_entrada)):
        B += n_entrada[i] * comp_entrada[i][0]
        C += n_entrada[i] * comp_entrada[i][1]
        O += n_entrada[i] * comp_entrada[i][2]
        N += n_entrada[i] * comp_entrada[i][3]

A_est = ((B*2 + C/2) - O)/0.42
D = N + A_est*0.79*2; 



n_salida_est = np.append(n_salida_est, np.array([B, C/2, D, A_est]))



A_real = A_est/phi

#print(n_salida_est)

#------------------------------------------------------------



#Cantidad total de O:

O_real = A_real*2*0.21 + O #Oxígeno total

#print(O_real)

# B CO2 + c CO + + D H2O + E H2 + F OH + G O2 + h NO + I N2

#Balance de O : Ot = 2B + c + D + F

#Ecuaciones de equilibrio
def k1(T): #H2O-OH+1/2H2
    Temps = [1000, 1500, 2000, 2100, 2200, 2500]
    k = [-11.280, -6.344, -3.776, -3.434, -3.091, -2.27]
    p = lagrange(Temps, k)
    n = p(T)
    #n=-11.280*(T==1000) - 6.344*(T==1500)-3.776*(T==2000)-2.27*(T==2500)
    #print(n)
    return 10**n

def k2(T): #CO2+H2-CO+H2O
    Temps = [1000, 1500, 2000, 2100, 2200, 2500]
    k = [-0.159, 0.4035, 0.656, 0.688, 0.316, 0.784]
    p = lagrange(Temps, k)
    n = p(T)
    #n=-0.159*(T==1000)+0.4035*(T==1500)+0.656*(T==2000)+0.784*(T==2500)
    return 10**n

def k3(T): #1/2*O2+1/2*N2-NO
    Temps = [1000, 1500, 2000, 2100, 2200, 2500]
    k = [-4.062, -2.501, -1.699, -1.586, -1.484, -1.227]
    p = lagrange(Temps, k)
    n = p(T)
    #n=-4.062*(T==1000)-2.501*(T==1500)-1.699*(T==2000)-1.227*(T==2500)
    return 10**n

def k4(T): #CO2-CO+1/2*O2
    Temps = [1000, 1500, 2000, 2100, 2200, 2500]
    k = [-10.221, -5.36, -2.884, -2.539, -2.226, -1.44]
    p = lagrange(Temps, k)
    n = p(T)
    #n=-10.221*(T==1000)-5.36*(T==1500)-2.884*(T==2000)-1.44*(T==2500)
    return 10**n


def equations(p):

    """
    Se usan las 4 ecuaciones de balance de masa y 4 de equilibrio para 
    encontrar cada coeficiente, los coeficientes corresponden a:
        B1: CO2
        c1: CO
        D1: H2O
        E1: H2
        F1: OH
        G1: O2
        h1: NO
        I1: N2
    """
    f = (P/101.325) # Presión del sistema sobre la de referencia
    B1, c1, D1, e1, F1, G1, h1, I1 = p
    n = B1+c1+D1+e1+F1+G1+h1+I1# Moles totales en productos
    eqn1 = B1+c1-C_coef # Balance de carbonos
    eqn2 = 2*D1+2*e1+F1-H_coef #Balance Hidrógeno
    eqn3 = 2*B1+c1+D1+F1+ 2*G1 +h1-O_coef # Balnce de oxígenos
    eqn4 = h1 + 2*I1-N_coef # Balance de nitrógenos
    
    eqn5 = F1**2*e1*(f/n)-(k1(T))**2*D1**2 # #H2O-OH+1/2H2
    eqn6 = c1*D1-k2(T)*B1*e1 #CO2+H2-CO+H2O
    eqn7 = h1**2 - (k3(T))**2 * G1 * I1 #1/2*O2+1/2*N2-NO
    eqn8 = c1**2 * G1 * (f/n) - (k4(T))**2 * B1**2 #CO2-CO+1/2*O2

    return (eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8)




presion = np.arange(50, 300, 50)  # [kPa]
T_produc = np.arange(1000, 3000, 500)

                       # T  P n1 n2 n3 n4 n5 n6 n7 n8
comp_matrix = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
for i in T_produc: 
    for j in presion:
        T = i
        P = j
        
        B2, c2, D2, E2, F2, G2, h2, I2=  fsolve(equations, (1,1,2,1,0.1,0.1,1,8))
        result = [B2, c2, D2, E2, F2, G2, h2, I2]
        comp_matrix = np.append(comp_matrix, [[T, P, B2, c2, D2, E2, F2, G2, h2, I2]], axis=0)
        
        # print('CO2=',B2)
        # print('CO=',c2)
        # print('H2O=',D2)
        # print('H2=',E2)
        # print('OH=',F2)
        # print('O2=',G2)
        # print('NO=',h2)
        # print('N2=',I2)
        # print('------------------------- \n')


comp_matrix = np.delete(comp_matrix, 0, 0)        
comp_matrix = abs(comp_matrix)

h_com_reac = np.dot(n_entrada, prop.hf_reactivos())

                  #   [T, P, h_com]
h_com_pro = np.array([[0, 0, 0]])
deltaH_pro = np.array([[0, 0, 0]])

for i in range(len(comp_matrix)):
    h_com_i = np.dot(comp_matrix[i, 2:], prop.hf_productos())
    h_com_pro = np.append(h_com_pro, [[comp_matrix[i, 0], comp_matrix[i, 1],
                                       h_com_i]], axis=0)
    
    deltaH_i = np.dot(comp_matrix[i, 2:], prop.deltaHp(comp_matrix[i, 0]))
    deltaH_pro = np.append(deltaH_pro, [[comp_matrix[i, 0], comp_matrix[i, 1],
                                       deltaH_i]], axis=0)
    
h_com_pro = np.delete(h_com_pro, 0, 0)  
deltaH_pro = np.delete(deltaH_pro, 0, 0)  


H_reac = h_com_reac

H_pro = h_com_pro[:, 2] + deltaH_pro[:, 2]

# Calor de combustion: 
Q_com = H_pro - H_reac


LHV = h_com_pro[:, 2] - h_com_reac

eficiencia_com = []

for i in range(len(Q_com)):
    eta_i = Q_com[i]/LHV[i]
    eficiencia_com.append(eta_i)

#Eficiencia de la combustion
eficiencia_com = np.array(eficiencia_com)


T_adiabatica = np.array([])
for i in range(len(comp_matrix)):
    
    RHS = (h_com_reac - h_com_pro[i, 2])
    T = prop.calcT(RHS, comp_matrix[i, 2:])
    T_adiabatica = np.append(T_adiabatica, [T])
 


T_adiabaticai = [np.mean(T_adiabatica), np.mean(T_adiabatica), np.mean(T_adiabatica),
                np.mean(T_adiabatica), np.mean(T_adiabatica)]



T_adiabatica_f = np.array([[0, 0, 0, 0, 0]])
T_adiabatica_f = np.append(T_adiabatica_f, [T_adiabaticai], axis=0)



epsilon = 1E-3
delta = 100
print(len(presion))


while abs(delta) >= epsilon:
    comp_matrix_adia = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

    for j in range(len(presion)):
        T = T_adiabaticai[j]
        P = presion[j]
        B2, c2, D2, E2, F2, G2, h2, I2 =  fsolve(equations, (1,1,2,1,0.1,0.1,1,8))
        result = np.array([[T, P, B2, c2, D2, E2, F2, G2, h2, I2]])
        comp_matrix_adia = np.append(comp_matrix_adia, result, axis=0)
    
    comp_matrix_adia = np.delete(comp_matrix_adia, 0, 0)  
    
    h_com_pro = np.array([[0, 0, 0]])
    
    for i in range(len(comp_matrix_adia)):
        h_com_i = np.dot(comp_matrix_adia[i, 2:], prop.hf_productos())
        h_com_pro = np.append(h_com_pro, [[comp_matrix_adia[i, 0], comp_matrix_adia[i, 1],
                                           h_com_i]], axis=0)
         
    h_com_pro = np.delete(h_com_pro, 0, 0)  
    
    
    T_adiabaticai = np.array([])
    
    for i in range(len(comp_matrix_adia)):
        
        RHS = (h_com_reac - h_com_pro[i, 2])
        Ti = prop.calcT(RHS, comp_matrix_adia[i, 2:])
        T_adiabaticai = np.append(T_adiabaticai, [Ti])
    T_adiabatica_f = np.append(T_adiabatica_f, [T_adiabaticai], axis=0)
    delta = T_adiabatica_f[-1][0] - T_adiabatica_f[-2][0]
    
print(T_adiabatica_f)
        
        
        
plt.figure(1)
plt.plot(presion, T_adiabatica_f[-1], 'ro')
plt.xlabel('Presion [kPa]')
plt.ylabel('Temperatura de llama adiabatica [K]')
plt.title('Temperatura de llama adiabatica para cada presión')
plt.grid()
plt.show()

xlist = presion #Vector de valores x
ylist = T_produc#vector de valores y
X, Y = np.meshgrid(xlist, ylist)#Crea la malla x-y



"""
for i in range(8):
    Z = # valor de composición
    #print(Z)

    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(X, Y, Z)
    fig.colorbar(cp) # Add a colorbar to a plot
       
    #plt.grid(linewidth=1.5,color='k')
    #ax.set_title('Filled Contours Plot')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title('Temperature contours for '+str(self.nx*self.ny)+' cells')
    plt.show()
figure(2)
"""
        
#     T_adiabatica = np.array([])
#     for i in range(len(result)):
        
#         RHS = (h_com_reac - h_com_pro[i, 2])
#         T = prop.calcT(RHS, comp_matrix[i, 2:])
#         T_adiabatica = np.append(T_adiabatica, [T])
     


#     T_adiabatica = np.mean(T_adiabatica)

# print(comp_matrix_adia)
        
     
        
     
###############################################################################





# presion = np.arange(50, 300, 50)  # [kPa]
# T_produc = np.arange(1000, 3000, 500)





        

