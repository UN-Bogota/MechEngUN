#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 15:53:23 2024

@author: julian
"""
import scipy 
import numpy as np
import pandas as pd
import cantera as ct
import matplotlib.pyplot as plt
from ic_f import ic_fun
from calculos_presion_21_3_2024 import *

def normal(x):
    x=x/sum(x)
    return x
def entropy(pk,qk):
    D = sum(pk * np.log(pk / qk))
    return D
def data(load,prueba,fuel='B8'):
    d1=dataframes_cop[(dataframes_cop['load']==load)*(dataframes_cop['fuel']==fuel)*(dataframes_cop['prueba']==prueba)]
    P_p=d1.groupby('angle')[['P_abs_cam','rpm']].mean()
    a_p=P_p.index
    P_pe=P_p.to_numpy()
    rpm=P_pe.mean(axis=0)[1]
    P_pa=P_p.to_numpy()[:,0]
    return a_p,P_pa,rpm
#a_p es el ángulo de los datos tomados escogido y corregido para estar entre -360 y 360
#P_pa es el vector de la presión de los datos tomados
def entropy_map(rpm,m_diesel,a_p,P_pa,rv=15.6,p_out=0.765e5, in_val_c=3.2e-6,ang_in_op=-18.,ang_in_clo=225.,out_cal_c=3.6e-6, ang_out_op=522.,ang_out_clo=18.,inj_op=250.,inj_clo=265.):
    """
    (rpm,m_diesel,a_p,P_pa,rv=14.58,p_out=0.72e5, in_val_c=1.e-6,ang_in_op=-18.,ang_in_clo=198.,out_cal_c=1.e-6, ang_out_op=522.,ang_out_clo=18.,inj_op=250.,inj_clo=265.)
    Parameters
    ----------
    rpm : rpm del motor en las condiciones dadas [rpm]
    m_diesel : masa de Diesel inyectada por ciclo ]kg]
    a_p : Vector con todos los ángulos de los datos tomados, normalizados entre -360° y 360°
    P_pa : Vector de presión para los ángulos P_pa [no importa unidades]
    rv : Realción de compresión. The default is 14.58.
    p_out : Presión de salida [Pa].The default is 0.72e5.
    in_val_c : Coeficiente de flujo a la entrada. The default is 1.e-6.
    ang_in_op : ángulo de apertura válvulas de admisión [°] The default is -18.
    ang_in_clo : ángulo de cierre válvulas de admisión [°]  The default is 198..
    out_cal_c : Coeficiente de flujo a la entrada.  The default is 1.e-6.
    ang_out_op : ángulo de apertura válvulas de admisión [°] The default is 522..
    ang_out_clo : ángulo de cierre válvulas de admisión [°]  The default is 18..
    inj_op :  ángulo de apertura inyector de Diesel [°]  The default is 250..
    inj_clo : ángulo de cierre inyector de Diesel[°] The default is 265..

    Returns
    -------
    KL : Entropía relativa promedio

    """
    KLi=[]
    for ki in range(len(m_diesel)):
        # KLi=entropy_map(rev_train[k],m_d_train[k],a_p_train[k],P_pa_train[k],j)  
        # KLj.append(KLi)
        #print(ki)
        try:
            t,a,P=ic_fun(rpm[ki],m_diesel[ki],rv,p_out,in_val_c,ang_in_op,ang_in_clo,out_cal_c, ang_out_op,ang_out_clo,inj_op,inj_clo)
        except: ct.CanteraError
        else:
            inda=np.arange(0,len(a),1)
            i2pick=np.random.choice(inda,720,replace=False)
            ind=[]
            for i in i2pick:
                new=abs(a_p[ki]-a[i])<1.
                ind.append(np.where(new == 1)[0][0])
            Ps=P_pa[ki][ind]/sum(P_pa[ki])
            Qs=P[i2pick]/sum(P)
            #print(type(list(P)))
            #print(entropy(Ps,Qs))
            #print(Ps.shape,Qs.shape)
            #print(scipy.stats.entropy(Ps,Qs))
            KLi.append(scipy.stats.entropy(Ps,Qs))
            #print(KLi)
            
    KL=np.mean(KLi)
    return KL 
#def entropy_plot():
    
    
load_train=np.arange(1,7)
a_p_train=[]
P_pa_train=[]
rev_train=[]
for i in range(len(load_train)):
    a,p,rev=data(load_train[i],2)
    a_p_train.append(a)
    P_pa_train.append(p)
    rev_train.append(rev)

m_d_train=[9.69033636e-06, 8.91993552e-06, 8.63599677e-06, 9.44264069e-06,1.03165939e-05, 1.08600980e-05]


KL=[]
#entropy_map(rev_train,m_d_train,a_p_train,P_pa_train)
#input(' done')
x_rv=np.linspace(14,16.5,15)
for j in x_rv:
    KLj=[]
    n=[]
    KLj=entropy_map(rev_train,m_d_train,a_p_train,P_pa_train,j)
    KL.append(np.mean(KLj))
    print(j)

plt.plot(x_rv,KL)
plt.title('Relación de compresión')
plt.xlabel('rv')
plt.ylabel('Entropía relativa')
plt.vlines(x_rv[np.argmin(KL)],np.min(KL)-0.01,np.max(KL),'r','--')
# # Presión de salida 
# KL_pout=[]
# x_pout=np.linspace(7.5e5,8e5,10)
# for j in x_pout:
#     KLj=[]
#     n=[]
#     KLj=entropy_map(rev_train,m_d_train,a_p_train,P_pa_train,p_out=j)
#     KL_pout.append(np.mean(KLj))
#     print(j)

# plt.scatter(x_pout/1e6,KL_pout)
# plt.title('Presión de salida')
# plt.xlabel('P_out [kPa]')
# plt.ylabel('Entropía relativa')
#plt.vlines(x_pout[np.argmin(KL_pout)],np.min(KL_pout)-0.01,np.max(KL_pout),'r','--')
# #Coeficiente de flujo a la entrada
# KL_cinval=[]
# x_cinval=np.arange(2e-6,4e-6,2e-7)
# for j in x_cinval:
#     KLj=[]
#     n=[]
#     KLj=entropy_map(rev_train,m_d_train,a_p_train,P_pa_train,in_val_c=j)
#     KL_cinval.append(np.mean(KLj))
#     print(j)

# plt.plot(x_cinval,KL_cinval)
# plt.title('Coeficiente de las válvulas a la entrada')
# plt.xlabel('C_val_in')
# plt.ylabel('Entropía relativa')
# plt.vlines(x_cinval[np.argmin(KL_cinval)],np.min(KL_cinval)-0.01,np.max(KL_cinval),'r','--')
#Ańgulo de apertura de válvulas de admisión
# KL_anvalin=[]
# x_anvalin=np.linspace(-25,-10,16)
# for j in x_anvalin:
#     KLj=[]
#     n=[]
#     KLj=entropy_map(rev_train,m_d_train,a_p_train,P_pa_train,ang_in_op=j)
#     KL_anvalin.append(np.mean(KLj))
#     print(j)

# plt.plot(x_anvalin,KL_anvalin)
# plt.title('Ángulo de apertura de válvulas de admisión')
# plt.xlabel('Ángulo')
# plt.ylabel('Entropía relativa')
# plt.vlines(x_anvalin[np.argmin(KL_anvalin)],np.min(KL_anvalin)-0.01,np.max(KL_anvalin),'r','--')

#Ańgulo de cierre de válvulas de admisión
# KL_anvalinc=[]
# x_anvalinc=np.linspace(190,240,25)
# for j in x_anvalinc:
#     KLj=[]
#     n=[]
#     KLj=entropy_map(rev_train,m_d_train,a_p_train,P_pa_train,ang_in_clo=j)
#     KL_anvalinc.append(np.mean(KLj))
#     print(j)

# plt.plot(x_anvalinc,KL_anvalinc)
# plt.title('Ángulo de cierre de válvulas de admisión')
# plt.xlabel('Ángulo')
# plt.ylabel('Entropía relativa')
# plt.vlines(x_anvalinc[np.argmin(KL_anvalinc)],np.min(KL_anvalinc)-0.01,np.max(KL_anvalinc),'r','--')

# KL_coutval=[]
# x_coutval=np.arange(2e-6,4e-6,2e-7)
# for j in x_coutval:
#     KLj=[]
#     n=[]
#     KLj=entropy_map(rev_train,m_d_train,a_p_train,P_pa_train,out_cal_c=j)
#     KL_coutval.append(np.mean(KLj))
#     print(j)

# plt.plot(x_coutval,KL_coutval)
# plt.title('Coeficiente de las válvulas a la salida')
# plt.xlabel('C_val_out')
# plt.ylabel('Entropía relativa')
# plt.vlines(x_coutval[np.argmin(KL_coutval)],np.min(KL_coutval)-0.01,np.max(KL_coutval),'r','--')

#Ańgulo de apertura de válvulas de exosto
KL_anvalout=[]
x_anvalout=np.linspace(480,530,25)
for j in x_anvalout:
    KLj=[]
    n=[]
    KLj=entropy_map(rev_train,m_d_train,a_p_train,P_pa_train,ang_out_op=j)
    KL_anvalout.append(np.mean(KLj))
    print(j)

plt.plot(x_anvalout,KL_anvalout)
plt.title('Ángulo de apertura de válvulas de exosto')
plt.xlabel('Ángulo')
plt.ylabel('Entropía relativa')
plt.vlines(x_anvalout[np.argmin(KL_anvalout)],np.min(KL_anvalout)-0.01,np.max(KL_anvalout),'r','--')

#Ańgulo de cierre de válvulas de exosto
KL_anvaloutc=[]
x_anvaloutc=np.linspace(-10,20,25)
for j in x_anvaloutc:
    KLj=[]
    n=[]
    KLj=entropy_map(rev_train,m_d_train,a_p_train,P_pa_train,ang_in_op=j)
    KL_anvaloutc.append(np.mean(KLj))
    print(j)

plt.plot(x_anvaloutc,KL_anvaloutc)
plt.title('Ángulo de cierre de válvulas de exosto')
plt.xlabel('Ángulo')
plt.ylabel('Entropía relativa')
plt.vlines(x_anvaloutc[np.argmin(KL_anvaloutc)],np.min(KL_anvaloutc)-0.01,np.max(KL_anvalioutc),'r','--')
