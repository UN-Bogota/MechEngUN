#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 07:53:34 2023

@author: bojack
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
P_1_Tad = pd.read_csv('equlibrium_result_at_P_101.325_T_ad.csv')
P_1_Tad = np.array(P_1_Tad)

P_5_Tad = pd.read_csv('equlibrium_result_at_P_506.625_T_ad.csv')
P_5_Tad = np.array(P_1_Tad)

P_10_Tad = pd.read_csv('equlibrium_result_at_P_1013.25_T_ad.csv')
P_10_Tad = np.array(P_10_Tad)

P_1_T_1000 = pd.read_csv('equlibrium_result_at_T_1000_P_101.325.csv')
P_1_T_1000 = np.array(P_1_T_1000)

P_1_T_1500 = pd.read_csv('equlibrium_result_at_T_1500_P_101.325.csv')
P_1_T_1500 = np.array(P_1_T_1500)

P_1_T_2000 = pd.read_csv('equlibrium_result_at_T_2000_P_101.325.csv')
P_1_T_2000 = np.array(P_1_T_2000)

P_5_T_1000 = pd.read_csv('equlibrium_result_at_T_1000_P_506.625.csv')
P_5_T_1000 = np.array(P_5_T_1000)

P_5_T_1500 = pd.read_csv('equlibrium_result_at_T_1500_P_506.625.csv')
P_5_T_1500 = np.array(P_5_T_1500)

P_5_T_2000 = pd.read_csv('equlibrium_result_at_T_2000_P_506.625.csv')
P_5_T_2000 = np.array(P_5_T_2000)

P_10_T_1000 = pd.read_csv('equlibrium_result_at_T_1000_P_1013.25.csv')
P_10_T_1000 = np.array(P_10_T_1000)

P_10_T_1500 = pd.read_csv('equlibrium_result_at_T_1500_P_1013.25.csv')
P_10_T_1500 = np.array(P_10_T_1500) 

P_10_T_2000 = pd.read_csv('equlibrium_result_at_T_2000_P_1013.25.csv')
P_10_T_2000 = np.array(P_10_T_2000)

#------------------------------------------------------------------------------

fig = plt.figure(1, figsize=(10, 5))

especie = 9
plt.plot(P_1_Tad[1:, 0].astype(float), P_1_Tad[1:, especie+1].astype(float), label = '$X_{'+str(P_1_Tad[0, especie+1])+'}$ - $T_{ad}$')
plt.plot(P_1_T_1000[1:, 0].astype(float), P_1_T_1000[1:, especie].astype(float), label = '$X_{'+str(P_1_T_1000[0, especie])+'}$ - $T=1000$ K')
plt.plot(P_1_T_1500[1:, 0].astype(float), P_1_T_1500[1:, especie].astype(float), label = '$X_{'+str(P_1_T_1500[0, especie])+'}$ - $T=1500$ K')
plt.plot(P_1_T_2000[1:, 0].astype(float), P_1_T_2000[1:, especie].astype(float), label = '$X_{'+str(P_1_T_2000[0, especie])+'}$ - $T=2000$ K')
    
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('Fracción molar en equilibrio a P= 1 atm para diferentes temperaturas')
plt.ylabel('Fracción molar')
plt.xlabel('Radio de equivalencia en la entrada $\phi$')
plt.grid(color='k', linestyle=':')
plt.show()

fig.savefig("nombre_imagen.png", dpi=600)