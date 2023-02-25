# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 12:51:39 2022

@author: Bojack Horseman
"""
import ThProperties as Th
from matplotlib import pyplot as plt
delta = 1000
epsilon = 10
Li = 300
Lu = 5000
Ti = (Li+Lu)/2
x = 9.1139
while abs(delta) >= epsilon:
    delta = 874830 - ((0.113*(Th.cal_property(Ti, 'CO2', 'h') - Th.cal_property(298.15, 'CO2', 'h')) +
                           0.018*(Th.cal_property(Ti, 'CO', 'h') - Th.cal_property(298.15, 'CO', 'h')) +
                           6E-4*(Th.cal_property(Ti, 'NO', 'h') - Th.cal_property(298.15, 'NO', 'h')) +
                           1.195*(Th.cal_property(Ti, 'N2', 'h') - Th.cal_property(298.15, 'N2', 'h')) +
                           0.08*(Th.cal_property(Ti, 'O2', 'h') - Th.cal_property(298.15, 'O2', 'h')))*x +
                          2.162*(Th.cal_property(Ti, 'H2O', 'h') - Th.cal_property(298.15, 'H2O', 'h')))

    if delta < 0:
        Lu = Ti
    else:
        Li = Ti
        
    Ti = (Li+Lu)/2

print('El delta de parada es: ', delta)
print('La temperatura para la cual se cumple la igualdad es: ', Ti)


# 4. suponga varias eficiencias y encuentre las temperaturas:
eta = [0, 0.15, 0.35, 0.4, 0.6, 0.8, 0.9, 1]

LHV = -862770
Hr = -71263
Hcom_p = -946090
temperatures = []
for i in eta:
    Ti = 300
    delta = 1000
    while delta >= epsilon and Ti <= 5000:
        delta = abs((i*LHV+Hr-Hcom_p) - ((0.113*(Th.cal_property(Ti, 'CO2', 'h') - Th.cal_property(298.15, 'CO2', 'h')) +
                               0.018*(Th.cal_property(Ti, 'CO', 'h') - Th.cal_property(298.15, 'CO', 'h')) +
                               6E-4*(Th.cal_property(Ti, 'NO', 'h') - Th.cal_property(298.15, 'NO', 'h')) +
                               1.195*(Th.cal_property(Ti, 'N2', 'h') - Th.cal_property(298.15, 'N2', 'h')) +
                               0.08*(Th.cal_property(Ti, 'O2', 'h') - Th.cal_property(298.15, 'O2', 'h')))*x +
                              2.162*(Th.cal_property(Ti, 'H2O', 'h') - Th.cal_property(298.15, 'H2O', 'h'))))

        Ti += 0.1
    temperatures.append(Ti)

print('Las temperaturas son: ', temperatures)
print('Para cada eta: ', eta)
plt.plot(eta, temperatures, 'ro')
plt.xlabel('Eficiencia')
plt.ylabel('Temperatura [K]')
plt.title('Eficiencia vs T')
plt.grid()
plt.show()
