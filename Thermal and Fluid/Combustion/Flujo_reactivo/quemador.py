#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 17:43:46 2023

@author: bojack
"""


def MBtu_hr_in2_to_W_mm2(value):
    
    cte_MBtu_to_w = 293071.07017222
    MBtu_hr_in2_to_W_in2 = value*cte_MBtu_to_w
    cte_inch2__to_mm2 = 645.16
    MBtu_hr_in2_to_W_mm2 = MBtu_hr_in2_to_W_in2/cte_inch2__to_mm2
    
    return MBtu_hr_in2_to_W_mm2


h_h2o_25 = 104.89 #[kJ/kg]
h_h2o_100 = 419.04 #[kJ/kg]
LHV_fuel = 50016 # [kJ/kg]
rho_h2o = 997 #[kg/m³]
vol = 2 #L
tiempo = 5 #minutos
eff = 0.3
porcentaje_air_pri = 0.45 # por la tabla 8.25
A_F_estoi = 17.1913


masa = rho_h2o*vol*0.001

Energy = masa*(h_h2o_100-h_h2o_25)

# La energia total teniendo en cuenta la eficiencia:
    
E_total = Energy/eff

# Potencia que se debe entregar: 
    
P = E_total/(tiempo*60) # kW



max_GIR = MBtu_hr_in2_to_W_mm2(12.5)
A_tot = P*1000/max_GIR # [mm2]



massFlux_fuel = P/LHV_fuel # [kg/s]

massFlux_air_pri = porcentaje_air_pri*A_F_estoi*massFlux_fuel


rho_mean = 1.083952

Q_tot = (massFlux_fuel+massFlux_air_pri)/rho_mean

S = (1-porcentaje_air_pri)/(porcentaje_air_pri+(0.21/2))





