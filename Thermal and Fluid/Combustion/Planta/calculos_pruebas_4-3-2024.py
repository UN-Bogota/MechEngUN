# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 21:40:43 2024

@author: Seforeros
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import MultipleLocator


def perform_calc(df):
    
    rho_s, P_s, T_s = 0.6569, 14.696, 25 + 273.15 # Condiciones estándar g/L, psia, K 
    
    R_CH4 = 0.51828 # kJ/kg K
    nc = 2.0 # Number of revolutions per work cycle
    Vd = 406e-6 # Displaced Volume m3
    rho_B10 = 838.00 # kg/m3
    LHV_B10 = 46079.97 # kJ/kg
    LHV_GNV = 42188.30 # kJ/m3
    LHV_GNVm = 42619.0 #kJ/kg

    # Exponential association for the generator y = a(1-exp(-bx))
    a, b = 0.7560299, 0.005690096
        
    df['flujo_masico'] = ((df.m_0-df.m_60)/60 + (df.m_0-df.m_30)/30)/2
    df['potencia_elec'] = (df.Vi*df.Ii + df.Vc*df.Ic + df.Vd*df.Id)/1000 # kW
    df['eficiencia_gen'] = a*(1-np.exp(-b*1000*df.potencia_elec)) # -
    df['potencia_freno'] = df.potencia_elec/df.eficiencia_gen # kW
    df['bmep'] = (df.potencia_freno*nc)/(Vd*df.rpm/60) # kPa
    #df['flujo_vol']=df[['flujo_vol_1', 'flujo_vol_2']].mean(axis=1)#L/min
   
    df['energía_Diesel'] = (df.flujo_masico/1000)*2*LHV_B10/(df.rpm/60)
    if 'flujo_vol_1' in df.keys():    
        df['densidad'] = rho_s*(T_s/(df.T_ingreso+273.15))*(df.P_ingreso/P_s) # kg/m3 Con correccion del manual
    
        df['flujo_m1'] = (df.flujo_vol_1/60000)*df.densidad # kg/s
        df['flujo_m2'] = (df.flujo_vol_2/60000)*df.densidad # kg/s
    else:
        df['densidad'] = rho_s
        df['flujo_m1'] = (df.flujo_M_1/60000)*df.densidad # kg/s
        df['flujo_m2'] = (df.flujo_M_2/60000)*df.densidad # kg/s
    # df['flujo_m1'] = (df.P_ingreso*6.89476)*(df.flujo_vol_1/60000)/(R_CH4 * (df.T_ingreso + 273.15))
    # df['flujo_m2'] = (df.P_ingreso*6.89476)*(df.flujo_vol_2/60000)/(R_CH4 * (df.T_ingreso + 273.15))
    
    df['energía_GNV1'] = (df.flujo_m1)*2*LHV_GNVm/(df.rpm/60)
    df['energía_GNV2'] = (df.flujo_m2)*2*LHV_GNVm/(df.rpm/60)    
    df['energia_calorifica1'] = df['energía_Diesel'] + df['energía_GNV1']
    df['energia_calorifica2'] = df['energía_Diesel'] + df['energía_GNV2']
    df['energia_calorifica_avg'] = 100*df[['energia_calorifica1', 'energia_calorifica2']].mean(axis=1)
    df['energia_calorifica_std'] = 100*df[['energia_calorifica1', 'energia_calorifica2']].std(axis=1)

    df['sustitucion1'] =  df['energía_GNV1']/df['energia_calorifica1']
    df['sustitucion2'] =  df['energía_GNV2']/df['energia_calorifica2']
    df['sustitucion_avg'] = 100*df[['sustitucion1', 'sustitucion2']].mean(axis=1)
    df['sustitucion_std'] = 100*df[['sustitucion1', 'sustitucion2']].std(axis=1)
    
    df['eff_termica1'] = df.bmep*Vd/df.energia_calorifica1
    df['eff_termica2'] = df.bmep*Vd/df.energia_calorifica2
    df['eff_termica_avg'] = 100*df[['eff_termica1', 'eff_termica2']].mean(axis=1)
    df['eff_termica_std'] = 100*df[['eff_termica1', 'eff_termica2']].std(axis=1)

    df['flujo_energetico1'] = df.energia_calorifica1/2*df.rpm #kJ/min
    df['flujo_energetico2'] = df.energia_calorifica2/2*df.rpm #kJ/min
    df['flujo_energetico_avg'] = df[['flujo_energetico1', 'flujo_energetico2']].mean(axis=1)
    df['flujo_energetico_std'] = df[['flujo_energetico1', 'flujo_energetico2']].std(axis=1)
    
    df['consumo_especifico1'] = (df['energía_GNV1']/LHV_GNVm*1000+df['energía_Diesel']/LHV_B10*1000)/(df.bmep*Vd)*3600 #g/kWh
    df['consumo_especifico2'] = (df['energía_GNV2']/LHV_GNVm*1000+df['energía_Diesel']/LHV_B10*1000)/(df.bmep*Vd)*3600 #g/kWh
    df['consumo_especifico_avg'] = df[['consumo_especifico1', 'consumo_especifico2']].mean(axis=1)
    df['consumo_especifico_std'] = df[['consumo_especifico1', 'consumo_especifico2']].std(axis=1)

    return df



# Lectura del df  -------------------------------------------------------------

filename = 'pruebas_lab.xlsx'

df0 = pd.read_excel(filename, 'CH4_2', skiprows= 4, nrows= (23-6+1), usecols='A:S')

#df0 = df0.drop(['m_0', 'm_30', 'm_60'], axis=1)

#df = df0.groupby(['mapa', 'carga']).mean().reset_index()
#df_std = df0.groupby(['mapa', 'carga']).std().reset_index()


df1 = pd.read_excel(filename, 'CH4_1', skiprows= 4, nrows= (41-6+1), usecols='A:Q')
#df0['flujo_masico'] = ((df0.m_0-df0.m_60)/60 + (df0.m_0-df0.m_30)/30 + (df0.m_30-df0.m_60)/30)/3

#df0 = df0.drop(['m_0', 'm_30', 'm_60'], axis=1)

df1_ = df1.groupby(['mapa', 'carga']).mean().reset_index()
#df1_std = df1.groupby(['mapa', 'carga']).std().reset_index()

# Calculos  -------------------------------------------------------------------

df = perform_calc(df0)
df_ = perform_calc(df1_)


dual_df = df[df['mapa'].isin(['Dual 1 CH4', 'Dual 2 CH4'])]
diesel = df[df['mapa'].isin(['Diesel'])]
dual1 = df[df['mapa'].isin(['Dual 1 CH4'])]
dual2 = df[df['mapa'].isin(['Dual 2 CH4'])]

dual_df_ = df_[df_['mapa'].isin(['Dual 1 CH4', 'Dual 2 CH4'])]
diesel_ = df_[df_['mapa'].isin(['Diesel'])]
dual1_ = df_[df_['mapa'].isin(['Dual 1 CH4'])]
dual2_ = df_[df_['mapa'].isin(['Dual 2 CH4'])]


w, h = 10, 6 # tamaño de las figuras

texto = '_____ Prueba 3/3/2023 \n\n _ _ _ Prueba 21/3/2023 '

# Fig1 - Porcentaje de sustitucion --------------------------------------------

fig, ax1 = plt.subplots(figsize=(w, h))
fig1, fig11, fig12 = dual_df, dual1, dual2 # cambiar
fig1_, fig11_, fig12_ = dual_df_, dual1_, dual2_ # cambiar

sns.lineplot(data=fig1, x='carga', y='sustitucion_avg', hue='mapa', palette=['blue', 'red'])
#ax1.fill_between(fig11['carga'], fig11.sustitucion_avg - fig11.sustitucion_std, fig11.sustitucion_avg + fig11.sustitucion_std, alpha=0.2, color='blue')
#ax1.fill_between(fig12['carga'], fig12.sustitucion_avg - fig12.sustitucion_std, fig12.sustitucion_avg + fig12.sustitucion_std, alpha=0.2, color='red')
sns.lineplot(data=fig1_, x='carga', y='sustitucion_avg', hue='mapa', linestyle='dashed', palette=['blue', 'red'], legend=None)

#ax1.text(0.24*w, 3.6*h, texto, fontsize=10, ha='center')
ax1.set(ylabel = 'Sustitución de $CH_4$ [ % ]')
ax2 = ax1.twiny()
sns.lineplot(data=fig1, x='bmep', y='sustitucion_avg', visible=False)
plt.show()


# Fig2 - Eficiencia termica ---------------------------------------------------

fig, ax1 = plt.subplots(figsize=(w, h))
fig1, fig11, fig12 = df, dual1, dual2 # cambiar
fig1_, fig11_, fig12_ = df_, dual1_, dual2_ # cambiar

sns.lineplot(data=fig1, x='carga', y='eff_termica_avg', hue='mapa', palette=['black', 'blue', 'red'])
#ax1.fill_between(fig11['carga'], fig11.eff_termica_avg - fig11.eff_termica_std, fig11.eff_termica_avg + fig11.eff_termica_std, alpha=0.2, color='blue')
#ax1.fill_between(fig12['carga'], fig12.eff_termica_avg - fig12.eff_termica_std, fig12.eff_termica_avg + fig12.eff_termica_std, alpha=0.2, color='red')
sns.lineplot(data=fig1_, x='carga', y='eff_termica_avg', hue='mapa', linestyle='dashed', palette=['black', 'blue', 'red'], legend=None)

#ax1.text(0.24*w, 3.6*h, texto, fontsize=10, ha='center')
ax1.set(ylabel = 'Eficiencia térmica al freno ($\eta_{th,b}$) [ % ]')
ax2 = ax1.twiny()
sns.lineplot(data=fig1, x='bmep', y='eff_termica_avg', visible=False)
plt.show()


# Fig3 - Flujo energetico -----------------------------------------------------

fig, ax1 = plt.subplots(figsize=(w, h))
fig1, fig11, fig12 = df, dual1, dual2 # cambiar
fig1_, fig11_, fig12_ = df_, dual1_, dual2_ # cambiar

sns.lineplot(data=fig1, x='carga', y='flujo_energetico_avg', hue='mapa', palette=['black', 'blue', 'red'])
#ax1.fill_between(fig11['carga'], fig11.flujo_energetico_avg - fig11.flujo_energetico_std, fig11.flujo_energetico_avg + fig11.flujo_energetico_std, alpha=0.2, color='blue')
#ax1.fill_between(fig12['carga'], fig12.flujo_energetico_avg - fig12.flujo_energetico_std, fig12.flujo_energetico_avg + fig12.flujo_energetico_std, alpha=0.2, color='red')
sns.lineplot(data=fig1_, x='carga', y='flujo_energetico_avg', hue='mapa', linestyle='dashed', palette=['black', 'blue', 'red'], legend=None)

#ax1.text(0.24*w, 3.6*h, texto, fontsize=10, ha='center')
ax1.set(ylabel = 'Flujo energético [ kJ / min ]')
ax2 = ax1.twiny()
sns.lineplot(data=fig1, x='bmep', y='flujo_energetico_avg', visible=False)
plt.show()


# Fig4 - Consumo específico de combustible ------------------------------------

fig, ax1 = plt.subplots(figsize=(w, h))
fig1, fig11, fig12 = df, dual1, dual2 # cambiar
fig1_, fig11_, fig12_ = df_, dual1_, dual2_ # cambiar

sns.lineplot(data=fig1, x='carga', y='consumo_especifico_avg', hue='mapa', palette=['black', 'blue', 'red'])
#ax1.fill_between(fig11['carga'], fig11.consumo_especifico_avg - fig11.consumo_especifico_std, fig11.consumo_especifico_avg + fig11.consumo_especifico_std, alpha=0.2, color='blue')
#ax1.fill_between(fig12['carga'], fig12.consumo_especifico_avg - fig12.consumo_especifico_std, fig12.consumo_especifico_avg + fig12.consumo_especifico_std, alpha=0.2, color='red')
sns.lineplot(data=fig1_, x='carga', y='consumo_especifico_avg', hue='mapa', linestyle='dashed', palette=['black', 'blue', 'red'], legend=None)

#ax1.text(2.4, 10, texto, fontsize=10, ha='center')
ax1.set(ylabel = 'Consumo específico de combustible ($sfc$) [ g / kW-h ]')
ax2 = ax1.twiny()
sns.lineplot(data=fig1, x='bmep', y='consumo_especifico_avg', visible=False)
plt.show()


# Fig5 - Potencia al freno ----------------------------------------------------

fig, ax1 = plt.subplots(figsize=(w, h))
fig1 = df # cambiar
fig1_ = df_ # cambiar

sns.lineplot(data=fig1, x='carga', y='potencia_freno', hue='mapa', palette=['black', 'blue', 'red'])
sns.lineplot(data=fig1_, x='carga', y='potencia_freno', hue='mapa', linestyle='dashed', palette=['black', 'blue', 'red'], legend=None)

#ax1.text(2.4, 3.1, texto, fontsize=10, ha='center')
ax1.set(ylabel = 'Potencia al freno ($P_b$) [ kW ]')
ax2 = ax1.twiny()
sns.lineplot(data=fig1, x='bmep', y='potencia_freno', visible=False)
plt.show()

#Caracterización inyector
#sustitución=13.033 ln(t_iny)-1.751
df_tot=pd.concat([df,df_],axis=0)
sns.regplot(data=df_tot[df_tot['mapa']!='Diesel'],x='t_inyeccion',y='sustitucion_avg',logx=True,ci=95)