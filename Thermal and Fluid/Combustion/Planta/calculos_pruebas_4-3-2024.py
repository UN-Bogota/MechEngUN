# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 21:40:43 2024

@author: Seforeros
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class ExcelReader:
    def __init__(self, filename, sheet_name, start_row, end_row, columns_to_read):
        self.filename = filename
        self.sheet_name = sheet_name
        self.start_row = start_row
        self.end_row = end_row
        self.columns_to_read = columns_to_read

    def read_data(self):
        try:
            self.df = pd.read_excel(self.filename, self.sheet_name, skiprows=self.start_row, nrows=self.end_row - self.start_row + 1, usecols=self.columns_to_read)
            return self.df
        except FileNotFoundError:
            print(f"Error: File '{self.filename}' not found.")
            return None
        except Exception as e:
            print(f"Error: {e}")
            return None
        
# Create a class to represent the DataFrame as an object
class df_to_object:
    def __init__(self, df):
        self.rpm = df['rpm'] #rpm
        self.Vi = df['Vi'] #V
        self.Vc = df['Vc'] #V
        self.Vd = df['Vd'] #V
        self.Ii = df['Ii'] #A
        self.Ic = df['Ic'] #A
        self.Id = df['Id'] #A
        self.T_ingreso = df['T_ingreso'] #°C
        self.P_ingreso = df['P_ingreso'] #psi
        self.flujo_vol_1 = df['flujo_vol_1'] #lpm
        self.flujo_vol_2 = df['flujo_vol_2'] #lpm
        self.t_inyeccion = df['t_inyeccion'] #ms
        self.flujo_masico = df['flujo_masico'] #g/s

def perform_calc(df, dual):
    """
    Args:
        df (DataFrame): Containing the data related.
        dual (Boolean): True if the process was released with dual fuel, False in other case.
    """

    nc = 2.0 # Number of revolutions per work cycle
    Vd = 406e-6 # Displaced Volume m3
    rho_B10 = 838.00 # kg/m3
    LHV_B10 = 46079.97 # kJ/kg
    LHV_GNV = 42188.30 # kJ/m3

    # Exponential association for the generator y = a(1-exp(-bx))
    a, b = 0.7560299, 0.005690096

    df.potencia_elec = (df.Vi*df.Ii + df.Vc*df.Ic + df.Vd*df.Id)/1000 # kW
    df.eficiencia_gen = a*(1-np.exp(-b*1000*df.potencia_elec)) # -
    df.potencia_freno = df.potencia_elec/df.eficiencia_gen # kW
    df.bmep = (df.potencia_freno*nc)/(Vd*df.rpm/60) # kPa

    if dual:
        df.energia_calorifica1 = (df.flujo_vol_1/1000/60*df.t_inyeccion/1000*LHV_GNV) + (df.flujo_masico/1000)*2*LHV_B10/(df.rpm/60)

        df.energia_calorifica2 = (df.flujo_vol_2/1000/60*df.t_inyeccion/1000*LHV_GNV) + (df.flujo_masico/1000)*2*LHV_B10/(df.rpm/60)

        df.sustitucion1 = (df.flujo_vol_1/1000/60*df.t_inyeccion/1000*LHV_GNV)/df.energia_calorifica1
        df.sustitucion2 = (df.flujo_vol_2/1000/60*df.t_inyeccion/1000*LHV_GNV)/df.energia_calorifica2
    else: 
        df.energia_calorifica1 = (df.flujo_masico/1000)*2*LHV_B10/(df.rpm/60)

        df.energia_calorifica2 = (df.flujo_masico/1000)*2*LHV_B10/(df.rpm/60)

    df.eff_termica1 = df.bmep*Vd/df.energia_calorifica1
    df.eff_termica2 = df.bmep*Vd/df.energia_calorifica2

    return df

# Calculos --------------------------------------------------------------------

filename = 'pruebas_lab.xlsx'
#filename = '/workspaces/MechEngUN/Thermal and Fluid/Combustion/Planta/pruebas_lab.xlsx'

diesel = df_to_object(ExcelReader(filename, 'CH4_1', 43, 48, 'C:O').read_data())
dual1 = df_to_object(ExcelReader(filename, 'CH4_1', 50, 55, 'C:O').read_data())
dual2 = df_to_object(ExcelReader(filename, 'CH4_1', 57, 62, 'C:O').read_data())

diesel = perform_calc(diesel, False)
dual1 = perform_calc(dual1, True)
dual2 = perform_calc(dual2, True)

# Graficas

    # Eficiencia termica ------------------------------------------------------

plt.figure(figsize=(8*1.5, 6*1.5))

plt.plot(diesel.bmep, 100*diesel.eff_termica1, color = 'black', label='Diesel')

dual1_eff_avg = 100*(dual1.eff_termica1 + dual1.eff_termica2)/2
dual2_eff_avg = 100*(dual2.eff_termica1 + dual2.eff_termica2)/2

plt.plot(dual1.bmep, dual1_eff_avg, color = 'blue', label='Dual1')
plt.fill_between(dual1.bmep, 100*dual1.eff_termica1, 100*dual1.eff_termica2,
                 color='blue', alpha=0.2, label=None)

plt.plot(dual2.bmep, dual2_eff_avg, color = 'red', label='Dual2')
plt.fill_between(dual2.bmep, 100*dual2.eff_termica1, 100*dual2.eff_termica2,
                 color='red', alpha=0.2, label=None)

plt.xlabel('bmep [kPa]')
plt.ylabel('Eficiencia térmica [%]')

plt.legend(loc=2)
plt.grid(True)  # Add grid lines

# Show the plot
plt.show(block = True)


    # Sustitucion -------------------------------------------------------------
    
plt.figure(figsize=(8*1.5, 6*1.5))

dual1_sust_avg = 100*(dual1.sustitucion1 + dual1.sustitucion2)/2
dual2_sust_avg = 100*(dual2.sustitucion1 + dual2.sustitucion2)/2


plt.plot(dual1.bmep, dual1_sust_avg, color = 'blue', label='Dual1')
plt.fill_between(dual1.bmep, 100*dual1.sustitucion1, 100*dual1.sustitucion2,
                 color='blue', alpha=0.2, label=None)

plt.plot(dual2.bmep, dual2_sust_avg, color = 'red', label='Dual2')
plt.fill_between(dual2.bmep, 100*dual2.sustitucion1, 100*dual2.sustitucion2,
                 color='red', alpha=0.2, label=None)



plt.xlabel('bmep [kPa]')
plt.ylabel('Porcentaje de sustitución energética del Metano [%]')
#plt.ylim([0,50])
plt.legend(loc=2)
plt.grid(True)  # Add grid lines


plt.show()