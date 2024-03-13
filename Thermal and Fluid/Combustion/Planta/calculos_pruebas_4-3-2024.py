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
        dual (Boolean): True if the process was releades with dual fuel, False in other case.
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
        df.energia_calorifica1 = (df.flujo_vol_1/1000/60*df.t_iny/1000*LHV_GNV) + (df.flujo_masico/1000)*2*LHV_B10/(df.rpm/60)

        df.energia_calorifica2 = (df.flujo_vol_2/1000/60*df.t_iny/1000*LHV_GNV) + (df.flujo_masico/1000)*2*LHV_B10/(df.rpm/60)

        df.sustitucion1 = dual*(df.flujo_vol_1/1000/60*df.t_iny/1000*LHV_GNV)/df.energia_calorifica1
        df.sustitucion2 = dual*(df.flujo_vol_2/1000/60*df.t_iny/1000*LHV_GNV)/df.energia_calorifica2
    else: 
        df.energia_calorifica1 = (df.flujo_masico/1000)*2*LHV_B10/(df.rpm/60)

        df.energia_calorifica2 = (df.flujo_masico/1000)*2*LHV_B10/(df.rpm/60)

    df.eff_termica1 = df.bmep/df.energia_calorifica1
    df.eff_termica2 = df.bmep/df.energia_calorifica2
    return df


filename = 'pruebas_lab.xlsx'

diesel = df_to_object(ExcelReader(filename, 'CH4_1', 43, 48, 'C:O').read_data())
dual1 = df_to_object(ExcelReader(filename, 'CH4_1', 50, 55, 'C:O').read_data())
dual2 = df_to_object(ExcelReader(filename, 'CH4_1', 57, 62, 'C:O').read_data())

diesel1 = perform_calc(diesel, False)