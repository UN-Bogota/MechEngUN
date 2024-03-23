import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def list_files(folder_path):
    try:
        filenames = []
        for dirpath, dirnames, files in os.walk(folder_path):
            for name in files:
            
                filenames.append(os.path.join(dirpath, name))
            #print(files)
            break  # Only process the top-level directory
            
        return filenames
    except FileNotFoundError:
        return []  # Handle the case where the folder doesn't exist


def create_df(folder_list):
    # Create an empty dataframe to store the combined data
    combined_df = pd.DataFrame()
    
    col_names = np.array(['P_abs_adm', 'PMS', 'P_abs_cam', 'angle', 'volume'])
       
    for csv_file in folder_list:
        #print(f'file --- {csv_file[-8:]}')

        df = pd.read_csv(csv_file, names=col_names, header=0).T         
        df['fuel'] = csv_file[-8:-6]
        df['load'] = int(csv_file[-5:-4])
        df['rpm'] = int(csv_file[-13:-9])
        
        combined_df = pd.concat([combined_df, df], axis=0)  
        #print(f'{df}')

    return combined_df

# Start of the code ------------------------------------------------

local_storage_directory = 'Mar21-2024/' #change depending of your storage location
filenames = list_files(local_storage_directory)

dataframes = create_df(filenames)
#print(dataframes)

# caracterización de la planta
diametro = 0.086 #m
radio_ciguenal = 0.035  #m
lon_biela = 0.117  #m
carrera = 2*radio_ciguenal
r = 14.58


df['Temperatura']=df['volumen']*df['P_abs']/1000/(287)