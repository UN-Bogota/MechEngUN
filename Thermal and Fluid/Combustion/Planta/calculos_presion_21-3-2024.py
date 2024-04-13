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
        
        #print(csv_file)
        df = pd.read_csv(csv_file, names=col_names, header=0)         
        df['fuel'] = csv_file[14+len(local_storage_directory):16+len(local_storage_directory)]
        df['load'] = int(csv_file[17+len(local_storage_directory):18+len(local_storage_directory)])
        df['rpm'] = int(csv_file[9+len(local_storage_directory):13+len(local_storage_directory)])
        
        combined_df = pd.concat([combined_df, df], axis=0)  
        #print(f'{df}')

    return combined_df

# Start of the code ------------------------------------------------


local_storage_directory = 'Mar21-2024/' #change depending of your storage location
filenames = list_files(local_storage_directory)

dataframes1 = create_df(filenames)
dataframes1['prueba']=1
local_storage_directory = 'Abril 5 2024/'
filenames = list_files(local_storage_directory)

dataframes2 = create_df(filenames)
dataframes2['prueba']=2
dataframes=pd.concat([dataframes1,dataframes2],axis=0).reset_index()
#print(dataframes)

# caracterización de la planta
diametro = 0.086 #m
radio_ciguenal = 0.035  #m
lon_biela = 0.117  #m
carrera = 2*radio_ciguenal
r = 14.58

P_abs=dataframes.groupby(['fuel','load','rpm','prueba']).max()
df_plot = P_abs['P_abs_cam'].to_frame().reset_index()
sns.lineplot(data=df_plot,x='load',y='P_abs_cam',hue='fuel',style='prueba')

#df_plot = P_abs['an'].to_frame().reset_index()
#sns.lineplot(data=df_plot,x='load',y='P_abs_cam',hue='fuel')

sns.relplot(
    data=dataframes2, x="angle", y="P_abs_cam",
    col="load", hue="fuel",
    kind="scatter",col_wrap=2)

#sns.lineplot(data=dataframes,x='angle',y='P_abs_cam',hue='load',col='fuel')
#sns.lineplot(P_abs['B8'],P_abs['D1'])
df_plot['angle']=0
for i in range(len(df_plot)):
    tem=dataframes[(dataframes['fuel']==df_plot['fuel'][i])*(dataframes['load']==df_plot['load'][i])*(dataframes['P_abs_cam']==df_plot['P_abs_cam'][i])]['angle'].to_numpy()
    df_plot['angle'][i]=tem[0]
plt.xlim([-50,100])
# plt.figure()
# sns.lineplot(data=df_plot, x='load', y='angle', hue='fuel')
