import numpy as np

def log10_kp(name, temp):
    """
    This function returns the equilibrium constant, Kp, for a given reaction 
    and temperature using experimental data. The experimental data is stored 
    in a dictionary containing temperature values and corresponding 
    equilibrium constants for several reactions. 
    The function first retrieves the experimental data for the specified 
    reaction and then interpolates the equilibrium constant at the given 
    temperature using numpy's interp function. 
    The resulting log10(Kp) value is returned.

    Parameters
    ----------
    name (str): 
        The name of the reaction for which the equilibrium constant is to be calculated.
    temp (float): 
        The temperature (in K) at which the equilibrium constant is to be calculated.

    Returns
    -------
    y_interp : float 
    The log10(Kp) value for the specified reaction at the given temperature.

    """
    
    cte_data = {
        'temp': np.concatenate(([298, 500, 1000, 1200], np.arange(1600, 3600, 100))),
        'H2_to_2H': np.array([-71.224, -40.316, -17.292, -13.414, -8.532, -7.666, -6.896, -6.204, -5.580, -5.016, -4.502, -4.032, -3.600, -3.202, -2.836, -2.494, -2.178, -1.882, -1.606, -1.348, -1.106, -0.878, -0.664, -0.462]),
        'O2_to_2O': np.array([-81.208, -45.880, -19.614, -15.208, -9.684, -8.706, -7.836, -7.058, -6.356, -5.720, -5.142, -4.614, -4.130, -3.684, -3.272, -2.892, -2.536, -2.206, -1.898, -1.610, -1.340, -1.086, -0.846, -0.620]),
        'N2_to_2N': np.array([-159.600, -92.672, -43.056, -34.754, -24.350, -22.512, -20.874, -19.410, -18.092, -16.898, -15.810, -14.818, -13.908, -13.070, -12.298, -11.580, -10.914, -10.294, -9.716, -9.174, -8.664, -8.186, -7.736, -7.312]),
        'O2-N2_to_NO': np.array([-15.171, -8.783, -4.062, -3.275, -2.290, -2.116, -1.962, -1.823, -1.699, -1.586, -1.484, -1.391, -1.305, -1.227, -1.154, -1.087, -1.025, -0.967, -0.913, -0.863, -0.815, -0.771, -0.729, -0.690]),
        'H2O_to_H2-O2': np.array([-40.048, -22.886, -10.062, -7.899, -5.180, -4.699, -4.270, -3.886, -3.540, -3.227, -2.942, -2.682, -2.443, -2.224, -2.021, -1.833, -1.658, -1.495, -1.343, -1.201, -1.067, -0.942, -0.824, -0.712]),
        'H2O_to_OH-H2': np.array([-46.054, -26.130, -11.280, -8.811, -5.677, -5.124, -4.613, -4.190, -3.776, -3.434, -3.091, -2.809, -2.520, -2.270, -2.038, -1.823, -1.624, -1.438, -1.265, -1.103, -0.951, -0.809, -0.674, -0.547]),
        'CO2_to_CO-O2': np.array([-45.066, -25.025, -10.221, -7.764, -4.706, -4.169, -3.693, -3.267, -2.884, -2.539, -2.226, -1.940, -1.679, -1.440, -1.219, -1.015, -0.825, -0.649, -0.485, -0.332, -0.189, -0.054, 0.071, 0.190]),
        'CO2-H2_to_CO-H2O': np.array([-5.018, -2.139, -0.159, 0.135, 0.474, 0.530, 0.577, 0.619, 0.656, 0.688, 0.716, 0.742, 0.764, 0.784, 0.802, 0.818, 0.833, 0.846, 0.858, 0.869, 0.878, 0.888, 0.895, 0.902])
        }

    x_exp = cte_data.get('temp')
    y_exp = cte_data.get(name)

    y_interp = np.interp(temp, x_exp, y_exp)
    
    return y_interp

def kp_values(temp, reactions = []):
    
    if len(reactions) == 0:
        reactions = ['H2_to_2H', 'O2_to_2O', 'N2_to_2N', 'O2-N2_to_NO', 'H2O_to_H2-O2', 'H2O_to_OH-H2', 'CO2_to_CO-O2',
                     'CO2-H2_to_CO-H2O']
    
    kp_val = []

    for i in reactions:
        kp = log10_kp(i, temp)
        kp_val.append(np.power(10, kp))
        
    return kp_val