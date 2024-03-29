# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 14:21:53 2023

@author: seforeros
"""
import numpy as np
import math

Ru = 1.987207 # [cal/mol.K]

constants = {
    '1': np.array([[1.915E+14, 0.00, 1.644E+04],[5.481E+11, 0.39, -2.930E+02]]),
    '2': np.array([[5.080E+04, 2.67, 6.292E+03],[2.667E+04, 2.65, 4.880E+03]]),
    '3': np.array([[2.160E+08, 1.51, 3.430E+03],[2.298E+09, 1.40, 1.832E+04]]),
    '4': np.array([[2.970E+06, 2.02, 1.340E+04],[1.465E+05, 2.11, -2.904E+03]]),
    '5': np.array([[4.577E+19, -1.40, 1.044E+05],[1.146E+20, -1.68, 8.200E+02]]),
    '6': np.array([[4.515E+17, -0.64, 1.189E+05],[6.165E+15, -0.50, 0.000E+00]]),
    '7': np.array([[9.880E+17, -0.74, 1.021E+05],[4.714E+18, 1.00, 0.000E+00]]),
    '8': np.array([[1.912E+23, -1.83, 1.185E+05],[4.500E+22, -2.00, 0.000E+00]]),
    '9': np.array([[1.475E+12, 0.60, 0.000E+00],[3.090E+12, 0.53, 4.887E+04]]),
    '10': np.array([[1.660E+13, 0.00, 8.230E+02],[3.164E+12, 0.35, 5.551E+04]]),
    '11': np.array([[7.079E+13, 0.00, 2.950E+02],[2.027E+10, 0.72, 3.684E+04]]),
    '12': np.array([[3.250E+13, 0.00, 0.000E+00],[3.252E+12, 0.33, 5.328E+04]]),
    '13': np.array([[2.890E+13, 0.00, -4.970E+02],[5.861E+13, 0.24, 6.908E+04]]),
    '14': np.array([[4.634E+16, -0.35, 5.067E+04],[4.200E+14, 0.00, 1.198E+04]]),
    '15': np.array([[2.951E+14, 0.00, 4.843E+04],[3.656E+08, 1.14, -2.584E+03]]),
    '16': np.array([[2.410E+13, 0.00, 3.970E+03],[1.269E+08, 1.31, 7.141E+04]]),
    '17': np.array([[6.025E+13, 0.00, 7.950E+03],[1.041E+11, 0.70, 2.395E+04]]), #18
    '18': np.array([[9.550E+06, 2.00, 3.970E+03],[8.660E+03, 2.68, 1.856E+04]]), #19
    '19': np.array([[1.000E+12, 0.00, 0.000E+00],[1.838E+10, 0.59, 3.089E+04]])
    }
def k_values(T):
    
    '''
    INPUT: Temperature [K]
    
    OUTPUT: 19x2 array that contains the Elementary reaction rate coefficients
            kf (1st column) and kr (2nd column)
    '''
    
    K = np.empty((19, 2))
    
    for key in constants:
        
        forward = constants.get(key)[0]
        reverse = constants.get(key)[1]
        
        kf = forward[0]*(T**forward[1])*math.exp(-forward[2]/(Ru*T))
        kr = reverse[0]*(T**reverse[1])*math.exp(-reverse[2]/(Ru*T))
        
        K[int(key)-1] = [kf,kr]
     
    return K

def get_kc(T, reactions):
    
    K = []
    
    for key in reactions:
        
        forward = constants.get(key)[0]
        reverse = constants.get(key)[1]
        
        kf = forward[0]*(T**forward[1])*math.exp(-forward[2]/(Ru*T))
        kr = reverse[0]*(T**reverse[1])*math.exp(-reverse[2]/(Ru*T))
        
        kc = kf/kr
        K.append(kc)
    
    return K
    

def get_ef():
    correction = {
        '5': np.array([2.5, 12.0]),
        '6': np.array([2.5, 12.0]),
        '7': np.array([2.5, 12.0]),
        '8': np.array([0.73, 12.0]),
        '9': np.array([1.3, 14.0]),
        '15': np.array([2.5, 12.0]),
        }
    
    A = np.ones((6, 9))
    i = 0
    for key in correction:
        
        A[i,1] = correction[key][0]
        A[i,5] = correction[key][1]
        i += 1
    
    return A