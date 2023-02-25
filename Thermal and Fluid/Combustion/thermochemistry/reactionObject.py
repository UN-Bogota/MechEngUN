#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:35:50 2023

@author: bojack
"""

import numpy as np
import kp 


class Reaction:
    
    
    def __init__(self):
        
        # Parametros de entrada para crear el objeto
        self.n_entrada = np.array([])
        self.comp_entrada = np.array([[0, 0, 0, 0]])
        
        # número de moles y composicion de los productos de la reac. esteq.
        self.n_salida_est = np.array([])
        self.comp_salida_est = np.array([[1, 0, 2, 0],
                                         [0, 2, 1, 0],
                                         [0, 0, 0, 2]])
        
        
        # C H O N  -----> Elementos que componen los productos reales
        self.n_salida_real = np.array([])# C  H  O  N 
        self.comp_salida_real = np.array([
                                          [1, 0, 2, 0], # CO2
                                          [1, 0, 1, 0], # CO
                                          [0, 2, 1, 0], # H2O 
                                          [0, 2, 0, 0], # H2
                                          [0, 1, 1, 0], # OH
                                          [0, 0, 2, 0], # O2
                                          [0, 0, 0, 2], # N2
                                          [0, 0, 1, 1], # NO
                                          [3, 8, 0, 0], # C3H8
                                          [0, 1, 0, 0], # H
                                          [0, 0, 1, 0], # O
                                          [0, 0, 0, 1]  # N
                                          ])
        self.A_real = 1
        self.A_est = 1
        self.phi = self.A_est/self.A_real
        
        self.Carbon = 0
        self.Hydrogen = 0
        self.Oxygen = 0 
        self.Nitrogen = 0
        
    def addPhi(self, phi):
        self.phi = phi
    
    def addFuelSpecies(self, comp_entrada_i, n_i = 1):
        
        self.n_entrada = np.append(self.n_entrada , n_i)
        self.comp_entrada = np.append(self.comp_entrada , [comp_entrada_i], axis = 0)
        
        if all(self.comp_entrada[0] == np.array([0, 0, 0, 0])): 
            self.comp_entrada = np.delete(self.comp_entrada, 0, 0)
            
    def getTotalFuelElements(self):
        Carbon = 0
        Hydrogen = 0
        Oxygen = 0
        Nitrogen = 0
        for i in range(len(self.comp_entrada)):
            Carbon += self.n_entrada[i] * self.comp_entrada[i][0]
            Hydrogen += self.n_entrada[i] * self.comp_entrada[i][1]
            Oxygen += self.n_entrada[i] * self.comp_entrada[i][2]
            Nitrogen += self.n_entrada[i] * self.comp_entrada[i][3]
            
        return Carbon, Hydrogen, Oxygen, Nitrogen

    def getA_est(self):
        
        Carbon, Hydrogen, Oxygen, Nitrogen = self.getTotalFuelElements()
        self.A_est = ((Carbon*2 + Hydrogen/2) - Oxygen)/(2*0.21) # Balance de Oxigenos
        
        return self.A_est

    def getA_real(self):
        
        A_est = self.getA_est()
        self.A_real = A_est/self.phi
        
        return self.A_real
    
    def getTotalReacRealElements(self):
        A_real = self.getA_real()
        self.Carbon, self.Hydrogen, self.Oxygen, self.Nitrogen = self.getTotalFuelElements()
        self.Oxygen +=  A_real*2*0.21
        self.Nitrogen += A_real*0.79*2
        
        return self.Carbon, self.Hydrogen, self.Oxygen, self.Nitrogen
        
    def getEquations(self, vars, T, P = 101.325):
        
        """
        Se usan las 4 ecuaciones de balance de masa y 4 de equilibrio para 
        encontrar cada coeficiente, los coeficientes corresponden a:
                    
            B: CO2
            C: CO
            D: H2O
            E: H2
            F: OH
            G: O2
            H: N2
            I: NO
            J: C3H8
            K: H
            L: O
            M: N
            falta el NO2
        """
        
        B, C, D, E, F, G, H, I, J, K, L, M = vars
        n = B + C + D + E + F + G + H + I + J + K + L + M 
        
        f = (P/101.325)/n
        
        # Todas se igualan a 0
        if self.Carbon == 0:
            self.getTotalReacRealElements()
            
        C_balance = B + C + 3*J - self.Carbon
        H_balance = 2*D + 2*E + F + 8*J + K - self.Hydrogen
        O_balance = 2*B + C + D + F + 2*G + I + L - self.Oxygen 
        N_balance = 2*H + I + M - self.Nitrogen
        
        
        # ['H2_to_2H', 'O2_to_2O', 'N2_to_2N', 'O2-N2_to_NO', 'H2O_to_H2-O2', 
        #  'H2O_to_OH-H2', 'CO2_to_CO-O2','CO2-H2_to_CO-H2O']
        
        kp_list = kp.kp_values(T)
        
        eqn1 = (K**2) * f - kp_list[0]*E # H2_to_2H
        eqn2 = L**2 * f - kp_list[1]*G # O2_to_2O
        eqn3 = M**2 *f - kp_list[2]*H  # N2_to_2N
        eqn4 = I**2 - kp_list[3]**2 *(G*H) #  O2-N2_to_NO
        eqn5 = G*F**2*f - (kp_list[4]*D)**2
        eqn6 = F**2 * E * f - (kp_list[5]*D)**2
        eqn7 = C**2 * G * f - (kp_list[6]*B)**2
        eqn8 = C*D/(E*B)
        
        
        
        return [C_balance, H_balance, O_balance, N_balance, 
                eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8]
        
        
        
        