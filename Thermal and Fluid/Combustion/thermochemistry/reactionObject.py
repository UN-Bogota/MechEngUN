#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:35:50 2023

@author: bojack
"""

import numpy as np

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
                                          [0, 0, 0, 1], # N
                                          [0, 0, 0, 2], # N2
                                          [0, 0, 2, 1] # NO2
                                          ])
        self.A_real = 1
        self.A_est = 1
        self.phi = self.A_est/self.A_real
        
        self.realCarbon = 0
        self.realHydrogen = 0
        self.realOxygen = 0 
        self.realNitrogen = 0
        
    def addPhi(self, phi):
        self.phi = phi
    
    def addFuelSpecies(self, comp_entrada_i, n_i = 1):
        
        # [1, 4, 0, 0]
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
        
        
        
        
        