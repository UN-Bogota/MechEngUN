#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:35:50 2023

@author: bojack
"""

import numpy as np
import kp 
from scipy.optimize import fsolve

class Reaction:
    
    
    """
        The Reaction class defines the object for the chemical reaction of combustion.
        The object has methods to calculate the coefficients for the balanced equation, 
        and to calculate the mole fractions and mass fractions of the products.
        
        Attributes:
        -----------
        n_entrada : numpy array
            Number of moles of each fuel species.
        comp_entrada : numpy array
            Composition of each fuel species in terms of the elements C, H, O, and N.
        n_salida_est : numpy array
            Number of moles of each product species as determined by stoichiometry.
        comp_salida_est : numpy array
            Composition of each product species as determined by stoichiometry.
        n_salida_real : numpy array
            Number of moles of each product species as determined experimentally.
        comp_salida_real : numpy array
            Composition of each product species as determined experimentally.
        A_real : float
            The real stoichiometric coefficient.
        A_est : float
            The estimated stoichiometric coefficient.
        phi : float
            The fuel-to-air equivalence ratio.
        Carbon : float
            The total number of carbon atoms in the fuel.
        Hydrogen : float
            The total number of hydrogen atoms in the fuel.
        Oxygen : float
            The total number of oxygen atoms in the fuel.
        Nitrogen : float
            The total number of nitrogen atoms in the fuel.
        
        Methods:
        --------
        addPhi(phi: float)
            Adds the equivalence ratio to the Reaction object.
        addFuelSpecies(comp_entrada_i: list, n_i: int = 1)
            Adds a new fuel species to the Reaction object.
        getTotalFuelElements()
            Calculates the total number of each element in the fuel.
        getA_est()
            Calculates the estimated stoichiometric coefficient.
        getA_real()
            Calculates the real stoichiometric coefficient.
        getTotalReacRealElements()
            Calculates the total number of each element in the products.
        getEquations(vars: list, T: float = 1000, P: float = 101.325)
            Calculates the coefficients for the balanced equation.
    """
    
    def __init__(self):
        
        C_molarMass = 12.0107 # kg/kmol
        H_molarMass = 1.00794 # kg/kmol
        O_molarMass = 15.9994 # kg/kmol
        N_molarMass = 14.0067 # kg/kmol
        
        genMolarWeight = np.array([C_molarMass, H_molarMass, O_molarMass, N_molarMass])
        
        self.genMolarWeight = genMolarWeight.transpose() # general CHON molar weight
                
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
        
        self.Carbon = None
        self.Hydrogen = None
        self.Oxygen = None
        self.Nitrogen = None
        
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
        
    def getEquations(self, vars, T = 2000, P = 101.325):
        
        """
        Se usan las 4 ecuaciones de balance de masa y 8 de equilibrio para 
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
        if self.Carbon == None:
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
        eqn3 = M**2 * f - kp_list[2]*H  # N2_to_2N
        eqn4 = I**2 - kp_list[3]**2 * (G*H) #  O2-N2_to_NO
        eqn5 = G * E**2 * f - (kp_list[4]*D)**2 # -- H2O_to_H2-O2
        eqn6 = F**2 * E * f - (kp_list[5]*D)**2 # H2O_to_OH-H2
        eqn7 = C**2 * G * f - (kp_list[6]*B)**2 # CO2_to_CO-O2
        eqn8 = C*D - (E*B)*kp_list[7] # CO2-H2_to_CO-H2O
                
        return [C_balance, H_balance, O_balance, N_balance, 
                eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8]
        
    def solveSystem(self, CI):
        
        B, C, D, E, F, G, H, I, J, K, L, M =  fsolve(self.getEquations, CI)
        self.n_salida_real = [B, C, D, E, F, G, H, I, J, K, L, M]
        
        return self.n_salida_real
    
    def getReactmass(self):
        
        fuelMAss = np.matmul(self.n_entrada, (np.matmul(self.comp_entrada, 
                                                        self.genMolarWeight))) #kgFuel
        
        air = np.array([0, 0, self.A_real*0.21*2, self.A_real*0.79*2])
        airMass = np.matmul(air, self.genMolarWeight)
        totalReacMass = fuelMAss + airMass
        
        return totalReacMass
    
    def getProdMass(self):
        
        prodMAss = np.matmul(self.n_salida_real, (np.matmul(self.comp_salida_real, 
                                                        self.genMolarWeight))) #kgFuel
        
        return prodMAss
        
        