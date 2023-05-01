#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:35:50 2023

@author: bojack
"""

import numpy as np
import kp 
import kc
from scipy.optimize import fsolve
import properties as prop

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
        self.reactNames = []
        self.productNames = ['CO2', 'CO', 'H2O', 'H2', 'OH', 'O2', 'N2', 'NO', 'C3H8', 'H', 'O', 'N']
        self.solveWithEnergy = False
        self.prodTemperature = None
        self.prodPressure = None
        
    def addComp_salida_est(self, compMatrix):
        self.comp_salida_est = compMatrix
    
    def addComp_salida_real(self, compMatrix, names):
        self.comp_salida_real = compMatrix
        self.productNames = names
        
    def addPhi(self, phi):
        self.phi = phi
        
    def addFuelSpecies(self, comp_entrada_i, name, n_i = 1):
        
        self.n_entrada = np.append(self.n_entrada , n_i)
        self.comp_entrada = np.append(self.comp_entrada , [comp_entrada_i], axis = 0)
        self.reactNames.append(name)
        
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
        
    def addFirstLaw(self, boolean):
        self.solveWithEnergy = boolean
    #def getEquations(self, vars, solveWithEnergy = True, P = 101.325, Patm = 101.325):
        
    def addProductTemperature(self, temp):
        self.prodTemperature = temp
    
    def addProductPressure(self, P):
        self.prodPressure = P
        
    def getEquations(self, vars, Patm = 101.325):
        
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
        
        if self.solveWithEnergy:
            B, C, D, E, F, G, H, I, J, K, L, M, T = vars
        else:
            B, C, D, E, F, G, H, I, J, K, L, M = vars
        
        n = B + C + D + E + F + G + H + I + J + K + L + M  
        f = (self.prodPressure/Patm)/n
        
        if self.Carbon == None:
            self.getTotalReacRealElements()
        
        # Todas se igualan a 0
        
        # Mass balance:
        C_balance = B + C + 3*J - self.Carbon
        H_balance = 2*D + 2*E + F + 8*J + K - self.Hydrogen
        O_balance = 2*B + C + D + F + 2*G + I + L - self.Oxygen 
        N_balance = 2*H + I + M - self.Nitrogen
        
        # There are all the ecuations names, please select all that apply.
        
        # ['H2_to_2H', 'O2_to_2O', 'N2_to_2N', 'O2-N2_to_NO', 'H2O_to_H2-O2', 
        #  'H2O_to_OH-H2', 'CO2_to_CO-O2','CO2-H2_to_CO-H2O']
        
        try:
            kp_list = kp.kp_values(T)
        except:
            if self.prodTemperature == None:
                raise ValueError("Please add the product temperature")
                
            else:
                T = self.prodTemperature
                kp_list = kp.kp_values(T)
            
        
        # Equilibrium equations
        
        eqn1 = (K**2) * f - kp_list[0]*E # H2_to_2H
        eqn2 = L**2 * f - kp_list[1]*G # O2_to_2O
        eqn3 = M**2 * f - kp_list[2]*H  # N2_to_2N
        eqn4 = I**2 - kp_list[3]**2 * (G*H) #  O2-N2_to_NO
        eqn5 = G * E**2 * f - (kp_list[4]*D)**2 # -- H2O_to_H2-O2
        eqn6 = F**2 * E * f - (kp_list[5]*D)**2 # H2O_to_OH-H2
        eqn7 = C**2 * G * f - (kp_list[6]*B)**2 # CO2_to_CO-O2
        eqn8 = C*D - (E*B)*kp_list[7] # CO2-H2_to_CO-H2O
        
        # Energy equations
        
        if self.solveWithEnergy:
            
            n_salida = np.array([B, C, D, E, F, G, H, I, J, K, L, M])
            
            h_com_reac = np.dot(self.n_entrada, prop.hf_reactivos(self.reactNames)) # Entalpía de combustion reactivos
            h_com_pro = np.dot(n_salida, prop.hf_productos(self.productNames)) #Entalpía de combustión productos
                        
            # Falta agregarle el delta de entalpía a los reactivos para cuando ingresan al 
            # reactor a una T diferente a t ambiente.
        
            # The left hand side:
                
            LHS = (h_com_reac - h_com_pro)
            
            # El termino de calbio de la entalpía prod por la deferencia de T:
            
            # Garantizar que los nombres coincidan con los n
            deltaH_pro = np.dot(n_salida, prop.deltaH(T, self.productNames))
            
            energyEq = LHS - deltaH_pro
            
            return [C_balance, H_balance, O_balance, N_balance, 
                    eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, energyEq]
        else:
            return [C_balance, H_balance, O_balance, N_balance, 
                    eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8]
        
    def getH2Equations(self, vars, Patm = 101.325):
        
        """
        Se usan las 4 ecuaciones de balance de masa y 8 de equilibrio para 
        encontrar cada coeficiente, los coeficientes corresponden a:
                    
            B: H
            C: H2
            D: O
            E: O2
            F: OH
            G: H2O
            H: HO2
            I: H2O2
            J: N2
            
        """
        
        if self.solveWithEnergy:
            B, C, D, E, F, G, H, I, J, T = vars
        else:
            B, C, D, E, F, G, H, I, J = vars
        
        n = B + C + D + E + F + G + H + I + J
        
        f = (self.prodPressure/Patm)/n
        
        # Todas las ecuaciones se igualan a 0
        
        # Mass balance:
            
        self.getTotalReacRealElements()
        
        H_balance = B + 2*C + F + 2*G + H + 2*I - self.Hydrogen
        O_balance = D + 2*E + F + G + 2*H + 2*I - self.Oxygen 
        N_balance = 2*J - self.Nitrogen
        
        # There are all the ecuations names, please select all that apply.
        
        # ['H2_to_2H', 'O2_to_2O', 'N2_to_2N', 'O2-N2_to_NO', 'H2O_to_H2-O2', 
        #  'H2O_to_OH-H2', 'CO2_to_CO-O2','CO2-H2_to_CO-H2O']
        
        reactions = ['H2_to_2H', 'O2_to_2O', 'H2O_to_H2-O2', 'H2O_to_OH-H2']
        
        #newReactions = ['17', '18'] # Del mecanismo de reacción tomamos la eq 18 y 19
        #newReactions = ['13', '16'] # Del mecanismo de reacción tomamos la eq 13 y 17
        newReactions = ['13', '17'] # Del mecanismo de reacción tomamos la eq 13 y 18
        try:
            kp_list = kp.kp_values(T, reactions)
            kc_list = kc.get_kc(T, newReactions)
            
        except:
            if self.prodTemperature == None:
                raise ValueError("Please add the product temperature")
                
            else:
                T = self.prodTemperature
                kp_list = kp.kp_values(T, reactions)
                kc_list = kc.get_kc(T, newReactions)
            
        
        # Equilibrium equations
        
        eqn1 = (B**2) * f - kp_list[0]*C # H2_to_2H
        eqn2 = D**2 * f - kp_list[1]*E # O2_to_2O
        eqn3 = E * (C**2) * f - (kp_list[2]*G)**2 # -- H2O_to_H2-O2
        eqn4 = F**2 * C * f - (kp_list[3]*G)**2 # H2O_to_OH-H2
        
        # Cinetic Equations
        
        # eqn5 = C * H - (kc_list[0] * I * B) # H2O2-H_to_H2-HO2
        # eqn6 = F * H - (kc_list[1] * I * D) # H2O2-O_to_OH-HO2
        
        eqn5 = G * E - (kc_list[0] * H * F) 
        eqn6 = C * H - (kc_list[1] * I * B) 
        # Energy equations
        
        if self.solveWithEnergy:
            n_salida = np.array([B, C, D, E, F, G, H, I, J])
            
            h_com_reac = np.dot(self.n_entrada, prop.hf_reactivos(self.reactNames)) # Entalpía de combustion reactivos
            h_com_pro = np.dot(n_salida, prop.hf_productos(self.productNames)) #Entalpía de combustión productos
            
            # Falta agregarle el delta de entalpía a los reactivos para cuando ingresan al 
            # reactor a una T diferente a t ambiente.
        
            # The left hand side:
            deltaH_reac = np.dot(self.n_entrada, prop.deltaH(298.15, self.reactNames))    
            
            LHS = (h_com_reac + deltaH_reac - h_com_pro)
            
            # El termino de calbio de la entalpía prod por la deferencia de T:
            
            # Garantizar que los nombres coincidan con los n
            deltaH_pro = np.dot(n_salida, prop.deltaH(T, self.productNames))
            energyEq = LHS - deltaH_pro
            
            return [H_balance, O_balance, N_balance, 
                    eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, energyEq]
        else:
            return [H_balance, O_balance, N_balance, 
                    eqn1, eqn2, eqn3, eqn4, eqn5, eqn6]
        
    def solveSystem(self, CI):
        
        if self.solveWithEnergy:
            B, C, D, E, F, G, H, I, J, K, L, M, T =  fsolve(self.getEquations, CI)
            self.n_salida_real = [B, C, D, E, F, G, H, I, J, K, L, M]
            return [B, C, D, E, F, G, H, I, J, K, L, M, T]
            
        else:
            B, C, D, E, F, G, H, I, J, K, L, M =  fsolve(self.getEquations, CI)
            self.n_salida_real = [B, C, D, E, F, G, H, I, J, K, L, M]
        
            return self.n_salida_real
        
    def solveH2System(self, CI):
        
        if self.solveWithEnergy:
            B, C, D, E, F, G, H, I, J, T =  fsolve(self.getH2Equations, CI, maxfev=10000)
            self.n_salida_real = [B, C, D, E, F, G, H, I, J]
            return [B, C, D, E, F, G, H, I, J, T]
            
        else:
            B, C, D, E, F, G, H, I, J =  fsolve(self.getH2Equations, CI, maxfev=10000)
            self.n_salida_real = [B, C, D, E, F, G, H, I, J]
        
            return self.n_salida_real
    
    def getReactmass(self):
        
        fuelMAss = np.matmul(self.n_entrada, (np.matmul(self.comp_entrada, 
                                                        self.genMolarWeight))) #kgFuel
        
        air = np.array([0, 0, self.A_real*0.21*2, self.A_real*0.79*2])
        airMass = np.matmul(air, self.genMolarWeight)
        totalReacMass = fuelMAss + airMass
        
        return totalReacMass
    
    def getProdMass(self):
        
        if self.solveWithEnergy:
            prodMAss = np.matmul(self.n_salida_real, (np.matmul(self.comp_salida_real, 
                                                                  self.genMolarWeight))) #kgFuel
        else:
            prodMAss = np.matmul(self.n_salida_real, (np.matmul(self.comp_salida_real, 
                                                                  self.genMolarWeight))) #kgFuel
        return prodMAss
        
        
        