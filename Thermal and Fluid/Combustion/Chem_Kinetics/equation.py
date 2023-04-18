# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:18:54 2023

@author: seforeros
"""
import numpy as np
from resources import *
import properties as prop
import ThProperties as Th
class kinetics():
    def __init__(self, nup, nupp):
        
        self.nup = nup
        self.nupp = nupp
        #self.k = k
        self.reactorType = None
        self.effiency = None
        
    def get_M(self, concentrations):
        """Recibe eficiencias (lista longitud especiaes) y concentraciones (incógnitas)
        Retorna lista con las M
        """
        cT=concentrations.T
        M = self.effiency@cT
        return M.T
    def setT(self,T):
        self.T=T
    def setReactorType(self, reactotType):
        self.reactorType = reactotType
        
    def setEfficiency(self, eficiencia):
        
        self.effiency = eficiencia

    def getDiffEq(self,t, concentrations, omega):
        lab=['H','H2', 'O', 'O2', 'OH' ,'H2O', 'N2', 'HO2', 'H2O2']
        if self.reactorType == 'constant-volume':
            T=[concentrations[-2],concentrations[-1]]
            self.T=T[0]
            concentrations=concentrations[0:-2]
        concentrations = concentrations/sum(concentrations)
        concentrations = concentrations.reshape(1, len(concentrations))
        
        con_M = self.get_M(concentrations)
        #print(con_M[0])
        Mr=[con_M[0][4],con_M[0][5]]
        self.k=k_values(self.T,Mr)
        
        R=8.3144
        concentrations = np.concatenate((concentrations, con_M), axis=1)
        reac, conc = self.nup.shape
        
        q = np.array([[0]])
        
        for i in range(reac):
            prod_f = concentrations  ** self.nup[i, :]
            prod_r = concentrations  ** self.nupp[i, :]
            
            prod_f = np.prod(prod_f)
            prod_r = np.prod(prod_r)
            
            q_i = self.k[i, 0] * prod_f - self.k[i, 1] * prod_r
            
            q = np.append(q, [[q_i]], axis = 0)
        q = np.delete(q, 0, 0)
        self.nu = self.nupp-self.nup
        omega1 = q.T @ self.nu
        omega2=omega1[0]
        #omega = omega[0,0:9]
        for i in range(9):
            omega[i]=omega2[i]
        if self.reactorType == 'constant-volume':
            #print(len(omega),len(prop.hf_reactivos(lab)))
            h_i = np.dot(omega[0:9], prop.hf_reactivos(lab)) # Entalpía de combustion reactivos
            cp_v=np.array([Th.cal_property(self.T, i, 'cp') for i in lab])
            #print(len(concentrations[0]),len(cp_v))
            omega[9]=(R*T[0]*sum(omega[0:9])-h_i)/(np.dot(concentrations[0][0:9],(cp_v-R)))
            omega[10]=R*T[0]*sum(omega[0:9])+R*omega[9]*sum(concentrations[0][0:9])
        #print(omega)
        #return omega #np array
    
    
    def getEcuations(self,t,concentrations,omega):
        
        if self.reactorType == None:
            
            return self.getDiffEq(self,t,concentrations,omega)
            
        elif self.reactorType == 'constant-volume':
            
            omega = self.getDiffEq(self,t,concentrations,omega)
            pass
            # acoplat la ecuacion de temperatura a la ecuación de las omega
            
    # Anadir una fucnion que calcule las condiciones iniciales en función del phi
    
            
            