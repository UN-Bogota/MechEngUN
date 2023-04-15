# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:18:54 2023

@author: seforeros
"""
import numpy as np

class kinetics():
    def __init__(self, nup, nupp, k):
        
        self.nup = nup
        self.nupp = nupp
        self.k = k
        self.reactorType = None
        self.effiency = None
        
    def get_M(self, concentrations):
        """Recibe eficiencias (lista longitud especiaes) y concentraciones (incógnitas)
        Retorna lista con las M
        """
        cT=concentrations.T
        M = self.effiency@cT
        return M.T
    
    def setReactorType(self, reactotType):
        self.reactorType = reactotType
        
    def setEfficiency(self, eficiencia):
        
        self.effiency = eficiencia

    def getDiffEq(self,t, concentrations,omega):
        
        concentrations = concentrations.reshape(1, len(concentrations))
        
        con_M = self.get_M(concentrations)
        
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
        omega = q.T @ self.nu
        
        omega = omega[0,0:9]
        print(omega)
        #return omega #np array
    
    
    def getEcuations(self,t,concentrations,omega):
        
        if self.reactorType == None:
            
            return self.getDiffEq(self,t,concentrations,omega)
            
        elif self.reactorType == 'constant-volume':
            
            omega = self.getDiffEq(self,t,concentrations,omega)
            pass
            # acoplat la ecuacion de temperatura a la ecuación de las omega
            
    # Anadir una fucnion que calcule las condiciones iniciales en función del phi
    
            
            