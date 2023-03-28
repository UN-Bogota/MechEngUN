# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:18:54 2023

@author: seforeros
"""
import numpy as np

class kinetics():
    def __init__(self, nup, nupp, k):
        
        pass
    
    def get_M(self, concentrations):
        """Recibe eficiencias (lista longitud especiaes) y concentraciones (incógnitas)
        Retorna lista con las M
        """
        pass

    def getDiffEq(self, concentrarions):
        
        reac, conc = self.nup.shape
        
        q = np.array([])
        for i in range(reac):
            prod_f = concentrarions  ** self.nup[i, :]
            prod_r = concentrarions  ** self.nupp[i, :]
            
            prod_f = np.prod(prod_f)
            prod_r = np.prod(prod_r)
            
            q_i = self.k[i, 0] * prod_f - self.k[i, 1] * prod_r
            
            # Hacer el append a q de q_i
            
            