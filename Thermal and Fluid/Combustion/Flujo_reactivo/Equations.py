#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 18:15:21 2023

@author: julian
"""
import numpy as np
n=10
dx=L/(n-1)
especies=['','','','','','','']
cells=n*(3*len(especies)+1)
A=np.zeros([cells])
C=np.zeros([cells,1])
Av=np.zeros([len(especies)])
for i in range(especies):
    for j in range(especies):
            Av[i,i]+=-X[i]*X[j]/D[i,j]
            Av[i,j]+=X[i]*X[j]/D[i,j]
    
