#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 18:05:43 2023

@author: bojack
"""

from reactionObject import Reaction

# Creo el objeto
metanoReac = Reaction()

# Se definen las especies del combustible



methane = [1, 4, 0, 0]
phi = 0.9
metanoReac.addFuelSpecies(methane)
metanoReac.addPhi(phi)


CI = [1, 1, 2, 1, 0.1, 0.1, 1, 8, 1, 1, 1, 1]

metanoReac.solveSystem(CI)