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

