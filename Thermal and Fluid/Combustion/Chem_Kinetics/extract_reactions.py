# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 15:01:33 2023

@author: seforeros
"""
import re
import numpy as np

constants = {}

with open('h2_v1b_mech.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith('   '):
            reaction = re.findall(r'\w+\+\w+.*\w+\s+\d+\.\d+E[+-]\d+\s+\d+\.\d+', line)
            if reaction:
                reaction = reaction[0].split()
                constants[reaction[0]] = {'A': reaction[2], 'Ea': reaction[3]}

print(constants)


