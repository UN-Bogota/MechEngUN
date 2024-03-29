# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 12:08:07 2022

@author: Bojack Horseman
"""


import numpy as np


def cal_property(T, name, prop):
    
    if 290 <= T <= 1000:
        if name == 'CO':
            # h_est = 8669.0
            # h_form = -110530.0
            a1 = 0.03262451E+02
            a2 = 0.15119409E-02
            a3 = -0.03881755E-04
            a4 = 0.05581944E-07
            a5 = -0.02474951E-10
            a6 = -0.14310539E+05
            a7 = 0.04848897E+02
        elif name == 'CO2':
            # h_est = 9364.0
            # h_form = -393520.0
            """
            a1 = 0.24007797E+01
            a2 = 0.87350957E-02
            a3 = -0.66070878E-05
            a4 = 0.20021861E-08
            a5 = 0.63274039E-15
            a6 = -0.48377527E+05
            a7 = 0.96951457E+01
            """
            a1 = 0.02275724E+02
            a2 = 0.09922072E-01
            a3 = -0.10409113E-04
            a4 = 0.06866686E-07
            a5 = -0.02117280E-10
            a6 = -0.04837314E+06
            a7 = 0.10188488E+02
        elif name == 'H2':
            # h_est = 8468.0
            # h_form = 0.0
            a1 = 0.03298124E+02
            a2 = 0.08249441E-02
            a3 = -0.08143015E-05
            a4 = -0.09475434E-09
            a5 = 0.04134872E-11
            a6 = -0.10125209E+04
            a7 = -0.03294094E+02
        elif name == 'H':
            # h_form = 218000.0
            a1 = 0.02500000E+02
            a2 = 0.00000000E+00
            a3 = 0.00000000E+00
            a4 = 0.00000000E+00
            a5 = 0.00000000E+00
            a6 = 0.02547162E+06
            a7 = -0.04601176E+01
        elif name == 'OH':
            # h_form = 39460.0
            # h_est = 9188.0
            a1 = 0.03637266E+02
            a2 = 0.01850910E-02
            a3 = -0.16761646E-05
            a4 = 0.02387202E-07
            a5 = -0.08431442E-11
            a6 = 0.03606781E+05
            a7 = 0.13588605E+01
        elif name == 'H2O':
            # h_est = 9904.0
            # h_form = -241820.0
            a1 = 0.03386842E+02
            a2 = 0.03474982E-01
            a3 = -0.06354696E-04
            a4 = 0.06968581E-07
            a5 = -0.02506588E-10
            a6 = -0.03020811E+06
            a7 = 0.02590232E+02
        elif name == 'N2':
            # h_form = 0.0
            # h_est = 8669.0
            a1 = 0.03298677E+02
            a2 = 0.14082404E-02
            a3 = -0.03963222E-04
            a4 = 0.05641515E-07
            a5 = -0.02444854E-10
            a6 = -0.10208999E+04
            a7 = 0.03950372E+02
        elif name == 'N':
            # h_form = 472680.0
            # h_est = 8669.0
            a1 = 0.02503071E+02
            a2 = -0.02180018E-03
            a3 = 0.05420529E-06
            a4 = -0.05647560E-09
            a5 = 0.02099904E-12
            a6 = 0.05609890E+06
            a7 = 0.04167566E+02
        elif name == 'NO':
            # h_form = 88850.0  # Verificar
            a1 = 0.03376541E+02
            a2 = 0.12530634E-02
            a3 = -0.03302750E-04
            a4 = 0.05217810E-07
            a5 = -0.02446262E-10
            a6 = 0.09817961E+05
            a7 = 0.05829590E+02
        elif name == 'NO2':
            # h_form = 55565.0  # verificar
            a1 = 0.02670600E+02
            a2 = 0.07838500E-01
            a3 = -0.08063864E-04
            a4 = 0.06161714E-07
            a5 = -0.02320150E-10
            a6 = 0.02896290E+05
            a7 = 0.11612071E+02
        elif name == 'O2':
            # h_form = 0.0
            # h_est = 8682.0
            a1 = 0.03212936E+02
            a2 = 0.11274864E-02
            a3 = -0.05756150E-05
            a4 = 0.13138773E-08
            a5 = -0.08768554E-11
            a6 = -0.10052490E+04
            a7 = 0.06034737E+02
        elif name == 'O':
            # h_form = 429170.0
            # h_est = 6852.0
            a1 = 0.02946428E+02
            a2 = -0.16381665E-02
            a3 = 0.02421031E-04
            a4 = -0.16028431E-08
            a5 = 0.03890696E-11
            a6 = 0.02914764E+06
            a7 = 0.02963995E+02
        elif name == 'C3H8':
            # h_est = 0
            # h_form = 0.0
            a1 = 0
            a2 = 0
            a3 = 0
            a4 = 0
            a5 = 0
            a6 = 0
            a7 = 0
            
    if 1000 < T <= 5000:
        if name == 'CO':
            # h_est = 8669.0
            # h_form = -110530.0
            a1 = 0.03025078E+02
            a2 = 0.14426885E-02
            a3 = -0.05630827E-05
            a4 = 0.10185813E-09
            a5 = -0.06910951E-13
            a6 = -0.14268350E+05
            a7 = 0.06108217E+02
        elif name == 'CO2':
            # h_est = 9364.0
            # h_form = -393520.0
            a1 = 0.04453623E+02
            a2 = 0.03140168E-01
            a3 = -0.12784105E-05
            a4 = 0.02393996E-08
            a5 = -0.16690333E-13
            a6 = -0.04896696E+06
            a7 = -0.09553959E+01
        elif name == 'H2':
            # h_est = 8468.0
            # h_form = 0.0
            a1 = 0.02991423E+02
            a2 = 0.07000644E-02
            a3 = -0.05633828E-06
            a4 = -0.09231578E-10
            a5 = 0.15827519E-14
            a6 = -0.08350340E+04
            a7 = -0.13551101E+01
        elif name == 'H':
            # h_form = 218000.0
            a1 = 0.02500000E+02
            a2 = 0.00000000E+00
            a3 = 0.00000000E+00
            a4 = 0.00000000E+00
            a5 = 0.00000000E+00
            a6 = 0.02547162E+06
            a7 = -0.04601176E+01
        elif name == 'OH':
            # h_est = 9188.0
            # h_form = 39460.0
            a1 = 0.02882730E+02
            a2 = 0.10139743E-02
            a3 = -0.02276877E-05
            a4 = 0.02174683E-09
            a5 = -0.05126305E-14
            a6 = 0.03886888E+05
            a7 = 0.05595712E+02
        elif name == 'H2O':
            # h_est = 9904.0
            # h_form = -241820.0
            a1 = 0.02672145E+02
            a2 = 0.03056293E-01
            a3 = -0.08730260E-05
            a4 = 0.12009964E-09
            a5 = -0.06391618E-13
            a6 = -0.02989921E+06
            a7 = 0.06862817E+02
        elif name == 'N2':
            # h_est = 8669.0
            # h_form = 0.0
            a1 = 0.02926640E+02
            a2 = 0.14879768E-02
            a3 = -0.05684760E-05
            a4 = 0.10097038E-09
            a5 = -0.06753351E-13
            a6 = -0.09227977E+04
            a7 = 0.05980528E+02
        elif name == 'N':
            # h_est = 8669.0
            # h_form = 472680.0
            a1 = 0.02450268E+02
            a2 = 0.10661458E-03
            a3 = -0.07465337E-06
            a4 = 0.01879652E-09
            a5 = -0.10259839E-14
            a6 = 0.05611604E+06
            a7 = 0.04448758E+02
        elif name == 'NO':
            # h_form = 88850.0  # Verificar
            a1 = 0.03245435E+02
            a2 = 0.12691383E-02
            a3 = -0.05015890E-05
            a4 = 0.09169283E-09
            a5 = -0.06275419E-13
            a6 = 0.09800840E+05
            a7 = 0.06417293E+02
        elif name == 'NO2':
            # h_form = 55565.0  # verificar
            a1 = 0.04682859E+02
            a2 = 0.02462429E-01
            a3 = -0.10422585E-05
            a4 = 0.01976902E-08
            a5 = -0.13917168E-13
            a6 = 0.02261292E+05
            a7 = 0.09885985E+01
        elif name == 'O2':
            a1 = 0.03697578E+02
            a2 = 0.06135197E-02
            a3 = -0.12588420E-06
            a4 = 0.01775281E-09
            a5 = -0.11364354E-14
            a6 = -0.12339301E+04
            a7 = 0.03189165E+02
        elif name == 'O':
            a1 = 0.02542059E+02
            a2 = -0.02755061E-03
            a3 = -0.03102803E-07
            a4 = 0.04551067E-10
            a5 = -0.04368051E-14
            a6 = 0.02923080E+06
            a7 = 0.04920308E+02
        elif name == 'C3H8':
            a1 = 0
            a2 = 0
            a3 = 0
            a4 = 0
            a5 = 0
            a6 = 0
            a7 = 0
    
    if 298 <= T <= 2000:
        
        if name == 'HO2':
            
            A = 26.00960
            B = 34.85810	
            C = -16.30060
            D = 3.110441
            E = -0.018611
            F = -7.140991
            G = 250.7660
            H = 2.092001
            
        elif name == 'H2O2':
            if 298 <= T <= 1500:
                A = 34.25667
                B = 55.18445
                C = -35.15443
                D = 9.087440
                E = -0.422157
                F = -149.9098
                G = 257.0604
                H = -136.1064
            
            else: 
                A = 0
                B = 0
                C = 0
                D = 0
                E = 0
                F = 0
                G = 0
                H = 0        
        
    if 2000 <= T <= 6000:
        
        if name == 'HO2':
            
            A = 45.87510
            B = 8.814350
            C = -1.636031
            D = 0.098053
            E = -10.17380
            F = -26.90210
            G = 266.5260
            H = 2.092001
            
        elif name == 'H2O2':
            
            A = 0
            B = 0
            C = 0
            D = 0
            E = 0
            F = 0
            G = 0
            H = 0
                
    Ru = 8.31446261815324  # kJ/kmol K
    if prop == 'cp':
        
        if name == 'H2O2' or name == 'HO2':
            
            t = T/1000
            cp = A  + B * t  + C * (t**2) + D * (t**3) + E / (t**2)
            
        else: 
            cp = Ru * (a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4)
        return cp
    
    elif prop == 'h':
        
        if name == 'H2O2' or name == 'HO2':
            
            t = T/1000
            h = (A * t + B * (t**2) / 2 + C * (t**3) / 3 + D * (t**4) / 4 - E / t + F - H) * 1000
            
        else: 
            h = Ru * T * (a1 + ((a2/2) * T) + ((a3/3) * (T**2)) + (a4/4) * (T**3) + (a5/5) * (T**4) + (a6/T))
            
        return h
    
    elif prop == 's':
        s = Ru * (a1*np.log(T) + a2*T + a3/2 * T**2 + a4/3 * T**3 + a5/4 * T**4 + a7)
        return s

