class Element(object):

    def __init__(self, name, T, n):
        self.h_form = None
        self.h_est = None
        self.a1 = None
        self.a2 = None
        self.a3 = None
        self.a4 = None
        self.a5 = None
        self.a6 = None
        self.a7 = None
        self.name = name
        self.temperature = T
        self.Ru = 8.31446261815324  # kJ/kmol K

    def get_coefficients(self):
        if 300 <= self.temperature <= 1000:
            if self.name == 'CO':
                h_est = 8669
                h_form = -110530
                a1 = 0.03262451E+02
                a2 = 0.15119409E-02
                a3 = -0.03881755E-04
                a4 = 0.05581944E-07
                a5 = -0.02474951E-10
                a6 = -0.14310539E+05
                a7 = 0.04848897E+02
            elif self.name == 'CO2':
                h_est = 9364
                h_form = -393520
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
            elif self.name == 'H2':
                h_est = 8468
                h_form = 0
                a1 = 0.03298124E+02
                a2 = 0.08249441E-02
                a3 = -0.08143015E-05
                a4 = -0.09475434E-09
                a5 = 0.04134872E-11
                a6 = -0.10125209E+04
                a7 = -0.03294094E+02
            elif self.name == 'H':
                h_form = 218000
                a1 = 0.02500000E+02
                a2 = 0.00000000E+00
                a3 = 0.00000000E+00
                a4 = 0.00000000E+00
                a5 = 0.00000000E+00
                a6 = 0.02547162E+06
                a7 = -0.04601176E+01
            elif self.name == 'OH':
                h_form = 39460
                h_est = 9188
                a1 = 0.03637266E+02
                a2 = 0.01850910E-02
                a3 = -0.16761646E-05
                a4 = 0.02387202E-07
                a5 = -0.08431442E-11
                a6 = 0.03606781E+05
                a7 = 0.13588605E+01
            elif self.name == 'H2O':
                h_est = 9904
                h_form = -241820
                a1 = 0.03386842E+02
                a2 = 0.03474982E-01
                a3 = -0.06354696E-04
                a4 = 0.06968581E-07
                a5 = -0.02506588E-10
                a6 = -0.03020811E+06
                a7 = 0.02590232E+02
            elif self.name == 'N2':
                h_form = 0
                h_est = 8669
                a1 = 0.03298677E+02
                a2 = 0.14082404E-02
                a3 = -0.03963222E-04
                a4 = 0.05641515E-07
                a5 = -0.02444854E-10
                a6 = -0.10208999E+04
                a7 = 0.03950372E+02
            elif self.name == 'N':
                h_form = 472680
                h_est = 8669
                a1 = 0.02503071E+02
                a2 = -0.02180018E-03
                a3 = 0.05420529E-06
                a4 = -0.05647560E-09
                a5 = 0.02099904E-12
                a6 = 0.05609890E+06
                a7 = 0.04167566E+02
            elif self.name == 'NO':
                h_form = 88850 # Verificar
                a1 = 0.03376541E+02
                a2 = 0.12530634E-02
                a3 = -0.03302750E-04
                a4 = 0.05217810E-07
                a5 = -0.02446262E-10
                a6 = 0.09817961E+05
                a7 = 0.05829590E+02
            elif self.name == 'NO2':
                h_form = 55565 # verificar
                a1 = 0.02670600E+02
                a2 = 0.07838500E-01
                a3 = -0.08063864E-04
                a4 = 0.06161714E-07
                a5 = -0.02320150E-10
                a6 = 0.02896290E+05
                a7 = 0.11612071E+02
            elif self.name == 'O2':
                h_form = 0
                h_est = 8682
                a1 = 0.03212936E+02
                a2 = 0.11274864E-02
                a3 = -0.05756150E-05
                a4 = 0.13138773E-08
                a5 = -0.08768554E-11
                a6 = -0.10052490E+04
                a7 = 0.06034737E+02
            elif self.name == 'O':
                h_form = 429170
                h_est = 6852
                a1 = 0.02946428E+02
                a2 = -0.16381665E-02
                a3 = 0.02421031E-04
                a4 = -0.16028431E-08
                a5 = 0.03890696E-11
                a6 = 0.02914764E+06
                a7 = 0.02963995E+02

        if 1000 < self.temperature <= 5000:
            if self.name == 'CO':
                h_est = 8669
                h_form = -110530
                a1 = 0.03025078E+02
                a2 = 0.14426885E-02
                a3 = -0.05630827E-05
                a4 = 0.10185813E-09
                a5 = -0.06910951E-13
                a6 = -0.14268350E+05
                a7 = 0.06108217E+02
            elif self.name == 'CO2':
                h_est = 9364
                h_form = -393520
                a1 = 0.04453623E+02
                a2 = 0.03140168E-01
                a3 = -0.12784105E-05
                a4 = 0.02393996E-08
                a5 = -0.16690333E-13
                a6 = -0.04896696E+06
                a7 = -0.09553959E+01
            elif self.name == 'H2':
                h_est = 8468
                h_form = 0
                a1 = 0.02991423E+02
                a2 = 0.07000644E-02
                a3 = -0.05633828E-06
                a4 = -0.09231578E-10
                a5 = 0.15827519E-14
                a6 = -0.08350340E+04
                a7 = -0.13551101E+01
            elif self.name == 'H':
                h_form = 218000
                a1 = 0.02500000E+02
                a2 = 0.00000000E+00
                a3 = 0.00000000E+00
                a4 = 0.00000000E+00
                a5 = 0.00000000E+00
                a6 = 0.02547162E+06
                a7 = -0.04601176E+01
            elif self.name == 'OH':
                h_est = 9188
                h_form = 39460
                a1 = 0.02882730E+02
                a2 = 0.10139743E-02
                a3 = -0.02276877E-05
                a4 = 0.02174683E-09
                a5 = -0.05126305E-14
                a6 = 0.03886888E+05
                a7 = 0.05595712E+02
            elif self.name == 'H2O':
                h_est = 9904
                h_form = -241820
                a1 = 0.02672145E+02
                a2 = 0.03056293E-01
                a3 = -0.08730260E-05
                a4 = 0.12009964E-09
                a5 = -0.06391618E-13
                a6 = -0.02989921E+06
                a7 = 0.06862817E+02
            elif self.name == 'N2':
                h_est = 8669
                h_form = 0
                a1 = 0.02926640E+02
                a2 = 0.14879768E-02
                a3 = -0.05684760E-05
                a4 = 0.10097038E-09
                a5 = -0.06753351E-13
                a6 = -0.09227977E+04
                a7 = 0.05980528E+02
            elif self.name == 'N':
                h_est = 8669
                h_form = 472680
                a1 = 0.02450268E+02
                a2 = 0.10661458E-03
                a3 = -0.07465337E-06
                a4 = 0.01879652E-09
                a5 = -0.10259839E-14
                a6 = 0.05611604E+06
                a7 = 0.04448758E+02
            elif self.name == 'NO':
                h_form = 88850  # Verificar
                a1 = 0.03245435E+02
                a2 = 0.12691383E-02
                a3 = -0.05015890E-05
                a4 = 0.09169283E-09
                a5 = -0.06275419E-13
                a6 = 0.09800840E+05
                a7 = 0.06417293E+02
            elif self.name == 'NO2':
                h_form = 55565  # verificar
                a1 = 0.04682859E+02
                a2 = 0.02462429E-01
                a3 = -0.10422585E-05
                a4 = 0.01976902E-08
                a5 = -0.13917168E-13
                a6 = 0.02261292E+05
                a7 = 0.09885985E+01
            elif self.name == 'O2':
                h_form = 0
                h_est = 8682
                a1 = 0.03697578E+02
                a2 = 0.06135197E-02
                a3 = -0.12588420E-06
                a4 = 0.01775281E-09
                a5 = -0.11364354E-14
                a6 = -0.12339301E+04
                a7 = 0.03189165E+02
            elif self.name == 'O':
                h_form = 429170
                h_est = 6852
                a1 = 0.02542059E+02
                a2 = -0.02755061E-03
                a3 = -0.03102803E-07
                a4 = 0.04551067E-10
                a5 = -0.04368051E-14
                a6 = 0.02923080E+06
                a7 = 0.04920308E+02
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        self.a5 = a5
        self.a6 = a6
        self.a7 = a7
        self.h_est = h_est
        self.h_form = h_form

    def get_enthalpy(self):
        h = self.Ru * self.temperature * (self.a1 + ((self.a2 / 2) * self.temperature) +
                                          ((self.a3 / 3) * (self.temperature ** 2)) +
                                          (self.a4 / 4) * (self.temperature ** 3) +
                                          (self.a5 / 5) * (self.temperature ** 4) + (self.a6 / self.temperature))
        return h

    def get_entropy(self):
        s = self.Ru * (self.a1 * np.log(self.temperature) + self.a2 * self.temperature +
                       self.a3 / 2 * self.temperature ** 2 + self.a4 / 3 * self.temperature ** 3 +
                       self.a5 / 4 * self.temperature ** 4 + self.a7)
        return s

    def get_cp(self):
        cp = self.Ru * (self.a1 + self.a2 * self.temperature + self.a3 * self.temperature ** 2 +
                        self.a4 * self.temperature ** 3 + self.a5 * self.temperature ** 4)
        return cp



        
