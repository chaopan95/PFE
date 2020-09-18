#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:31:13 2020

@author: panchao
"""

import os
from numpy import zeros
from math import sqrt, exp, log
import matplotlib.pyplot as plt


def iniparser_load(input_file_name):
    """
    Parse input file
    """
    ini = {'model': 'Seyeux_2010', 'save_history': 'True',
           'output_file_name': 'c_solver_corrosion_evolution_test_1.output',
           'temperature_in_K': 603, 'pH_temperature': 7.2,
           'time_in_seconds': 36000, 'x_Cr': 0.3, 'x_Fe': 0.3, 'x_Ni': 0.3,
           'T': 603, 'pH': 7.2, 'alpha': 0.5, 'DV': 0.5, 'F_0f': 10,
           'F_0mf': 0.2, 'F_0fs': 0.2, 'CtotM_mf': 4.14e28, 'CtotI': 2.07e28,
           'CmFe_mf': 40.6e14, 'CMNi_mf': 40.6e14, 'CVCr_mf': 40.6e14,
           'CtotM_fs': 4.14e28, 'CmFe_fs': 40.6e7, 'CMNi_fs': 40.6e7,
           'CVCr_fs': 40.6e20, 'D_vO': 1e-19, 'D_mCr': 0, 'D_ICr': 0,
           'D_lCr': 1e-21, 'DG1': 100000, 'DG2': -100000, 'DG3': -100000,
           'DG4': -100000, 'DG5': -100000, 'DG6': 100000, 'DG7': -100000,
           'DG8': -100000, 'DG9': -100000, 'DG10': -100000, 'DG11': 100000,
           'DG12': -100000, 'DG13': -100000, 'DG14': -100000,
           'decay_length': 1, 'charge_number': 3, 'dissol_order': 1,
           'dissol_preexp': 0, 'dissol_Ea': 100000, 'q_vO': 2, 'q_lCr': 3,
           'q_mCr': 3, 'gamma_H': 1}

    with open(input_file_name) as f:
        contents = [line.strip() for line in f.readlines()]
    #header = contents[0][1:-1]
    for content in contents[1:]:
        if ' = ' not in content:
            continue
        assignment = content.split(' = ')
        try:
            value = float(assignment[1])
        except ValueError:
            value = assignment[1]
        ini[assignment[0]] = value

    assert ini['temperature_in_K'] > 0,\
        'Parameter temperature_in_K must be > 0'
    assert ini['pH_temperature'] > 0, 'Parameter pH_temperature must be > 0'
    assert ini['x_Cr'] > 0, 'Parameter x_Cr must be > 0'
    assert ini['time_in_seconds'] > 0, 'Parameter time_in_seconds must be > 0'
    assert 0 <= ini['gamma_H'] <= 1, 'Parameter gamma_H must be [0, 1]'
    assert ini['model'] in ['Seyeux_2010', 'Leistner_2012', 'Voyshnis_2014'],\
        'model does not belong to [Seyeux_2010, Leistner_2012, Voyshnis_2014]'
    return ini



class Oxydation:
    def __init__(self,
                 model='Voyshnis_2014',
                 save_history='True',
                 output_file_name='c_solver_corrosion_evolution_test_8.output',
                 temperature_in_K=600,
                 pH_temperature=7.2,
                 time_in_seconds=7776000,
                 # Cr, Fe, Ni mole fraction
                 x_Cr=0.32, x_Fe=0.10, x_Ni=0.58,
                 # Temperature in degree K
                 T=600,
                 # water chemistry
                 pH=7.2, 
                 # interface polarisability between 0 and 1
                 alpha=0.5,
                 # potential change with respect to a stationnary state (V)
                 DV=0.01,
                 # potential drop in the film, at metal-film interface, 
                 # at film-solution interface (V)
                 F_0f=0.1, F_0mf=0.001, F_0fs=0.3,
                 # total cation concentration in oxide at the metal-film
                 # interface
                 CtotM_mf=4.14e28,
                 # Interstice concentration in Cr2O3 oxide
                 CtotI=2.07e28,
                 # Fe3+ concentration in oxide at the metal-film interface
                 CmFe_mf=40.6e14,
                 # Ni2+ concentration in oxide at the metal-film interface
                 CMNi_mf=40.6e14,
                 # Cr3+ vacancy concentration in oxide at the metal-film
                 # interface
                 CVCr_mf=40.6e14,
                 # total cation concentration at oxide-solution interface
                 CtotM_fs=4.14e28,
                 # Fe3+ concentration at oxide-solution interface
                 CmFe_fs=40.6e7,
                 # Ni2+ concentration at oxide-solution interface
                 CMNi_fs=40.6e7,
                 # Cr3+ vacancy concentration at oxide-solution interface
                 CVCr_fs=40.6e20,
                 # diffusion coefficient of oxygen vacancies, chromium ions,
                 # chromium ions in interstices, Cr3+ cation
                 D_vO=5.3e-23, D_mCr=7.0e-24, D_ICr=1.0e-20, D_lCr=1e-21,
                 # Gibbs enrgy of formation of reaction 1-14
                 DG1=-15000, DG2=-90000, DG3=85000, DG4=-30000,
                 DG5=10000, DG6=-85000, DG7=10000, DG8=-100000,
                 DG9=100000, DG10=-1000, DG11=30000, DG12=-6000,
                 DG13=10000, DG14=-3000,
                 decay_length=5e-8,
                 # number of electrons transferred during dissolution reaction
                 charge_number=3,
                 # dissolution
                 dissol_order=1, dissol_preexp=0.84e17, dissol_Ea=5.4e4,
                 # electronic charge of oxygen vacancy, Cr3+ cation, Cr reseau
                 q_vO=2, q_lCr=3, q_mCr=3,
                 # parameters of Voyshnis 2014 for H
                 gamma_H=0.2,
                 # oxide molar volume in m³/mol
                 Omega=2.92e-3,
                 # Avogadro number (must be equal to 6.022e23 ?)
                 Nv=6.022e23,
                 # Faraday constant
                 F=9.65e4,
                 # Perfect gas constant in J/mol/K
                 R=8.314472,
                 # electronic charge in C
                 e_charge=1.602e-19,
                 # Boltzmann's constant in J/K
                 kB=1.38e-23,
                 # oxide density in g/cm³
                 rho_oxide=5.22,
                 # oxide molar mass in g/mol
                 M_oxide=151.9904,
                 # Cr stoichiometry in oxide
                 stoich_Cr=2,
                 ):
        """
        # Gibbs enrgy of formation of reaction 1-14
        # 1:  (Cr_M -> Cr3+_ox + 3e- + 3/2 V_o with V_o = oxygen vacancy)
        # 2:  (Cr_M + V_Cr -> Cr3+_ox + 3e)
        # 3:  (Cr_M + V_I -> I_Cr + 3e)
        # 4:  (Fe_M + V_Cr -> Fe3+_ox + 3e)
        # 5:  (Fe_M + V_I -> I_Fe + 3e)
        # 6:  (Ni_M + 2/3 V_Cr -> Ni2+_ox + 2e-)
        # 7:  (Ni_M + V_I -> I_Ni + 2e)
        # 8:  (V_o + H2O -> 2 H+ + O_ox)
        # 9:  (Cr3+_ox + 3/2 H2O -> 3/2 O_ox + Cr3+_ox + V_Cr)
        # 10: (I_Cr + 3/2 H2O -> 3/2 O_ox + Cr3+_ox + 3H+ + V_I)
        # 11: (Fe3+_ox -> V_Cr + Fe3+_aq)
        # 12: (I_Fe -> V_I + Fe3+_aq)
        # 13: (M_Ni -> 2/3 V_Cr + Ni2+_aq)
        # 14: (I_Ni -> V_I + Ni2+_aq)
        """
        self.model = model
        self.save_history = save_history
        self.output_file_name = output_file_name
        self.temperature_in_K = temperature_in_K
        self.pH_temperature = pH_temperature
        self.time_in_seconds = time_in_seconds
        self.x_Cr, self.x_Fe, self.x_Ni = x_Cr, x_Fe, x_Ni
        self.T, self.pH, self.alpha, self.DV = T, pH, alpha, DV
        self.F_0f, self.F_0mf, self.F_0fs = F_0f, F_0mf, F_0fs
        self.CtotM_mf = CtotM_mf
        self.CtotI = CtotI
        self.CmFe_mf, self.CMNi_mf = CmFe_mf, CMNi_mf
        self.CVCr_mf = CVCr_mf
        self.CtotM_fs, self.CmFe_fs = CtotM_fs, CmFe_fs
        self.CMNi_fs, self.CVCr_fs = CMNi_fs, CVCr_fs
        self.D_vO, self.D_mCr, self.D_ICr = D_vO, D_mCr, D_ICr
        self.D_lCr = D_lCr
        self.DG1, self.DG2, self.DG3, self.DG4 = DG1, DG2, DG3, DG4
        self.DG5, self.DG6, self.DG7, self.DG8 = DG5, DG6, DG7, DG8
        self.DG9, self.DG10, self.DG11 = DG9, DG10, DG11
        self.DG12, self.DG13, self.DG14 = DG12, DG13, DG14
        self.decay_length = decay_length
        self.charge_number = charge_number
        self.dissol_order, self.dissol_preexp = dissol_order, dissol_preexp
        self.dissol_Ea = dissol_Ea
        self.q_vO, self.q_lCr, self.q_mCr = q_vO, q_lCr, q_mCr
        self.gamma_H = gamma_H
        self.Omega, self.Nv, self.F, self.R = Omega, Nv, F, R
        self.e_charge, self.kB, self.rho_oxide = e_charge, kB, rho_oxide
        self.M_oxide, self.stoich_Cr = M_oxide, stoich_Cr

    def print_parameters(self):
        print('Parameters:')
        print("                             model =", self.model)
        print("                  temperature_in_K =", self.temperature_in_K)
        print("                    pH_temperature =", self.pH_temperature)
        print("                              x_Cr =", self.x_Cr)
        print("                              x_Fe =", self.x_Fe)
        print("                              x_Ni =", self.x_Ni)
        print("                   time_in_seconds =", self.time_in_seconds)
        print("                      save_history =", self.save_history)
        print("                  output_file_name =", self.output_file_name)
        print("                             alpha =", self.alpha)
        print("                                DV =", self.DV)
        print("                              F_0f =", self.F_0f)
        print("                             F_0mf =", self.F_0mf)
        print("                             F_0fs =", self.F_0fs)
        print("                              D_vO =", self.D_vO)
        print("                             D_mCr =", self.D_mCr)
        print("                             D_ICr =", self.D_ICr)
        print("                               DG1 =", self.DG1)
        print("                               DG8 =", self.DG8)
        print("                               DG2 =", self.DG2)
        print("                               DG4 =", self.DG4)
        print("                               DG6 =", self.DG6)
        print("                               DG9 =", self.DG9)
        print("                              DG11 =", self.DG11)
        print("                              DG13 =", self.DG13)
        print("                              DG10 =", self.DG10)
        print("                              DG12 =", self.DG12)
        print("                              DG14 =", self.DG14)
        print("                               DG3 =", self.DG3)
        print("                               DG5 =", self.DG5)
        print("                               DG7 =", self.DG7)
        print("                          CtotM_mf =", self.CtotM_mf)
        print("                             CtotI =", self.CtotI)
        print("                      decay_length =", self.decay_length)
        print("                     charge number =", self.charge_number)
        print("                order, dissolution =", self.dissol_order)
        print("preexponential factor, dissolution =", self.dissol_preexp)
        print("    activation energy, dissolution =", self.dissol_Ea)
        print("                           gamma_H =", self.gamma_H)

    def set_user_parameters(self, isPrint=False):
        """
        """
        self.T = self.temperature_in_K
        self.pH = self.pH_temperature

        ratio = (self.F_0f + self.DV*(1-self.alpha)) / self.F_0f

        print("ratio caracterising the square root shape of the curve",
              "(when close to 1) = {}".format(ratio))

        #******************* CALCULATE LUMPED PARAMETERS *******************#
        # JvO
        K1 = (self.q_vO*self.F*self.D_vO)/(self.R*self.T) *\
            exp(((-1*2./3*self.DG1) + (2.*self.F*self.F_0mf))/(self.R*self.T) +
            2./3.*log(self.x_Cr))
        P1 = K1 * (self.F_0f+self.DV)
        P2 = -1. * K1 * self.DV * self.alpha

        # JmCr
        C_VCr_mf = exp((self.DG2-3.*self.F*self.F_0mf)/(self.R*self.T) -
                       log(self.x_Cr))

        C_MFe_mf = exp((-1.*self.DG4+3.*self.F*self.F_0mf)/
                       (self.R*self.T)) * self.x_Fe * C_VCr_mf

        C_MNi_mf = exp((-1.*self.DG6+2.*self.F*self.F_0mf)/
                       (self.R*self.T)) * self.x_Ni * pow(C_VCr_mf, (2./3.))

        C_MCr_mf = self.CtotM_mf*self.Omega/self.Nv - (1.*C_MNi_mf +
                                                       1.*C_MFe_mf + C_VCr_mf)

        P3 = ((self.q_mCr*self.F*self.D_mCr)/
              (self.R*self.T)) * C_MCr_mf * (self.F_0f+self.DV)

        P4 = -1.*((self.q_mCr*self.F*self.D_mCr)/(self.R*self.T)) *\
            C_MCr_mf * self.alpha * self.DV

        # JICr
        eterm_1 = self.x_Cr*exp((-1.*self.DG3+3.*self.F*self.F_0mf)/
                                self.R/self.T)

        eterm_2 = self.x_Fe*exp((-1.*self.DG5+3.*self.F*self.F_0mf)/
                                self.R/self.T)

        eterm_3 = self.x_Ni*exp((-1.*self.DG7+2.*self.F*self.F_0mf)/
                                self.R/self.T)

        C_VI_mf = self.CtotI*self.Omega/self.Nv /\
            (1. + eterm_1 + eterm_2 + eterm_3)

        C_ICr_mf = eterm_1 * C_VI_mf

        P5 = ((self.q_mCr*self.F*self.D_ICr)/(self.R*self.T)) * C_ICr_mf *\
            (self.F_0f+self.DV)

        P6 = -1.*((self.q_mCr*self.F*self.D_ICr)/(self.R*self.T)) *\
            C_ICr_mf * self.alpha * self.DV

        P = P1 + P3 + P5
        Q = P2 + P4 + P6

        # Jdissol
        H_plus = exp(-1.*2.303*self.pH)

        Sdiss = self.charge_number*self.e_charge/self.kB/self.T *\
            self.alpha * self.DV

        Rdiss = self.dissol_preexp *\
            exp((-1.*self.dissol_Ea/self.kB/self.Nv/self.T +
                 self.charge_number*self.e_charge*self.F_0fs/self.kB/self.T))
        if isPrint:
            print("K1 = %.12e"%K1)
            print("P1 = %.12e"%P1)
            print("P2 = %.12e"%P2);
            print("P3 = %.12e"%P3);
            print("P4 = %.12e"%P4);
            print("P5 = %.12e"%P5);
            print("P6 = %.12e"%P6);
            print("P = %.12e"%P);
            print("Q = %.12e"%Q);
            print("Rdiss = %.12e"%Rdiss);
            print("Sdiss = %.12e"%Sdiss);
        return

    def recopie_param_xi(self, x, opt):
        """
        Si opt[i]= 1 x[i] est optimisé
        Si opt[i]= 0 x[i] n'est pas optimisé
        """
        # lacunes d'oxygene
        x[0], opt[0] = self.F_0fs, 1
        x[1], opt[1] = self.F_0f, 1
        x[2], opt[2] = self.F_0mf, 1
        x[3], opt[3] = self.alpha, 1

        x[4], opt[4] = self.DG1, 1
        x[5], opt[5] = self.D_vO, 1
        x[6], opt[6] = self.DG8, 1

        x[7], opt[7] = self.q_vO, 0

        x[8], opt[8] = self.CmFe_mf, 0
        x[9], opt[9] = self.CMNi_mf, 0
        x[10], opt[10] = self.CVCr_mf, 0
        x[11], opt[11] = self.CtotM_fs, 0
        x[12], opt[12] = self.CmFe_fs, 0
        x[13], opt[13] = self.CMNi_fs, 0
        x[14], opt[14] = self.CVCr_fs, 0
        # Cr intersticiels
        x[15], opt[15] = self.q_lCr, 0
        x[16], opt[16] = self.D_lCr, 0
        x[21], opt[17] = self.CtotM_mf, 0

        # Cr reseau
        x[22], opt[18] = self.q_mCr, 0
        x[23], opt[19] = self.D_mCr, 0

        return x, opt

    def recopie_xi_param(self, x):
        """
        """
        self.F_0fs = x[0]
        self.F_0f = x[1]
        self.F_0mf = x[2]

        # lacunes d'oxygene
        self.alpha = x[3]
        self.DG1 = x[4]
        self.D_vO = x[5]
        self.DG8 = x[6]
        self.q_vO = x[7]

        self.CmFe_mf = x[8]
        self.CMNi_mf = x[9]
        self.CVCr_mf = x[10]
        self.CtotM_fs = x[11]
        self.CmFe_fs = x[12]
        self.CMNi_fs = x[13]
        self.CVCr_fs = x[14]

        # Cr intersticiels
        self.q_lCr = x[15]
        self.D_lCr = x[16]
        self.CtotM_mf = x[17]

        # Cr reseau
        self.q_mCr = x[18]
        self.decay_length = x[19]
        self.DV = x[20]
        return

    def sauver(self, nom, n, t, y):
        """
        """
        with open(nom, 'w') as f:
            f.write('t[s], X_calc[m], X_exp[m]\n')
            for i in range(n+1):
                if i%3600 == 0:
                    f.write('%.12e, %.12e, %.12e\n'%(t[i], y[i],
                                                     self.fit_exp(t[i])))
            f.write('\n')

    def f_(self, x):
        """
        A decreqsing function reflecting the variation of the potential drop
        at the metal/oxide interface during the oxide growth
        """
        ratio = x/self.decay_length;
        return exp(-ratio)

    def fit_exp(self, t):
        """
        """
        # conversion en heures
        tt = t/3600.0
        val = 13.6*log(0.079*tt+1);

        # conversion en Angstrom
        return val*1e-9

    def time_integration_RK2_Seyeux_2010(self, h, n, t, y):
        """
        To confirm:
        in C code, for(i=0;i<=n; i++)
        """
        for i in range(n):
            k1 = h * self.f_Seyeux_2010(t[i], y[i]);
            # first estimation of intermediate time step value
            yi_plus_one_half = y[i] + k1/2;
            y[i+1] = y[i] + h*(self.f_Seyeux_2010(t[i]+h/2, yi_plus_one_half))
            t[i+1] = t[i] + h
            if yi_plus_one_half > y[i+1]:
                y[i+1] = y[i] + h*(self.f_Seyeux_2010(t[i]+h/2, y[i]))
                # second estimation of intermediate time step value
                yi_plus_one_half = 0.5 * (y[i] + y[i+1])
                y[i+1] = y[i] + h*self.f_Seyeux_2010(t[i]+h/2,
                                                     yi_plus_one_half)
        print(y[0], y[1])
        return t, y

    def f_Seyeux_2010(self, t, x):
        """
        """
        val = (self.Omega/self.Nv) * (self.J_vO_Seyeux_2010(x))
        return val

    def J_vO_Seyeux_2010(self, x):
        """
        oxygen vacancy flux
        """
        # new expression reviewed with Antoine Seyeux
        term_1 = (self.q_vO * self.F) * (self.F_0f + self.DV *
                 (1 - self.alpha * self.f_(x))) / (self.R * self.T * x)

        term_2 = exp((-2/3 * self.DG1 + 2 * self.F * self.F_0mf)/
                     (self.R * self.T) + 1.5 * log(self.x_Cr))

        term_3 = exp(self.q_vO * self.F *
                     (self.F_0f + self.DV * (1 - self.alpha * self.f_(x))) /
                     (self.R * self.T))

        term_4 = -exp(((self.DG8 - 2 * self.F *
                        (self.F_0fs + self.DV * self.alpha * self.f_(x)))/
                      (self.R * self.T)) - 4.606 * self.pH)

        term_5 = exp(self.q_vO * self.F *
                     (self.F_0f + self.DV * (1 - self.alpha * self.f_(x))) /
                     (self.R * self.T)) - 1

        val = (self.Nv / self.Omega) * term_1 * self.D_vO * (term_2 * term_3 -
                                                             term_4) / term_5

        return val

    def time_integration_RK2_Leistner_2012(self, h, n, t, y, J_v0):
        """
        # for(i=0;i<=n; i++)
        """
        J_v0[0] = self.J_vO_Leistner_2012(y[0])
        for i in range(n):
            k1 = h*self.f_Leistner_2012(y[i], y[i])
            # first estimation of intermediate time step value
            yi_plus_one_half = y[i] + k1/2
            y[i+1] = y[i] + h*(self.f_Leistner_2012(y[i], yi_plus_one_half))
            t[i+1] = t[i] + h
            if (yi_plus_one_half > y[i+1]):
                y[i+1] = y[i] + h*(self.f_Leistner_2012(y[i], y[i]))
                # second estimation of intermediate time step value
                yi_plus_one_half = 0.5 * (y[i] + y[i+1])
                y[i+1] = y[i] + h* self.f_Leistner_2012(y[i], yi_plus_one_half)
            J_v0[i+1] = self.J_vO_Leistner_2012(y[i+1])
        return t, y, J_v0

    def time_integration_RK2_Voyshnis_2014(self, h, n, t, y, J_v0,
                                           J_MCr, J_ICr, J_H):
        """
        """
        #print("time integration with RK2 method of Voyshnis 2014 model")
        J_v0[0] = self.J_vO_Leistner_2012(y[0])
        J_ICr[0] = self.J_ICr_Leistner_2012(y[0])
        J_MCr[0] = self.J_MCr_Leistner_2012(y[0])
        J_H[0] = 0.

        for i in range(n):
            k1 = h * self.f_Leistner_2012(y[i], y[i])
            # first estimation of intermediate time step value
            yi_plus_one_half = y[i] + k1/2
            y[i+1] = y[i] + h * (self.f_Leistner_2012(y[i], yi_plus_one_half))
            t[i+1] = t[i] + h
            if yi_plus_one_half > y[i+1]:
                y[i+1] = y[i] + h * (self.f_Leistner_2012(y[i], y[i]))
                # second estimation of intermediate time step value
                yi_plus_one_half = 0.5 * (y[i] + y[i+1])
                y[i+1] = y[i] + h * (self.f_Leistner_2012(y[i],
                                                          yi_plus_one_half))
            J_v0[i+1] = 2. * self.J_vO_Leistner_2012(y[i+1])
            J_ICr[i+1] = 3. * self.J_ICr_Leistner_2012(y[i+1])
            J_MCr[i+1] = 3. * self.J_MCr_Leistner_2012(y[i+1])
            J_H[i+1] = self.gamma_H * (J_v0[i+1] + J_ICr[i+1] + J_MCr[i+1])
        #print("t = %f, J_H = %e"%(t[i+1], J_H[i+1]))
        return t, y, J_v0, J_MCr, J_ICr, J_H

    def f_Leistner_2012(self, x_0, x):
        """
        """
        # test on x_0, initial value of thickness because x could be to
        # big because of RK2 integration scheme
        if x_0 > 1e-10:
            j_diss = self.J_dissol_Leistner_2012(x);
        else:
            j_diss = 0.;

        val = (self.Omega/self.Nv)*(self.J_vO_Leistner_2012(x)+
               self.J_MCr_Leistner_2012(x)+self.J_ICr_Leistner_2012(x)-j_diss)

        # /* comparison with approximated solution with lumped parameters */
        f_approx = self.approx_sol(x_0,x)

        if val < 0:
            print("val < 0. => x_0=%e, x=%e, JvO=%e, JmCr=%e, Jdiss=%e"%(x_0,
                  x, self.J_vO_Leistner_2012(x), self.J_MCr_Leistner_2012(x),
                  j_diss))
        return val

    def approx_sol(self, x_0, x):
        """
        """
        # JvO
        K1 = (self.q_vO*self.F*self.D_vO)/(self.R*self.T) *\
        exp(((-1.*2./3.*self.DG1) + (2.*self.F*self.F_0mf))/(self.R*self.T) +
            2./3.*log(self.x_Cr))

        P1 = K1 * (self.F_0f + self.DV)
        P2 = -1. * K1 * self.DV * self.alpha

        # JmCr
        C_VCr_mf = exp((self.DG2-3.*self.F*self.F_0mf)/(self.R*self.T) -
                       log(self.x_Cr))

        C_MFe_mf = exp((-1.*self.DG4+3.*self.F*self.F_0mf)/(self.R*self.T)) *\
        self.x_Fe * C_VCr_mf

        C_MNi_mf = exp((-1.*self.DG6+2.*self.F*self.F_0mf)/(self.R*self.T)) *\
        self.x_Ni * pow(C_VCr_mf,(2./3.))

        C_MCr_mf = self.CtotM_mf*self.Omega/self.Nv -\
        (1.*C_MNi_mf + 1.*C_MFe_mf + C_VCr_mf)

        P3 = ((self.q_mCr*self.F*self.D_mCr)/(self.R*self.T)) * C_MCr_mf *\
        (self.F_0f+self.DV)

        P4 = -1.*((self.q_mCr*self.F*self.D_mCr)/(self.R*self.T)) *\
        C_MCr_mf * self.alpha * self.DV

        # JICr
        eterm_1 = self.x_Cr*exp((-1.*self.DG3+3.*self.F*self.F_0mf)/
                                self.R/self.T)

        eterm_2 = self.x_Fe*exp((-1.*self.DG5+3.*self.F*self.F_0mf)/
                                self.R/self.T)

        eterm_3 = self.x_Ni*exp((-1.*self.DG7+2.*self.F*self.F_0mf)/
                                self.R/self.T)

        C_VI_mf = self.CtotI*self.Omega/self.Nv /\
        (1. + eterm_1 + eterm_2 + eterm_3)

        C_ICr_mf = eterm_1 * C_VI_mf

        P5 = ((self.q_mCr*self.F*self.D_ICr)/(self.R*self.T)) * C_ICr_mf *\
        (self.F_0f+self.DV)

        P6 = -1.*((self.q_mCr*self.F*self.D_ICr)/(self.R*self.T)) *\
        C_ICr_mf * self.alpha * self.DV

        P = P1 + P3 + P5
        Q = P2 + P4 + P6

        # Jdissol
        H_plus = exp(-1.*2.303*self.pH)

        Sdiss = self.charge_number*self.e_charge/self.kB/self.T *\
        self.alpha * self.DV

        Rdiss = self.dissol_preexp * H_plus *\
        exp((-1.*self.dissol_Ea/self.kB/self.Nv/self.T +
             self.charge_number*self.e_charge*self.F_0fs/self.kB/self.T))

        if x_0 > 1e-10:
            j_diss = (self.Omega/self.Nv) * Rdiss*exp(Sdiss*self.f_(x))
        else:
            j_diss = 0.

        f_approx = (P + Q*self.f_(x))/x - j_diss

        return f_approx

    def J_dissol_Leistner_2012(self, x):
        """
        film dissolution rate
        """
        H_plus = exp(-2.303*self.pH)
        F_fs = self.F_0fs + self.alpha * self.DV * exp(-1.*(x / 
                                                            self.decay_length))

        # 1/m2/s
        val = self.dissol_preexp * exp(-1.*(self.dissol_Ea/self.kB/self.Nv/
                                            self.T - self.charge_number*
                                            self.e_charge*self.Nv*F_fs/
                                            self.kB/self.Nv/self.T)) *\
              pow(H_plus, self.dissol_order)
        return val

    def J_vO_Leistner_2012(self, x):
        """
        Calculate oxygen vacancy flux
        """
        term_1 = (self.q_vO*self.F) *\
        (self.F_0f + self.DV*((1.-self.alpha*self.f_(x)))) /\
        (self.R*self.T*x)

        term_2 = exp(((-2./3.)*self.DG1+2.*self.F*self.F_0mf)/
                     (self.R*self.T) + (2./3.)*log(self.x_Cr))

        term_3 = exp(self.q_vO*self.F*
                     (self.F_0f+self.DV*((1.-self.alpha*self.f_(x)))) /
                     (self.R*self.T))

        term_4 = exp(((self.DG8-2.*self.F*
                       (self.F_0fs+self.DV*self.alpha*self.f_(x))) /
                     (self.R*self.T))-4.606*self.pH)

        term_5 = exp(self.q_vO*self.F*
                     (self.F_0f+self.DV*((1.-self.alpha*self.f_(x)))) /
                     (self.R*self.T) ) - 1.

        result = (self.Nv/self.Omega) * term_1 * self.D_vO *\
        (term_2 * term_3- term_4)/term_5

        return result

    def J_MCr_Leistner_2012(self, x):
        """
        Chrome ion flux through metallic positions
        """
        # calculate concentrations at m/f interface */
        C_VCr_mf = exp((self.DG2-3.*self.F*self.F_0mf)/
                       (self.R*self.T) - log(self.x_Cr))

        C_MFe_mf = exp((-1.*self.DG4+3.*self.F*self.F_0mf)/
                       (self.R*self.T)) * self.x_Fe * C_VCr_mf

        C_MNi_mf = exp((-1.*self.DG6+2.*self.F*self.F_0mf)/
                       (self.R*self.T)) * self.x_Ni *\
                       pow(C_VCr_mf, (2./3.))

        C_MCr_mf = self.CtotM_mf*self.Omega/self.Nv -\
        (1.*C_MNi_mf + 1.*C_MFe_mf + C_VCr_mf)

        # calculate concentrations at f/s interface */
        F_fs = self.F_0fs + self.alpha * self.DV * self.f_(x)
        C_VCr_fs = exp((-1.*self.DG9+3.*self.F*F_fs)/(self.R*self.T))

        # [mol/m^3]
        C_Ni2aq = self.rho_oxide*self.stoich_Cr/self.M_oxide/1.e-6*\
        self.x_Ni/self.x_Cr

        # [mol/m^3]
        C_Fe3aq = self.rho_oxide*self.stoich_Cr/self.M_oxide/1.e-6*\
        self.x_Fe/self.x_Cr

        C_Ni2aq = C_Ni2aq*self.Omega
        C_Fe3aq = C_Fe3aq*self.Omega

        C_MFe_fs = exp((self.DG11-3.*self.F*F_fs)/(self.R*self.T)) *\
        C_VCr_fs * C_Fe3aq

        C_MNi_fs = exp((self.DG13-2.*self.F*F_fs)/(self.R*self.T)) *\
        pow(C_VCr_fs, (2./3.)) * C_Ni2aq

        C_MCr_fs = self.CtotM_fs*self.Omega/self.Nv -\
        (1.*C_MNi_fs + 1.*C_MFe_fs + C_VCr_fs)

        # calculate flux across film */
        term_1 = (self.q_mCr*self.F) *\
        (self.F_0f + self.DV*((1.-self.alpha*self.f_(x)))) / (self.R*self.T*x)

        term_3 = exp(self.q_mCr*self.F*
                     (self.F_0f+self.DV*((1.-self.alpha*self.f_(x)))) /
                     (self.R*self.T))

        term_5 = exp(self.q_mCr*self.F*
                     (self.F_0f+self.DV*((1.-self.alpha*self.f_(x)))) /
                     (self.R*self.T)) - 1.

        val = (self.Nv/self.Omega) * term_1 * self.D_mCr *\
        (C_MCr_mf * term_3 - C_MCr_fs)/term_5

        return val

    def J_ICr_Leistner_2012(self, x):
        """
        Chrome ion flux through interstitial positions
        """
        # calculate concentrations at m/f interface */
        eterm_1 = self.x_Cr*exp((-1.*self.DG3+3.*self.F*self.F_0mf)/
                                self.R/self.T)
        eterm_2 = self.x_Fe*exp((-1.*self.DG5+3.*self.F*self.F_0mf)/
                                self.R/self.T)
        eterm_3 = self.x_Ni*exp((-1.*self.DG7+2.*self.F*self.F_0mf)/
                                self.R/self.T)

        C_VI_mf = self.CtotI*self.Omega/self.Nv /\
        (1. + eterm_1 + eterm_2 + eterm_3)

        C_ICr_mf = eterm_1 * C_VI_mf
        C_IFe_mf = eterm_2 * C_VI_mf
        C_INi_mf = eterm_3 * C_VI_mf

        # calculate concentrations at f/s interface */
        F_fs = self.F_0fs + self.alpha * self.DV * self.f_(x)

        # [mol/m^3]
        C_Ni2aq = self.rho_oxide*self.stoich_Cr/self.M_oxide/1e-6*\
        self.x_Ni/self.x_Cr

        # [mol/m^3]
        C_Fe3aq = self.rho_oxide*self.stoich_Cr/self.M_oxide/1e-6*\
        self.x_Fe/self.x_Cr

        C_Ni2aq = C_Ni2aq*self.Omega
        C_Fe3aq = C_Fe3aq*self.Omega

        eterm_4 = exp((self.DG10-3.*self.F*F_fs)/self.R/self.T - 6.909*self.pH)
        eterm_5 = C_Fe3aq*exp((self.DG12-3.*self.F*F_fs)/self.R/self.T)
        eterm_6 = C_Ni2aq*exp((self.DG14-2.*self.F*F_fs)/self.R/self.T)

        C_VI_fs = self.CtotI*self.Omega/self.Nv /\
        (1. + eterm_4 + eterm_5 + eterm_6)

        C_ICr_fs = eterm_4 * C_VI_fs
        C_IFe_fs = eterm_5 * C_VI_fs
        C_INi_fs = eterm_6 * C_VI_fs

        # calculate flux across film */
        term_1 = (self.q_mCr*self.F) *\
        (self.F_0f + self.DV*((1.-self.alpha*self.f_(x)))) / (self.R*self.T*x)

        term_3 = exp(self.q_mCr*self.F*
                     (self.F_0f+self.DV*((1.-self.alpha*self.f_(x)))) /
                     (self.R*self.T))

        term_5 = exp(self.q_mCr*self.F*
                     (self.F_0f+self.DV*((1.-self.alpha*self.f_(x)))) /
                     (self.R*self.T)) - 1.

        val = (self.Nv/self.Omega) * term_1 * self.D_ICr *\
        (C_ICr_mf * term_3 - C_ICr_fs)/term_5

        return val

    def load_params(self, param):
        """
        """
        param[0] = self.F_0fs
        param[1] = self.F_0f
        param[2] = self.F_0mf
        param[3] = self.alpha
        param[4] = self.DG1
        param[5] = self.D_vO
        param[6] = self.DG8
        param[7] = self.q_vO
        param[8] = self.CmFe_mf
        param[9] = self.CMNi_mf
        param[10] = self.CVCr_mf
        param[11] = self.CtotM_fs
        param[12] = self.CmFe_fs
        param[13] = self.CMNi_fs
        param[14] = self.CVCr_fs
        param[15] = self.q_lCr
        param[16] = self.D_lCr
        param[17] = self.CtotM_mf
        param[18] = self.q_mCr
        param[19] = self.decay_length
        param[20] = self.DV

        return param

    def sensitivity_analysis(self, h, n, t, y, J_v0, nom):
        """
        """
        # Number of parameters - 1
        nn= 21
        param = zeros(nn)
        y_old = zeros(n)
        y_new = zeros(n)
        sens = zeros(n)
        sensitivity = ''
        temp_name = ''
        param = self.load_params(param)
        t, y, J_v0 = self.time_integration_RK2_Leistner_2012(h, n, t, y, J_v0)

        for i in range(n):
            y_old[i] = y[i]

        temp_name += nom[:len(nom)-7]
        temp_name += '.sensitivity_analysis'
        sensitivity = temp_name

        with open(sensitivity, 'w') as f:
            for i in range(nn):
                # increase the value of parameter i by 10%
                param[i] = param[i]*1.1
                self.recopie_xi_param(param)
                self.time_integration_RK2_Leistner_2012(h, n, t, y, J_v0)
                # set all parameters back to original values
                param[i] = param[i]/1.1
                f.write('{}i'.format(i))
                f.write('\t')
                for j in range(n):
                    y_new[j] = y[j]
                    # //sensitivity of film width with respect to parameter i,
                    # at time j
                    sens[j] = (param[i]/y_old[j]) *\
                    (y_new[j] - y_old[j])/(param[i]*1.1 - param[i])

                    if j % 360 == 0:
                        f.write('%.3e'%sens[j])
                        f.write('\t')
                f.write('\n')

    def history_synthesis(self, n, t, y, step):
        """
        """
        string = ''
        result = 't[s], X_calc[m], X_exp[m]'

        for i in range(n+1):
            if i%step == 0:
                print(string, "%.12e, %.12e, %.12e"%(t[i], y[i],
                                                     self.fit_exp(t[i])))
                result += string
        return result

    def affiche(self, n, t, y):
        """
        """
        print("t[s], X_calc[m], X_exp[m]")
        for i in range(n+1):
            if i%3600 == 0:
                print("%.12e, %.12e, %.12e"%(t[i], y[i],
                                             self.fit_exp(t[i])))
        print()

    def write_output(self, nom, n, t, y, J_v0, J_H):
        """
        """
        with open(nom, 'w') as f:
            f.write('t[s], X_calc[m], X_exp[m], J_v0[1/m2/s], J_H[1/m2/s]\n')
            for i in range(n+1):
                if i%360 == 0:
                    f.write('%.12e, %.12e, %.12e, %.12e, %.12e\n'%(t[i], y[i],
                            self.fit_exp(t[i]), J_v0[i], J_H[i]))
            f.write('\n')

        # Write output for graphing
        temp_name = ''
        temp_name += nom[:len(nom)+7]
        graph_fig = temp_name
        temp_name += '.graphs'
        graph_fig = temp_name

        with open(graph_fig, 'w') as f:
            f.write('t[s], X_calc[m], X_exp[m], J_vO[1/m2/s], J_vO_x[1/m/s], ')
            f.write('f(X_calc)], J_dissol[1/m2/s], J_vO-J_dissol, ')
            f.write('J_MCr, J_ICr, J_H[/m2.s]')
            for i in range(n+1):
                if i%3600 == 0:
                    f.write('%.12e, %.12e, '%(t[i], y[i]))
                    f.write('%.12e, %.12e, '%(self.fit_exp(t[i]),
                                              self.J_vO_Leistner_2012(y[i])))
                    f.write('%.12e, %.12e, '%(y[i]*J_v0[i], self.f_(y[i])))
                    f.write('%.12e, '%self.J_dissol_Leistner_2012(y[i]))
                    f.write('%.12e, '%(self.J_vO_Leistner_2012(y[i])-
                                       self.J_dissol_Leistner_2012(y[i])))
                    f.write('%.12e, %.12e, '%(self.J_MCr_Leistner_2012(y[i]),
                                              self.J_ICr_Leistner_2012(y[i])))
                    f.write('%.12e\n'%J_H[i])
            f.write('\n')

    def write_potential(self, nom, n, y, t, d, F_0mf, F_0f, F_0fs, alpha, DV):
        """
        total concentration of H
        """
        temp_name = ''
        temp_name += nom[:len(nom)-7]
        potential = temp_name
        temp_name += '.potential'
        potential = temp_name

        with open(potential, 'w') as f:
            f.write('$t$[s], $\\phi^{\\circ}_{f/s}(d_{b})$, ')
            f.write('$\\phi_{f}(d_{b})$, $\\phi_{f/s}(x)$, $\\phi_{f}(x)$, ')
            f.write('$V(d_b)$, $V(x)$\n')
            for i in range(n+1):
                if i%360 == 0:
                    DF_fs = alpha * DV * exp(-1.*(y[i]/d))
                    DF_f = DV *(1.-alpha *  exp(-1.*(y[i]/d)))
                    F_fs = F_0fs + alpha * DV * exp(-1.*(y[i]/d))
                    F_f = F_0f + DV *(1.-alpha *  exp(-1.*(y[i]/d)))
                    V = F_0mf + (F_0f + DF_f) + (F_0fs + DF_fs)
                    V0 = F_0mf + F_0f + F_0fs
                    f.write('%.12e, %.12e, %.12e, %.12e, '%(t[i], F_0fs,
                                                            F_0f, F_fs))
                    f.write('%.12e, %.12e, %.12e\n'%(F_f, V0, V))
            f.write('\n')


    def pause_clavier(self):
        """
        """
        input("Appuyez sur une touche pour continuer")

    def affiche_param(self, x, NP):
        """
        """
        print("Parametres :")
        self.affiche_vect_scientifique(x, NP)
        self.affiche_vect(x, NP)
        print()

    def affiche_vect(self, x, n):
        """
        """
        # for (i=0;i<=n;i++)
        for i in range(n+1):
            print("%e "%x[i])
            if (i<n):
                print(", ")
        print()

    def affiche_vect_scientifique(self, x, n):
        """
        """
        # for (i=1;i<=n/2;i++):
        for i in range(int(n/2)+1):
            print("alpha%d=%e"%(i, x[2*(i-1)+2]))
            print("d%d=%e"%(i, x[2*(i-1)+1]))
        print()

    def ivector(self, n):
        """
        """
        # return (int *) malloc(n*sizeof(int))
        return zeros(n)

    def dvector(self, n):
        """
        """
        # return (double *) malloc(n*sizeof(double))
        return zeros(n)

    def dmatrix(self, row, col):
        """
        """
        return zeros((row, col))

    
    def optimisation(self):
        """
        """
        # Nombre total de paramètres dans le modèle de croissance
        nn= 19

        # Les NP+1 premiers parametres sont optimisees
        NP = 6

        # ALGO=1 : Plus grande pente
        # ALGO=2 : Fletcher-Reeves
        # ALGO=3 : Polak-Ribiere
        ALGO = 2

        algo = ["", "Plus grande pente", "Fletcher-Reeves", "Polak-Ribiere"]

        # Nombre maximum d'itérations
        MAXITER = 10

        # Pas des différences finies
        H = 1e-4

        # Tolérance sur le gradient
        GTOL = 1e-20

        # Tolérance sur l'erreur
        ETOL = 1e-19

        # Intervalle de temps de la simulation
        nb_heures = 10

        # Nombre de points d'intégration en s
        NINTG = 3600*(nb_heures)

        # Vecteur gradient et auxiliaire
        t = zeros(3600*(nb_heures+1))
        y = zeros(3600*(nb_heures+1))
        
        # print parameters
        print("Nombre de paramètres           = %d"%NP)
        print("Algorithme de descente         = %s"%algo[ALGO])
        print("Nombre maximum d'itérations    = %d"%MAXITER)
        print("Pas des différences finies     = %.1e"%H)
        print("Tolérance sur le gradient      = %.1e"%GTOL)
        print("Tolérance sur l'erreur         = %.1e"%ETOL)
        print("Nombre de points d'intégration = %d"%int(NINTG))

        # allocation des tableaux utilisés pour la résolution
        opt = self.ivector(nn+1)
        x = self.dvector(nn+1)
        z = self.dvector(nn+1)
        dx = self.dvector(nn+1)
        ddx = self.dvector(nn+1)

        # On commence par résoudre l'equation diff. avec les valeurs de départ
        # des paramètres et on sauvegarde le resultat dans un fichier
        # texte avec les valeurs attendues d'après le fit experimental
        t[0] = 0
        y[0] = 1e-19
        h = 1
        self.time_integration_RK2_Seyeux_2010(h, NINTG, t, y)
        self.affiche(NINTG, t, y)
        self.sauver("c_solver_corrosion_evolution_test_1.output", NINTG, t, y)

        # On calcule et on affiche la valeur de l'erreur et du gradient
        # pour les valeurs initiales   des paramètres
        self.recopie_param_xi(x, opt)
        self.affiche_param(x, NP);

        # test de erreur()
        print("Valeurs de l'erreur avec les paramètres de départ e=%e"%
              self.erreur(NINTG, h, x, t, y))

        # test de gradient()
        print("Valeurs du gradient avec les paramètres de départ g=%e"%
              self.gradient(H, NP, NINTG, h, x, dx, t, y))
        self.affiche_vect(dx, NP)

        # On optimise les paramètres avec la méthode du gradient conjugué
        # et on sauvegarde le resultat dans un fichier  texte avec les
        # valeurs attendues d'après le fit experimental
        self.gradient_conjugue(NP, ALGO, MAXITER, H, GTOL, ETOL, NINTG, x, dx,
                               ddx, t, y)
        self.affiche_param(x, NP)
        self.sauver("c_solver_corrosion_evolution_optimisation_1.output",
                    NINTG, t, y)

    def gradient_conjugue(self, NP, ALGO, MAXITER, H, GTOL, ETOL, NINTG,
                          x, dx, ddx, t, y):
        """
        """
        h = 1.
        e_old = 0
        print("debut gradient_conjugue")

        # La première direction, c'est le gradient [dEx, dEy, dEz]
        # g est le carré de sa norme
        g = self.gradient(NP, H, NINTG, h, x, dx, t, y)

        for n in range(1, MAXITER+1):
            # Récupère le minimum d'énergie dans la direction [dEx, dEy, dEz]
            e = self.descendre(NP, NINTG, h, x, dx, t, y)
            for i in range(NP+1):
                # Copie l'ancienne direction [dEx, dEy, dEz]
                # dans [dGx, dGy, dGz]
                ddx[i] = dx[i]
            # Calcule le nouveau gradient [dEx, dEy, dEz] et le carré de
            # sa norme
            h = self.gradient(NP, H, NINTG, h, x, dx, t, y)

            # Affiche [itération, énergie, gradient]
            print("n = %5d, e = %.10e, g = %e"%(n, e, sqrt(h/(NP))))
            self.affiche_param(x, NP)

            # ompare le gradient actuel au gradient terminal souhaité (GTOL)
            if sqrt(h/(NP))<GTOL:
                break

            # Pour Fletcher-Reeves on a tout : h = dE*dE
            if n > 1 and abs(e_old-e)<ETOL:
                print("Energie ne varie plus")
                break
            if ALGO == 3:
                # Pour Polak-Ribiere il faut calculer h = (dE-dG)*dG
                h = 0
                for i in range(NP+1):
                    h += (dx[i]-ddx[i])*dx[i]

            # Facteur de pondération entre ancienne direction et
            # nouveau gradient
            b = h/g

            # Il reste à calculer la nouvelle direction et le carré de sa norme
            g = 0

            for i in range(NP+1):
                dx[i] += b*ddx[i]
                g += dx[i]*dx[i]
            e_old = e.copy()
        print("fin gradient_conjugue")
        return e

    def gradient(self, NP, H, NINTG, h, x, dx, t, y):
        """
        """
        for i in range(NP+1):
            t_temp = x[i]
            x[i] = t_temp*(1+H)
            I1 = self.erreur(NINTG, h, x, t, y)
            x[i] = t_temp*(1-H)
            I2 = self.erreur(NINTG, h, x, t, y)
            x[i] = t_temp
            dx[i] = (I1-I2)/(2*H)
        g = 0
        for i in range(NP):
            g += dx[i]*dx[i]
        return g

    def erreur(self, NINTG, h, x, t, y):
        """
        """
        self.recopie_xi_param(x)
        # k = 1
        t[0]=0
        y[0]=1e-19
        h=1
        self.time_integration_RK2_Seyeux_2010(h, NINTG, t, y)

        # Calcul de l'integrale par la méthode des trapèzes
        I = self.f1(1, y) + self.f1(NINTG, y)
        for i in range(2, NINTG):
            for i in range(3600, NINTG, 3600):
                I = I + (2+2*(i%2))*self.f1(i, y)
        I = (h/3)*I
        return I

    def f1(self, i, y):
        """
        """
        return (y[i]-self.eqdf.fit_exp(i))*(y[i]-self.eqdf.fit_exp(i))*1e9

    def descendre(self, NP, NINTG, h, x, dx, t, y):
        """
        """
        # pas de départ de la descente
        R0 = 1e-6

        # La direction de descente est déjà dans [dEx, dEy, dEz]
        e = self.erreur(NINTG, h, x, t, y)

        # Recherche r, le pas de départ, en commençant à R0
        r = R0
        f = e

        # Tant que l'énergie augmente, il faut réduire le pas
        ### do-while in original code ###
        while e < f:
            # Avance
            self.progresser(NP, r, x, dx)
            e = self.erreur(NINTG, h, x, t, y)
            if e>f:
                # Recule
                self.progresser(NP, -r, x, dx)
                # Divise le pas par 2
                r /= 2
        ### do-while in original code ###

        if r < R0:
            # Informe que le pas initial est inférieur à R0
            print("r<R0 (%.1e) : %e"%(R0, r))
        # Le pas de départ est trouvé
        # Tant que l'énergie baisse on avance
        ### do-while in original code ###
        while e < f:
            f = e
            # Avance
            self.progresser(NP, r, x, dx)
            e = self.erreur(NINTG, h, x, t, y)
            # Si l'énergie baisse on augmente le pas
            if e < f:
                r *= 2
        ### do-while in original code ###

        # On est allé trop loin, il faut reculer
        self.progresser(NP, -r, x, dx)
        # Retourne la meilleure énergie (c'est f, pas e)
        return f

    def progresser(self, NP, r, x, dx):
        """
        """
        for i in range(NP+1):
            x[i] -= r*dx[i]*x[i]

    def run(self):
        '''
        '''
        NINTG = 36000
        t = zeros(NINTG+1)
        y = zeros(NINTG+1)
        J_v0 = zeros(NINTG+1)
        J_MCr = zeros(NINTG+1)
        J_ICr = zeros(NINTG+1)
        J_H = zeros(NINTG+1)
        t[0] = 0
        y[0] = 1e-12
        h = self.time_in_seconds/(1.*NINTG)
        if 'Seyeux_2010' == self.model:
            t, y = self.time_integration_RK2_Seyeux_2010(h, NINTG, t, y)
            if 'True' == self.save_history:
                print(self.history_synthesis(NINTG, t, y, 3600))
                self.sauver(self.output_file_name, NINTG, t, y)
            else:
                print("y(t = %.6e) = \n%.6e\n"%(t[NINTG], y[NINTG]))
    
        if 'Leistner_2012' in self.model:
            t, y, J_v0 = self.time_integration_RK2_Leistner_2012(h, NINTG,
                                                                 t, y, J_v0)
            if 'True' == self.save_history:
                print(self.history_synthesis(NINTG, t, y, 3600))
                self.write_output(self.output_file_name, NINTG,
                                  t, y, J_v0, J_H)
            else:
                print("y(t = %.6e) = \n%.6e"%(t[NINTG], y[NINTG], J_v0[NINTG]))
    
        if 'Voyshnis_2014' in self.model:
            t, y, J_v0, J_MCr, J_ICr, J_H =\
                self.time_integration_RK2_Voyshnis_2014(h, NINTG, t, y, J_v0,
                                                        J_MCr, J_ICr, J_H)
            return t, y
            '''
            if 'True' == self.save_history:
                self.history_synthesis(NINTG, t, y, 3600)
                self.write_output(self.output_file_name,
                                  NINTG, t, y, J_v0, J_H)
            '''


def plot(t, y):
    '''
    '''
    plt.plot(t[::3600][1:], y[::3600][1:],
             label='Thickness')
    plt.legend()
    plt.show()

def run():
    """
    """
    if not os.path.exists('Outputs/'):
        os.makedirs('Outputs/')
    os.chdir('Outputs/')
    code_name = 'c_solver_corrosion_evolution'

    # read parameters

    input_file_name = '../Inputs/c_solver_corrosion_evolution_test_8.input'
    ini = iniparser_load(input_file_name)
    oxydation = Oxydation(model=ini['model'],
                          save_history=ini['save_history'],
                          output_file_name=ini['output_file_name'],
                          temperature_in_K=ini['temperature_in_K'],
                          pH_temperature=ini['pH_temperature'],
                          time_in_seconds=ini['time_in_seconds'],
                          x_Cr=ini['x_Cr'], x_Fe=ini['x_Fe'], x_Ni=ini['x_Ni'],
                          T=ini['T'], pH=ini['pH'], alpha=ini['alpha'],
                          DV=ini['DV'], F_0f=ini['F_0f'], F_0mf=ini['F_0mf'],
                          F_0fs=ini['F_0fs'], CtotM_mf=ini['CtotM_mf'],
                          CtotI=ini['CtotI'], CmFe_mf=ini['CtotI'],
                          CMNi_mf=ini['CMNi_mf'], CVCr_mf=ini['CVCr_mf'],
                          CtotM_fs=ini['CtotM_fs'], CmFe_fs=ini['CmFe_fs'],
                          CMNi_fs=ini['CMNi_fs'], CVCr_fs=ini['CVCr_fs'],
                          D_vO=ini['D_vO'], D_mCr=ini['D_mCr'],
                          D_ICr=ini['D_ICr'], D_lCr=ini['D_lCr'],
                          DG1=ini['DG1'], DG2=ini['DG2'], DG3=ini['DG3'],
                          DG4=ini['DG4'], DG5=ini['DG5'], DG6=ini['DG6'],
                          DG7=ini['DG7'], DG8=ini['DG8'], DG9=ini['DG9'],
                          DG10=ini['DG10'], DG11=ini['DG11'], DG12=ini['DG12'],
                          DG13=ini['DG13'], DG14=ini['DG14'],
                          decay_length=ini['decay_length'],
                          charge_number=ini['charge_number'],
                          dissol_order=ini['dissol_order'],
                          dissol_preexp=ini['dissol_preexp'],
                          dissol_Ea=ini['dissol_Ea'], q_vO=ini['q_vO'],
                          q_lCr=ini['q_lCr'], q_mCr=ini['q_mCr'],
                          gamma_H=ini['gamma_H'])

    line = "******************************************************************"
    print("%s\ncode : %s\n\nBEGIN\n\n"%(line, code_name))
    
    oxydation.print_parameters()

    oxydation.set_user_parameters()

    NINTG = 36000
    t = zeros(NINTG+1)
    y = zeros(NINTG+1)
    J_v0 = zeros(NINTG+1)
    J_MCr = zeros(NINTG+1)
    J_ICr = zeros(NINTG+1)
    J_H = zeros(NINTG+1)
    t[0] = 0
    y[0] = 1e-12
    h = oxydation.time_in_seconds/(1.*NINTG)

    if 'Seyeux_2010' == oxydation.model:
        t, y = oxydation.time_integration_RK2_Seyeux_2010(h, NINTG, t, y)
        '''
        if 'True' == oxydation.save_history:
            print(oxydation.history_synthesis(NINTG, t, y, 3600))
            oxydation.sauver(oxydation.output_file_name, NINTG, t, y)
        else:
            print("y(t = %.6e) = \n%.6e\n"%(t[NINTG], y[NINTG]))
        '''

    if 'Leistner_2012' in oxydation.model:
        t, y, J_v0 = oxydation.time_integration_RK2_Leistner_2012(h, NINTG,
                                                                  t, y, J_v0)
        '''
        oxydation.write_potential(oxydation.output_file_name, NINTG, y, t,
                                  oxydation.decay_length, oxydation.F_0mf,
                                  oxydation.F_0f, oxydation.F_0fs,
                                  oxydation.alpha, oxydation.DV)
        if 'True' == oxydation.save_history:
            print(oxydation.history_synthesis(NINTG, t, y, 3600))
            oxydation.write_output(oxydation.output_file_name, NINTG,
                                   t, y, J_v0, J_H)
        else:
            print("y(t = %.6e) = \n%.6e"%(t[NINTG], y[NINTG], J_v0[NINTG]))
        oxydation.sensitivity_analysis(h, NINTG, t, y, J_v0,
                                       oxydation.output_file_name)
        '''

    if 'Voyshnis_2014' in oxydation.model:
        t, y, J_v0, J_MCr, J_ICr, J_H =\
            oxydation.time_integration_RK2_Voyshnis_2014(h, NINTG, t, y, J_v0,
                                                         J_MCr, J_ICr, J_H)
        
        plot(t/3600, y*1e9)
        oxydation.write_potential(oxydation.output_file_name, NINTG, y, t,
                                  oxydation.decay_length, oxydation.F_0mf,
                                  oxydation.F_0f, oxydation.F_0fs,
                                  oxydation.alpha, oxydation.DV)
        
        if 'True' == oxydation.save_history:
            print(oxydation.history_synthesis(NINTG, t, y, 3600))
            oxydation.write_output(oxydation.output_file_name,
                                   NINTG, t, y, J_v0, J_H)
        else:
            print("y(t = %.6e) = \n%.6e"%(t[NINTG], y[NINTG], J_v0[NINTG]))
        oxydation.sensitivity_analysis(h, NINTG, t, y, J_v0,
                                       oxydation.output_file_name)
        
    return t, y

def test(model='Voyshnis_2014'):
    Cr1 = 0.64
    Cr2 = 0.32
    
    oxy1 = Oxydation(x_Cr=Cr1, x_Ni=0.9-Cr1, model=model)
    oxy2 = Oxydation(x_Cr=Cr2, x_Ni=0.9-Cr2, model=model)
    print(oxy1.J_ICr_Leistner_2012(1e-12))
    print(oxy2.J_ICr_Leistner_2012(1e-12))
    return
    t, y1 = run(Oxydation(x_Cr=Cr1, x_Ni=0.9-Cr1, model=model))
    _, y2 = run(Oxydation(x_Cr=Cr1, x_Ni=0.9-Cr2, model=model))
    y1 = y1 * 1e9
    y2 = y2 * 1e9
    t = t/3600
    plt.plot(t[::3600][1:], y1[::3600][1:],
             label='Cr = {}'.format(Cr1))
    plt.plot(t[::3600][1:], y2[::3600][1:],
             label='Cr = {}'.format(Cr2))
    plt.xlabel('Time (h)')
    plt.ylabel('Thickness (nm)')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    #run()
    test(model='Voyshnis_2014')
