# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 08:50:29 2020

@author: Chao PAN

Simulation for RIS with model MIK (Modified Inverse Kirkendall)
Ternary Alloys Ee-Cr-Ni

This program calculates the amount of radiation induced segregation for a
ternary concentrated alloy. The formulation is based on the perks model
and is solved numerically using the gear subroutines.
"""


import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from time import time


class MIK:
    """
    Simulation for RIS with model MIK (Modified Inverse Kirkendall)
    Ternary Alloys Ee-Cr-Ni
    """
    def __init__(self,
                 CONCA=0.7, CONCB=0.21, CONCC=0.09,
                 R1=4., R2=18., RF=2018.,
                 N1=16, N2=14, N3=20,
                 ECOHA=-4.28, ECOHB=-4.21, ECOHC=-4.44,
                 EORDAB=0.003, EORDAC=-0.001, EORDBC=0.005,
                 EFA=1.4, EFB=1.6, EFC=1.79,
                 EMA=1.28, EMB=0.97, EMC=1.04,
                 EMIA=0.9, EMIB=0.9, EMIC=0.9,
                 WAV=1.8, WBV=3.2, WCV=1.0,
                 FAV=0.785, FBV=0.668, FCV=0.872, FI=0.660,
                 WAI=1.0, WBI=1.0, WCI=1.0,
                 SV=1.0,
                 EFGB=1.4,
                 TEMPC=360.,
                 DOSES=[0.994],
                 DISPRT=1.3999999999999999e-6,
                 DISL=1e14,
                 ETAV=1., ETAI=1.,
                 BIASV=1.0, BIASI=1.0,
                 EPS=1e-9,
                 NAT=9.1e28,
                 BOLTZ=8.6170002759899944e-5,
                 LAMBDA=3.4999999999999998e-10,
                 Z=12,
                 AL=1,
                 NUOV=1.5e13, NUOI=1.5e12,
                 SCFAC=9.9999997171806854E-010,
                 CAFRAC=None, CBFRAC=None, CCFRAC=None,
                 TFRAC=None, DFRAC=None, SFRAC=None,
                 Cr_profile=None):
        """
        Initiqlize all variables and set all parameters

        Initialise model: A=Fe, B=Cr, C=Ni

        CONCA, CONCB, CONCC: fractional concentration
            = 0.7, 0.21, 0.09

        R1, R2, RF: distance to grain boundary (nm)
            = 4.0, 18.0, 2018.0 (m)

        N1, N2, N3 : Number of points in mesh groups
            = 16, 14, 20

        ECOHA, ECOHB, ECOHC: Cohesive energies
            = -4.28, -4.10, -4.44 (eV)

        EORDAB, EORDAC, EORDBC: Ordering energies
            = 0.003, -0.001, 0.005 (eV)

        EFA, EFB, EFC: Pure element [vacancy] formation energies
            = 1.4, 1.6, 1.79 (eV)

        EMA, EMB, EMC: Pure element [vacancy] migration energies
            = 1.28, 0.97, 1.04 (eV)

        EMIA, EMIB, EMIC: Interstitial migration energies
            = 0.9, 0.9, 0.9 (eV)

        WAV, WBV, WCV: Relative vacancy jump frequency ratio
            = 1.8, 3.2, 1.0

        FAV, FBV, FCV, FI: Jump correlation factors
            = 0.785, 0.668, 0.872, 0.660

        WAI, WBI, WCI: Relative interstitial jump frequency ratio
            = 1.0, 1.0, 1.0

        SV: Vacancy formation enthalpy over BOLTZ
            = 1.0

        EFGB: Grain boundary formation energy
            = 1.4 (eV)

        TEMPC: Peak temperature
            = 360 (°C)
        
        DOSE: Dose
            = 0.994 (dpa)

        DISPRT: Peak displacement rate
            = 1.4e-6 (dpa/s)

        DISL: Peak dislocation density
            = 1e14 (m/m3)
        
        ETAV, ETAI: Vacancy/Interstitial production efficiency
            = 1.0, 1.0

        BIASV, BIASI: Dislocation bias for vacancy/interstitial
            = 1.0, 1.0

        EPS: Error tolerance
            = 1e-9

        NAT: Number of density
            = 9.1e28 (#/m3)

        BOLTZ: Boltzmann constant
            = 8.617e-5 (eV/K)

        LAMBDA: Jump distance
            = 3.5e-10 (m) 

        AL: Thermo factor
            = 1.0

        Z: Neighbor atoms
            = 12

        NUOV, NUOI: Debye frequencies
             = 1.5e13, 1.5e12 (/s)

        SCFAC: Nanometer
            = 1e-9 (m)

        CAFRAC: Fraction of max temperature
        CBFRAC, CCFRAC, TFRAC: Fraction of peak atoms
        DFRAC: Fraction of damage
        SFRAC: Fraction of max dislocation density
        """
        self.CONCA, self.CONCB, self.CONCC = CONCA, CONCB, CONCC
        self.ECOHA, self.ECOHB, self.ECOHC = ECOHA, ECOHB, ECOHC
        self.EORDAB, self.EORDAC, self.EORDBC = EORDAB, EORDAC, EORDBC
        self.EFA, self.EFB, self.EFC = EFA, EFB, EFC
        self.EMA, self.EMB, self.EMC = EMA, EMB, EMC
        self.EMIA, self.EMIB, self.EMIC = EMIA, EMIB, EMIC
        self.WAV, self.WBV, self.WCV = WAV, WBV, WCV
        self.FAV, self.FBV, self.FCV, self.FI = FAV, FBV, FCV, FI
        self.WAI, self.WBI, self.WCI = WAI, WBI, WCI
        self.SV = SV
        self.EFGB = EFGB
        self.TEMPC = TEMPC
        self.DOSES = DOSES
        self.DISPRT = DISPRT
        self.DISL = DISL
        self.ETAV, self.ETAI = ETAV, ETAI
        self.BIASV, self.BIASI = BIASV, BIASI
        self.EPS = EPS
        self.NAT = NAT
        self.BOLTZ = BOLTZ
        self.LAMBDA = LAMBDA
        self.Z = Z
        self.AL = AL
        self.NUOV, self.NUOI = NUOV, NUOI
        self.SCFAC = SCFAC

        # Space step
        self.NSTEP = N1+N2+N3

        # Pair interaction energy between like atoms
        self.EAA = ECOHA/(Z/2)
        self.EBB = ECOHB/(Z/2)
        self.ECC = ECOHC/(Z/2)
        # Pair interaction energy between unlike atoms
        self.EAB=0.5*(self.EAA+self.EBB)-EORDAB
        self.EAC=0.5*(self.EAA+self.ECC)-EORDAC
        self.EBC=0.5*(self.EBB+self.ECC)-EORDBC
        # Pair interaction energy for atoms and vacancy
        self.EAV = (ECOHA+EFA)/Z
        self.EBV = (ECOHB+EFB)/Z
        self.ECV = (ECOHC+EFC)/Z
        # Saddle point energy in the pure metal
        self.ESA = EMA+Z*(self.EAA+self.EAV)
        self.ESB = EMB+Z*(self.EBB+self.EBV)
        self.ESC = EMC+Z*(self.ECC+self.ECV)
        # Migration energy
        self.EA = np.zeros(self.NSTEP)
        self.EB = np.zeros(self.NSTEP)
        self.EC = np.zeros(self.NSTEP)

        # Geometry
        self.MESHSP = np.zeros(self.NSTEP)
        self.XVALUE = np.zeros(self.NSTEP)
        self.MESHSI = np.zeros(self.NSTEP)
        self.mesh(N1, N2, N3, R1, R2, RF, SCFAC)

        # Fraction of max temperature, peak atoms, max dislocation density,
        # damage
        self.CAFRAC = np.ones(self.NSTEP)
        self.CBFRAC = np.ones(self.NSTEP)
        self.CCFRAC = np.ones(self.NSTEP)
        self.TFRAC = np.ones(self.NSTEP)
        self.DFRAC = np.ones(self.NSTEP)
        self.SFRAC = np.ones(self.NSTEP)
        self.get_profile(CAFRAC, CBFRAC, CCFRAC, TFRAC, DFRAC, SFRAC)

        # Displacement of vacancy/interstitial
        self.DISPV = np.zeros(self.NSTEP)
        self.DISPV[:] = DISPRT*ETAV*self.DFRAC[:]
        self.DISPI = np.zeros(self.NSTEP)
        self.DISPI[:] = DISPRT*ETAI*self.DFRAC[:]
        self.DISLOC = np.zeros(self.NSTEP)
        self.DISLOC[:] = DISL*self.SFRAC[:]

        # Temperature
        self.TKT = np.zeros(self.NSTEP)
        self.TKT[:-1] = BOLTZ*(TEMPC+273)*self.TFRAC[:-1]

        # Concentrations
        self.CA = CONCA*self.CAFRAC[:]
        self.CB = CONCB*self.CBFRAC[:]
        self.CC = CONCC*self.CCFRAC[:]
        self.set_FeCrNi_profile(Cr_profile)
        self.CI = np.zeros(self.NSTEP)
        self.CVTHER = self.calculate_CVTHER()
        self.CV = self.CVTHER.copy()
        self.CERR = np.zeros(self.NSTEP)
        self.NA = np.zeros(self.NSTEP)
        self.NB = np.zeros(self.NSTEP)
        self.NC = np.zeros(self.NSTEP)
        self.NV = np.zeros(self.NSTEP)
        self.NI = np.zeros(self.NSTEP)

        # Pre-factor of vacancy jump frequence (/s))
        self.PREVA = NUOV*WAV*FAV
        self.PREVB = NUOV*WBV*FBV
        self.PREVC = NUOV*WCV*FCV
        # interstitial jump frequency (/s)
        self.NUIA = np.zeros(self.NSTEP)
        self.NUIA[:-1] = NUOI*WAI*FI*np.exp((-1*EMIA)/self.TKT[:-1])
        self.NUIB = np.zeros(self.NSTEP)
        self.NUIB[:-1] = NUOI*WBI*FI*np.exp((-1*EMIB)/self.TKT[:-1])
        self.NUIC = np.zeros(self.NSTEP)
        self.NUIC[:-1] = NUOI*WCI*FI*np.exp((-1*EMIC)/self.TKT[:-1])
        self.NUVA = np.zeros(self.NSTEP)
        self.NUVB = np.zeros(self.NSTEP)
        self.NUVC = np.zeros(self.NSTEP)

        # Diffusivities (m2/s)
        self.DAIO = np.zeros(self.NSTEP)
        self.DAIO[:-1] = 0.66667*self.NUIA[:-1]*LAMBDA**2
        self.DBIO = np.zeros(self.NSTEP)
        self.DBIO[:-1] = 0.66667*self.NUIB[:-1]*LAMBDA**2
        self.DCIO = np.zeros(self.NSTEP)
        self.DCIO[:-1] = 0.66667*self.NUIC[:-1]*LAMBDA**2
        self.DAI = np.zeros(self.NSTEP)
        self.DBI = np.zeros(self.NSTEP)
        self.DCI = np.zeros(self.NSTEP)
        self.DAV = np.zeros(self.NSTEP)
        self.DBV = np.zeros(self.NSTEP)
        self.DCV = np.zeros(self.NSTEP)
        # Diffusion coefficient (m2/s)
        self.DA = np.zeros(self.NSTEP)
        self.DB = np.zeros(self.NSTEP)
        self.DC = np.zeros(self.NSTEP)
        self.DV = np.zeros(self.NSTEP)
        self.DI = np.zeros(self.NSTEP)
        self.DIFV = np.zeros(self.NSTEP)
        self.DIFI = np.zeros(self.NSTEP)

        # Gradient of concentration
        self.GRADCA = np.zeros(self.NSTEP)
        self.GRADCB = np.zeros(self.NSTEP)
        self.GRADCC = np.zeros(self.NSTEP)
        self.GRADCV = np.zeros(self.NSTEP)
        self.GRADCI = np.zeros(self.NSTEP)

        # Coefficient of recombinaison
        self.RECA = np.zeros(self.NSTEP)
        self.RECB = np.zeros(self.NSTEP)
        self.RECC = np.zeros(self.NSTEP)
        self.RECOMB = np.zeros(self.NSTEP)
        self.INTSINK = np.zeros(self.NSTEP)
        self.VACSINK = np.zeros(self.NSTEP)
        self.VACSOUR = np.zeros(self.NSTEP)

        # Flux
        self.JA = np.zeros(self.NSTEP)
        self.JB = np.zeros(self.NSTEP)
        self.JC = np.zeros(self.NSTEP)
        self.JV = np.zeros(self.NSTEP)
        self.JI = np.zeros(self.NSTEP)
        self.JO = np.zeros(self.NSTEP)
        self.JA0 = 0.0
        self.JB0 = 0.0
        self.JC0 = 0.0
        self.DIVJA = np.zeros(self.NSTEP)
        self.DIVJB = np.zeros(self.NSTEP)
        self.DIVJC = np.zeros(self.NSTEP)
        self.DIVJV = np.zeros(self.NSTEP)
        self.DIVJI = np.zeros(self.NSTEP)

        # Initial value
        self.Y0 = np.zeros(5*self.NSTEP)
        self.Y0[:self.NSTEP] = self.CA
        self.Y0[self.NSTEP:2*self.NSTEP] = self.CB
        self.Y0[2*self.NSTEP:3*self.NSTEP] = self.CC
        self.Y0[3*self.NSTEP:4*self.NSTEP] = self.CV
        self.Y0[4*self.NSTEP:] = self.CI
        self.Y = self.Y0.copy()
        # dy/dt
        self.CADOT = np.zeros(self.NSTEP)
        self.CBDOT = np.zeros(self.NSTEP)
        self.CCDOT = np.zeros(self.NSTEP)
        self.CVDOT = np.zeros(self.NSTEP)
        self.CIDOT = np.zeros(self.NSTEP)
        self.YDOT = np.zeros(5*self.NSTEP)
        # Null list
        self.EBs = []
        self.ECs = []
        self.DBs = []
        self.DCs = []

    def get_profile(self, CAFRAC, CBFRAC, CCFRAC, TFRAC, DFRAC, SFRAC):
        """
        To receive an input profile

        Fraction of max temperature, peak atoms, damage, max dislocation
        density
        """
        NSTEP = self.NSTEP
        if CAFRAC is None:
            TFRAC = np.ones(NSTEP)
            TFRAC[-1] = 0
            CAFRAC = np.ones(NSTEP)
            CBFRAC = np.ones(NSTEP)
            CCFRAC = np.ones(NSTEP)
            DFRAC = np.ones(NSTEP)
            SFRAC = np.ones(NSTEP)
        else:
            _temp = np.zeros(NSTEP)
        
            _temp[:NSTEP-1] = TFRAC[:NSTEP-1]
            TFRAC = np.zeros(NSTEP)
            TFRAC[:NSTEP-1] = _temp[:NSTEP-1]
        
            _temp[:NSTEP] = CAFRAC[:NSTEP]
            CAFRAC = np.zeros(NSTEP)
            CAFRAC[:NSTEP] = _temp[:NSTEP]
        
            _temp[:NSTEP] = CBFRAC[:NSTEP]
            CBFRAC = np.zeros(NSTEP)
            CBFRAC[:NSTEP] = _temp[:NSTEP]
        
            _temp[:NSTEP] = CCFRAC[:NSTEP]
            CCFRAC = np.zeros(NSTEP)
            CCFRAC[:NSTEP] = _temp[:NSTEP]
        
            _temp[:NSTEP] = DFRAC[:NSTEP]
            DFRAC = np.zeros(NSTEP)
            DFRAC[:NSTEP] = _temp[:NSTEP]
        
            _temp[:NSTEP] = SFRAC[:NSTEP]
            SFRAC = np.zeros(NSTEP)
            SFRAC[:NSTEP] = _temp[:NSTEP]
        self.CAFRAC[:] = CAFRAC[:]
        self.CBFRAC[:] = CBFRAC[:]
        self.CCFRAC[:] = CCFRAC[:]
        self.TFRAC[:] = TFRAC[:]
        self.DFRAC[:] = DFRAC[:]
        self.SFRAC[:] = SFRAC[:]

    def set_FeCrNi_profile(self, Cr_profile):
        """
        This function is to produce a prifole Fe and a profile of Ni given a
        profile of Cr.
        """
        if Cr_profile is not None:
            xvalue = self.XVALUE*1e9
            for i in range(self.NSTEP):
                self.CB[i] = get_value_from_profile(xvalue[i], Cr_profile)*0.01
            self.CA[:] = 0.7
            self.CC[:] = 1 - self.CA - self.CB

    def calculate_CVTHER(self):
        """
        Perfect sink: defect concentrations at thermal equilibrium
        """
        SV = self.SV
        EFGB = self.EFGB
        NSTEP = self.NSTEP
        CVTHER = np.zeros(NSTEP)
        CVTHER[:-1] = np.exp(SV)*np.exp(-EFGB/self.TKT[:-1])
        CVTHER[-1] = CVTHER[-2]
        for i in range(1, NSTEP-1):
            CVTHER[i] = 0.5*(CVTHER[i]+CVTHER[i-1])
        return CVTHER

    def calculate_migration_energy(self):
        """
        Calculate migration energy of Fe, Cr, Ni
        """
        self.EA[:-1]=((self.ESA+self.ESA*self.NA[:-1]+self.ESB*self.NB[:-1]+
                       self.ESC*self.NC[:-1])/2)-\
                      ((self.Z*(self.NA[:-1]*self.EAA+self.NB[:-1]*self.EAB+
                                self.NC[:-1]*self.EAC+self.NV[:-1]*self.EAV))+
                       (self.Z*(self.NA[:-1]*self.EAV+self.NB[:-1]*self.EBV+
                                self.NC[:-1]*self.ECV)))

        self.EB[:-1]=((self.ESB+self.ESA*self.NA[:-1]+self.ESB*self.NB[:-1]+
                       self.ESC*self.NC[:-1])/2)-\
                      ((self.Z*(self.NA[:-1]*self.EAB+self.NB[:-1]*self.EBB+
                                self.NC[:-1]*self.EBC+self.NV[:-1]*self.EBV))+
                       (self.Z*(self.NA[:-1]*self.EAV+self.NB[:-1]*self.EBV+
                                self.NC[:-1]*self.ECV)))

        self.EC[:-1]=((self.ESC+self.ESA*self.NA[:-1]+self.ESB*self.NB[:-1]+
                       self.ESC*self.NC[:-1])/2)-\
                      ((self.Z*(self.NA[:-1]*self.EAC+self.NB[:-1]*self.EBC+
                                self.NC[:-1]*self.ECC+self.NV[:-1]*self.ECV))+
                       (self.Z*(self.NA[:-1]*self.EAV+self.NB[:-1]*self.EBV+
                                self.NC[:-1]*self.ECV)))

    def mesh(self, N1, N2, N3, R1, R2, RF, SCFAC):
        """
        In order to discretize PDE, it is necessary to have a well-mesh grid.
        R1, R2, RF denote a distance to a grain boundary (GB), then we have 3
        zones: [0, R1], [R1, R2], [R2, RF].
        N1, N2, N3 represent mesh size in the 3 zoens respectively.
        """
        self.MESHSP[:N1-1] = R1*SCFAC/N1
        self.MESHSP[N1-1:N1+N2-1] = (R2-R1)*SCFAC/N2
        self.MESHSP[N1+N2-1:N1+N2+N3-1] = (RF-R2)*SCFAC/N3

        for i in range(1, self.NSTEP):
            self.XVALUE[i] = self.XVALUE[i-1]+self.MESHSP[i-1]

    def calculate_dCv_dCi(self):
        '''
        Limite condition

        Increment per time step for vancy and Interstitial
        '''
        NAT = self.NAT
        BIASV = self.BIASV
        BIASI = self.BIASI
        self.CVDOT[0] = 0.0
        #self.CVDOT[0] = self.CVTHER[0]
        self.CIDOT[0] = 0.0
        self.CVDOT[1:] = -self.DIVJV[1:]/NAT-\
                         self.RECOMB[1:]*self.CV[1:]*self.CI[1:]-\
                         BIASV*self.VACSINK[1:]*self.CV[1:]+self.VACSOUR[1:]+\
                         self.DISPV[1:]
        self.CIDOT[1:] = -self.DIVJI[1:]/NAT-\
                         self.RECOMB[1:]*self.CV[1:]*self.CI[1:]-\
                         BIASI*self.INTSINK[1:]*self.CI[1:]+self.DISPI[1:]

    def fex(self, t, y):
        """
        Auxiliary function for solving PDE sytem

        fex = fy/dt, right-hand of PDEs
        """
        NSTEP = self.NSTEP
        LAMBDA = self.LAMBDA
        Z = self.Z
        AL = self.AL
        NAT = self.NAT

        self.Y[:] = y[:]
        self.CA[:] = self.Y[:NSTEP]
        self.CB[:] = self.Y[NSTEP:2*NSTEP]
        self.CC[:] = self.Y[2*NSTEP:3*NSTEP]
        self.CV[:] = self.Y[3*NSTEP:4*NSTEP]
        self.CI[:] = self.Y[4*NSTEP:]
        
        self.NA[:-1] = 0.5*(self.CA[1:]+self.CA[:-1])
        self.NB[:-1] = 0.5*(self.CB[1:]+self.CB[:-1])
        self.NC[:-1] = 0.5*(self.CC[1:]+self.CC[:-1])
        self.NV[:-1] = 0.5*(self.CV[1:]+self.CV[:-1])
        self.NI[:-1] = 0.5*(self.CI[1:]+self.CI[:-1])
        
        self.DAI[:-1] = self.DAIO[:-1]
        self.DBI[:-1] = self.DBIO[:-1]
        self.DCI[:-1] = self.DCIO[:-1]
        
        self.calculate_migration_energy()
        
        self.NUVA[:-1] = self.PREVA*np.exp(-self.EA[:-1]/self.TKT[:-1])
        self.NUVB[:-1] = self.PREVB*np.exp(-self.EB[:-1]/self.TKT[:-1])
        self.NUVC[:-1] = self.PREVC*np.exp(-self.EC[:-1]/self.TKT[:-1])
        # In coherence with Fortran code, NUIA, NUIB, NUIC are 0
        self.NUIA = np.zeros(self.NSTEP)
        self.NUIB = np.zeros(self.NSTEP)
        self.NUIC = np.zeros(self.NSTEP)
        
        self.DAV[:-1] = self.NUVA[:-1]*LAMBDA**2
        self.DBV[:-1] = self.NUVB[:-1]*LAMBDA**2
        self.DCV[:-1] = self.NUVC[:-1]*LAMBDA**2

        self.RECA[:-1] = (self.NUVA[:-1]+self.NUIA[:-1])*Z
        self.RECB[:-1] = (self.NUVB[:-1]+self.NUIB[:-1])*Z
        self.RECC[:-1] = (self.NUVC[:-1]+self.NUIC[:-1])*Z

        self.DIFV[:-1] = self.DAV[:-1]*self.NA[:-1]+\
                         self.DBV[:-1]*self.NB[:-1]+\
                         self.DCV[:-1]*self.NC[:-1]
        self.DIFI[:-1] = self.DAI[:-1]*self.NA[:-1]+\
                         self.DBI[:-1]*self.NB[:-1]+\
                         self.DCI[:-1]*self.NC[:-1]
        
        self.DIFV[-1] = self.DIFV[-2]
        self.DIFI[-1] = self.DIFI[-2]
        for i in range(1, NSTEP-1):
            self.DIFV[i] = 0.5*(self.DIFV[i]+self.DIFV[i-1])
            self.DIFI[i] = 0.5*(self.DIFI[i]+self.DIFI[i-1])        

        self.RECA[-1] = self.RECA[-2]
        self.RECB[-1] = self.RECB[-2]
        self.RECC[-1] = self.RECC[-2]
        self.CVTHER[-1] = self.CVTHER[-2]
        for i in range(1, NSTEP-1):
            self.RECA[i] = 0.5*(self.RECA[i]+self.RECA[i-1])
            self.RECB[i] = 0.5*(self.RECB[i]+self.RECB[i-1])
            self.RECC[i] = 0.5*(self.RECC[i]+self.RECC[i-1])
            self.CVTHER[i] = 0.5*(self.CVTHER[i]+self.CVTHER[i-1])

        self.JA0 = self.JB0 = self.JC0 = 0.0
        self.JA[-1] = self.JB[-1] = self.JC[-1] = 0.0

        self.GRADCA[:-1] = (self.CA[1:]-self.CA[:-1])/self.MESHSP[:-1]
        self.GRADCB[:-1] = (self.CB[1:]-self.CB[:-1])/self.MESHSP[:-1]
        self.GRADCC[:-1] = (self.CC[1:]-self.CC[:-1])/self.MESHSP[:-1]
        self.GRADCV[:-1] = (self.CV[1:]-self.CV[:-1])/self.MESHSP[:-1]
        self.GRADCI[:-1] = (self.CI[1:]-self.CI[:-1])/self.MESHSP[:-1]

        self.DA[:-1] = self.DAV[:-1]*self.NV[:-1]+self.DAI[:-1]*self.NI[:-1]
        self.DB[:-1] = self.DBV[:-1]*self.NV[:-1]+self.DBI[:-1]*self.NI[:-1]
        self.DC[:-1] = self.DCV[:-1]*self.NV[:-1]+self.DCI[:-1]*self.NI[:-1]
        self.DV[:-1] = self.DAV[:-1]*self.NA[:-1]+self.DBV[:-1]*self.NB[:-1]+\
                       self.DCV[:-1]*self.NC[:-1]
        self.DI[:-1] = self.DAI[:-1]*self.NA[:-1]+self.DBI[:-1]*self.NB[:-1]+\
                       self.DCI[:-1]*self.NC[:-1]

        self.EBs.append(self.EB)
        self.ECs.append(self.EC)
        self.DBs.append(self.DB)
        self.DCs.append(self.DC)

        self.JV[:-1] = NAT*(-self.DV[:-1]*self.GRADCV[:-1]+
                            self.NV[:-1]*AL*(self.DAV[:-1]*self.GRADCA[:-1]+
                                             self.DBV[:-1]*self.GRADCB[:-1]+
                                             self.DCV[:-1]*self.GRADCC[:-1]))
        self.JI[:-1] = NAT*(-self.DI[:-1]*self.GRADCI[:-1]-
                            self.NI[:-1]*AL*(self.DAI[:-1]*self.GRADCA[:-1]+
                                             self.DBI[:-1]*self.GRADCB[:-1]+
                                             self.DCI[:-1]*self.GRADCC[:-1]))
        self.JO[:-1] = self.JI[:-1]-self.JV[:-1]
        self.JA[:-1] = NAT*(-self.DA[:-1]*AL*self.GRADCA[:-1]+
                            self.NA[:-1]*(self.DAV[:-1]*self.GRADCV[:-1]-
                                          self.DAI[:-1]*self.GRADCI[:-1]))-\
                       self.JO[:-1]*self.NA[:-1]
        self.JB[:-1] = NAT*(-self.DB[:-1]*AL*self.GRADCB[:-1]+
                            self.NB[:-1]*(self.DBV[:-1]*self.GRADCV[:-1]-
                                          self.DBI[:-1]*self.GRADCI[:-1]))-\
                       self.JO[:-1]*self.NB[:-1]
        self.JC[:-1] = NAT*(-self.DC[:-1]*AL*self.GRADCC[:-1]+
                            self.NC[:-1]*(self.DCV[:-1]*self.GRADCV[:-1]-
                                          self.DCI[:-1]*self.GRADCI[:-1]))-\
                       self.JO[:-1]*self.NC[:-1]
        self.JV[:-1] = self.JV[:-1]-self.JO[:-1]*self.NV[:-1]
        self.JI[:-1] = self.JI[:-1]-self.JO[:-1]*self.NI[:-1]

        self.DIVJA[0] = 2.0*(self.JA[0]-self.JA0)/self.MESHSP[0]
        self.DIVJB[0] = 2.0*(self.JB[0]-self.JB0)/self.MESHSP[0]
        self.DIVJC[0] = 2.0*(self.JC[0]-self.JC0)/self.MESHSP[0]

        self.MESHSI[1:-1] = 0.5*(self.MESHSP[1:-1]+self.MESHSP[:-2])
        self.DIVJA[1:-1] = (self.JA[1:-1]-self.JA[:-2])/self.MESHSI[1:-1]
        self.DIVJB[1:-1] = (self.JB[1:-1]-self.JB[:-2])/self.MESHSI[1:-1]
        self.DIVJC[1:-1] = (self.JC[1:-1]-self.JC[:-2])/self.MESHSI[1:-1]
        self.DIVJV[1:-1] = (self.JV[1:-1]-self.JV[:-2])/self.MESHSI[1:-1]
        self.DIVJI[1:-1] = (self.JI[1:-1]-self.JI[:-2])/self.MESHSI[1:-1]

        self.DIVJA[-1] = 2.0*(self.JA[-1]-self.JA[-2])/self.MESHSP[-2]
        self.DIVJB[-1] = 2.0*(self.JB[-1]-self.JB[-2])/self.MESHSP[-2]
        self.DIVJC[-1] = 2.0*(self.JC[-1]-self.JC[-2])/self.MESHSP[-2]
        self.DIVJV[-1] = 2.0*(self.JV[-1]-self.JV[-2])/self.MESHSP[-2]
        self.DIVJI[-1] = 2.0*(self.JI[-1]-self.JI[-2])/self.MESHSP[-2]

        self.CADOT[:] = -self.DIVJA[:]/NAT
        self.CBDOT[:] = -self.DIVJB[:]/NAT
        self.CCDOT[:] = -self.DIVJC[:]/NAT

        self.RECOMB[:] = self.RECA[:]*self.CA[:]+self.RECB[:]*self.CB[:]+\
                         self.RECC[:]*self.CC[:]
        self.INTSINK[:] = self.DISLOC[:]*self.DIFI[:]
        self.VACSINK[:] = self.DISLOC[:]*self.DIFV[:]
        self.VACSOUR[:] = self.DISLOC[:]*self.DIFV[:]*self.CVTHER[:]

        self.calculate_dCv_dCi()

        self.Y[:NSTEP] = self.CA[:]
        self.Y[NSTEP:2*NSTEP] = self.CB[:]
        self.Y[2*NSTEP:3*NSTEP] = self.CC[:]
        self.Y[3*NSTEP:4*NSTEP] = self.CV[:]
        self.Y[4*NSTEP:] = self.CI[:]
        self.YDOT[:NSTEP] = self.CADOT[:]
        self.YDOT[NSTEP:2*NSTEP] = self.CBDOT[:]
        self.YDOT[2*NSTEP:3*NSTEP] = self.CCDOT[:]
        self.YDOT[3*NSTEP:4*NSTEP] = self.CVDOT[:]
        self.YDOT[4*NSTEP:] = self.CIDOT[:]
        return self.YDOT

    def solve(self, print_time_used=False, time_interval=[0, 1]):
        """
        Solve PDE system
        """
        t1 = time()
        if self.DOSES is None:
            t = [time_interval[0], time_interval[-1]]
            t_eval = time_interval
        else:
            t = [0, self.DOSES[-1]/self.DISPRT]
            t_eval = np.array(self.DOSES)/self.DISPRT
        y = np.array(self.CA.tolist() + self.CB.tolist() + self.CC.tolist() +
                     self.CV.tolist() + self.CI.tolist())
        sol = solve_ivp(self.fex, t, y,
                        method='LSODA',
                        atol=self.EPS,
                        t_eval=t_eval)
        t2 = time()
        assert sol.success, 'fail to solve'
        if print_time_used:
            print("Finished with {}s".format(round(t2-t1, 2)))
        return sol

    def output(self):
        """
        """
        NSTEP = self.NSTEP
        self.CA = self.Y[:NSTEP]
        self.CB = self.Y[NSTEP:2*NSTEP]
        self.CC = self.Y[2*NSTEP:3*NSTEP]
        self.CV = self.Y[3*NSTEP:4*NSTEP]
        self.CI = self.Y[4:NSTEP]
        self.CERR = 1-(self.CA+self.CB+self.CC)
        self.XOUT = self.XVALUE*1e9

        TEMP1 = (self.CA*np.exp(-self.XOUT/0.8452)).sum()/\
            np.exp(-self.XOUT/0.8452).sum()
        TEMP2 = (self.CB*np.exp(-self.XOUT/0.7474)).sum()/\
            np.exp(-self.XOUT/0.7474).sum()
        TEMP3 = (self.CC*np.exp(-self.XOUT/0.9472)).sum()/\
            np.exp(-self.XOUT/0.9472).sum()
        CASURF = TEMP1/(TEMP1+TEMP2+TEMP3)
        CBSURF = TEMP2/(TEMP1+TEMP2+TEMP3)
        CCSURF = TEMP3/(TEMP1+TEMP2+TEMP3)
        return CASURF, CBSURF, CCSURF

    def convolution(self, FWTM=2, a=5, coef_steps_gaussian=20,
                    is_plot=False):
        """
        Post processing

        (g * h)(x) = ∫ g(x')h(x-x')dx'

        FWTM (full width at tenth of maximum) = 2*√(2*ln10)*c
        where c is standard deviation
        
        Interval [-a, a] for gaussian fonction g(x) (x is distance in nm)

        coef_steps_gaussian is an coefficient of steps within 1nm for gaussian
        function.
        """
        dist = np.concatenate((-self.XVALUE[::-1], self.XVALUE[1:]),
                              axis=0)*1e9

        sigma = FWTM/(2.*np.sqrt(2*np.log(10)))
        mu = 0
        gaussian = lambda x: 1/(sigma * np.sqrt(2 * np.pi)) *\
            np.exp( - (x - mu)**2 / (2 * sigma**2) )

        x = np.linspace(-a, a, coef_steps_gaussian*a+1)
        x = dist
        gx = gaussian(x)
        gx = gx/gx.sum()
        hx_Cr = np.concatenate((self.CB[::-1], self.CB[1:]), axis=0)
        hx_Ni = np.concatenate((self.CC[::-1], self.CC[1:]), axis=0)
        ghx_Cr = np.convolve(gx, hx_Cr)
        ghx_Ni = np.convolve(gx, hx_Ni)
        assert ghx_Cr.shape[0] == ghx_Ni.shape[0], 'check convolution'
        center_index = int(ghx_Cr.shape[0]/2)
        gh_Cr = ghx_Cr[center_index-self.NSTEP+1:center_index+self.NSTEP]
        gh_Ni = ghx_Ni[center_index-self.NSTEP+1:center_index+self.NSTEP]

        if is_plot:
            fig, axs  = plt.subplots(1, 2, figsize=(13, 5))
            substances = ['Cr', 'Ni']
            hxs = [hx_Cr[15:-15], hx_Ni[15:-15]]
            ghs = [gh_Cr[15:-15], gh_Ni[15:-15]]
            dist = dist[15:-15]
            for ax, substance, hx, gh in zip(axs, substances, hxs, ghs):
                ax.plot(dist, hx, linewidth=1.0,
                    linestyle='solid', label='Before convolution')
                ax.plot(dist, gh, linewidth=1.0,
                    linestyle='solid', label='After convolution')
                ax.legend()
                ax.set_title(substance)
                ax.set_xlabel('Distance to grain boundary (nm)')
                ax.set_ylabel('Fractional concentration')
            if not os.path.exists('images/'):
                os.makedirs('images/')
            plt.savefig('images/convolution.jpg')
            plt.show()
        
        return ghx_Cr[center_index], ghx_Ni[center_index]

    def plot(self, ax, x, substance, err, sol, color, with_test,
             is_in_linear_scale, xticks, xlabels,
             labels=['Python', 'Fortran']):
        ''''''
        python, fortran = sol
        ax.plot(x, python, color=color, linewidth=1.0,
                linestyle='solid', label=labels[0])
        if with_test:
            ax.plot(x, fortran, color=color, linewidth=1.0,
                    linestyle='dotted', label=labels[1])
            title = '{} with max error rate {}% for {}, {}'.format(
                    substance, err[0],
                    round(python[err[1]], 3),
                    round(fortran[err[1]], 3))
        else:
            title = '{}'.format(substance)
        ax.set_title(title)
        ax.legend()
        if not is_in_linear_scale:
            ax.set_xticks(xticks)
            ax.set_xticklabels(xlabels, fontsize=12)
        ax.set_xlabel('Distance to grain boundary (nm)')
        ax.set_ylabel('Fractional concentration')

    def plot_result(self, is_in_linear_scale=True, with_test=False,
                    is_save=False, is_only_Cr=False):
        """
        To visiualize the contrations of Fe, Cr, Ni

        By default, we plot in a linear scale without test data
        """
        x = range(self.NSTEP)
        idx = list(x[::10]) + [x[-1]]
        xticks = np.array(x)[idx]
        xlabels = np.around(self.XVALUE*1e9, 2).astype(str)[idx]
        if is_in_linear_scale:
            x = self.XVALUE*1e9

        Fe_python, Cr_python, Ni_python = self.CA, self.CB, self.CC

        Fe_fortran = Cr_fortran = Ni_fortran = None
        Fe_err_idx = Cr_err_idx = Ni_err_idx = None
        Fe_error = Cr_error = Ni_error = None
        if with_test:
            Y = np.loadtxt("test/DOSE/{}.txt".format(self.DOSES[-1]))
            Fe_fortran, Cr_fortran, Ni_fortran = Y[:, 1], Y[:, 2], Y[:, 3]
            Fe_err_idx = abs(Fe_fortran-Fe_python).argmax()
            Cr_err_idx = abs(Cr_fortran-Cr_python).argmax()
            Ni_err_idx = abs(Ni_fortran-Ni_python).argmax()
            Fe_error = round(max(abs(Fe_fortran-Fe_python)/Fe_fortran)*100, 4)
            Cr_error = round(max(abs(Cr_fortran-Cr_python)/Cr_fortran)*100, 4)
            Ni_error = round(max(abs(Ni_fortran-Ni_python)/Ni_fortran)*100, 4)
            print('##########################################################')
            print('Compare the result between Python and Fortran code:')
            print('When DOSE = {}'.format(self.DOSES[-1]))
            print('Max error rate for Fe = {}%'.format(Fe_error))
            print('Max error rate for Cr = {}%'.format(Cr_error))
            print('Max error rate for Ni = {}%'.format(Ni_error))
            print('##########################################################')
        if is_only_Cr:
            _, axs = plt.subplots(1, 1, figsize=(8, 5))
        else:
            _, axs = plt.subplots(1, 3, figsize=(20, 5))
        substances = ['Fe', 'Cr', 'Ni']
        errs = [[Fe_error, Fe_err_idx],
                [Cr_error, Cr_err_idx],
                [Ni_error, Ni_err_idx]]
        sols = [[Fe_python, Fe_fortran],
                [Cr_python, Cr_fortran],
                [Ni_python, Ni_fortran]]
        colors = ['red', 'green', 'blue']

        if is_only_Cr:
            self.plot(axs, x, substances[1], errs[1], sols[1], colors[1],
                      with_test, is_in_linear_scale, xticks, xlabels)
        else:
            for ax, substance, err, sol, color in zip(axs, substances, errs,
                                                      sols, colors):
                self.plot(ax, x, substance, err, sol, color, with_test,
                          is_in_linear_scale, xticks, xlabels)
        if is_save:
            if not os.path.exists('images/'):
                os.makedirs('images/')
            plt.savefig('images/RIS.jpg')
        plt.show()

    def plot_EB_EC_DB_DC(self, token='EB'):
        '''
        Plot evolution of migration energy and diffusion coeffcient with
        time ans space
        '''
        if token =='EB':
            Ys = np.array(self.EBs)
        if token == 'EC':
            Ys = np.array(self.ECs)
        if token == 'DB':
            Ys = np.array(self.DBs)
        if token == 'DC':
            Ys = np.array(self.DCs)
        xlabels = np.around(self.XVALUE*1e9, 2).astype(str)
        _, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
        dists = range(0, Ys.shape[1], 20)
        for i in dists:
            ax1.plot(range(Ys.shape[0]), Ys[:, i],
                     label='Distance to GB = {} nm'.format(xlabels[i]))
        ax1.legend()
        ax1.set_title('Evolution of {} with No. of iterations'.format(token))
        ax1.set_xlabel('No. of iterations')
        if token[0] == 'D': ax1.set_ylabel('Diffusion coefficient')
        else: ax1.set_ylabel('Migration energy (eV)')
        x = range(self.NSTEP)
        idx = list(x[::10]) + [x[-1]]
        if token[0] == 'D':
            ax2.plot(x[:-1], Ys[0, :-1], color='violet',
                     label='Diffusion coefficient')
        else:
            ax2.plot(x[:-1], Ys[0, :-1], color='violet',
                     label='Migration energy')
        ax2.legend()
        ax2.set_xticks(np.array(x)[idx])
        ax2.set_xticklabels(xlabels[idx], fontsize=12)
        ax2.set_title('Evolution of {} with distance to GB'.format(token))
        ax2.set_xlabel('Distance to GB (nm)')
        if token[0] == 'D': ax2.set_ylabel('Diffusion coefficient')
        else: ax2.set_ylabel('Migration energy (eV)')
        plt.savefig('images/{}.jpg'.format(token))
        plt.show()

    def plot_E_D_1st_iter(self):
        '''
        '''
        self.fex([], self.Y0)
        EB = self.EBs[0]
        EC = self.ECs[0]
        DB = self.DBs[0]
        DC = self.DCs[0]
        _, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
        x = range(self.NSTEP)
        idx = list(x[::10]) + [x[-1]]
        xlabels = np.around(self.XVALUE*1e9, 2).astype(str)
        ax1.plot(x[:-1], EB[:-1], label='Cr')
        ax1.plot(x[:-1], EC[:-1], label='Ni')
        ax1.legend()
        ax1.set_xticks(np.array(x)[idx])
        ax1.set_xticklabels(xlabels[idx], fontsize=12)
        ax1.set_title('Migration energy (eV) at 1st iteration')
        ax1.set_xlabel('Distance to GB (nm)')
        ax2.plot(x[:-1], DB[:-1], label='Cr')
        ax2.plot(x[:-1], DC[:-1], label='Ni')
        ax2.legend()
        ax2.set_xticks(np.array(x)[idx])
        ax2.set_xticklabels(xlabels[idx], fontsize=12)
        ax2.set_title('Diffusion coefficient (m2/s) at 1st iteration')
        ax2.set_xlabel('Distance to GB (nm)')
        plt.savefig('images/E-D.jpg')
        plt.show()

    def plot_CV_CI(self, is_in_linear_scale=True):
        '''
        '''
        x = range(self.NSTEP)
        idx = list(x[::10]) + [x[-1]]
        xticks = np.array(x)[idx]
        xlabels = np.around(self.XVALUE*1e9, 2).astype(str)[idx]
        if is_in_linear_scale:
            x = self.XVALUE*1e9
        _, (ax1) = plt.subplots(1, 1)
        self.plot(ax=ax1, x=x, substance='Vacancy', err=None,
                  sol=[self.CV, []], color=None, with_test=False,
                  is_in_linear_scale=is_in_linear_scale,
                  xticks=xticks, xlabels=xlabels,
                  labels=['Vacancy', None])
        self.plot(ax=ax1, x=x, substance='Interstitial', err=None,
                  sol=[self.CI, []], color=None, with_test=False,
                  is_in_linear_scale=is_in_linear_scale,
                  xticks=xticks, xlabels=xlabels,
                  labels=['Interstitial', None])
        plt.title('Concentration of vacancy and interstitial')
        plt.savefig('images/CV_CI.jpg')
        plt.show()

    def compare_lc(self, sols):
        x = range(self.NSTEP)
        idx = list(x[::10]) + [x[-1]]
        xticks = np.array(x)[idx]
        xlabels = np.around(self.XVALUE*1e9, 2).astype(str)[idx]
        x = self.XVALUE*1e9
        _, (ax1) = plt.subplots(1, 1)
        self.plot(ax=ax1, x=x, substance='Cr', err=None,
                  sol=[sols[0], []], color=None, with_test=False,
                  is_in_linear_scale=True,
                  xticks=xticks, xlabels=xlabels,
                  labels=['Without sink strength at GB', None])
        self.plot(ax=ax1, x=x, substance='Cr', err=None,
                  sol=[sols[1], []], color=None, with_test=False,
                  is_in_linear_scale=True,
                  xticks=xticks, xlabels=xlabels,
                  labels=['With sink strength at GB', None])
        plt.title('Cr distribution')
        plt.savefig('images/Comparaison of limit conditions.jpg')
        plt.show()


def get_value_from_profile(x, profile):
    """
    This function is to acquire a corresponded value according an input in 
    a given profile with interpolation method, because the profile is discret
    and it's possible that the profle cannot cover x.

    Input:
        x: a value of distance
        profile: a curve of some a subtance, usually Cr
    Output:
        a value in the profile corresponds to x
    """
    dist = profile[:, 0]
    val = profile[:, 1]
    if x < dist[0]:
        k = (val[1] - val[0])/(dist[1] - dist[0])
        return val[0]
        return val[0] + k*(x - dist[0])
    elif x > dist[-1]:
        return val[-1]
    else:
        i = np.nonzero(dist<=x)[0][-1]
        j = np.nonzero(dist>=x)[0][0]
        if i == j:
            return val[j]
        else:
            k = (val[i] - val[j])/(dist[i] - dist[j])
            return val[j]+k*(x-dist[j])
 


if __name__ == "__main__":
    mik = MIK(DOSES=np.linspace(0, 0.994, 100), ECOHB=-4.21)
    mik.plot_E_D_1st_iter()

    sol = mik.solve(print_time_used=True)
    print(mik.CB[0])
    mik.plot_result(with_test=True, is_in_linear_scale=False, is_save=True,
                    is_only_Cr=False)
    mik.convolution(is_plot=True)
    mik.plot_CV_CI()
    '''
    mik.plot_EB_EC_DB_DC('EB')
    mik.plot_EB_EC_DB_DC('EC')
    mik.plot_EB_EC_DB_DC('DB')
    mik.plot_EB_EC_DB_DC('DC')
    '''
