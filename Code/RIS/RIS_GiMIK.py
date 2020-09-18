#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 10:58:38 2020

@author: panchao
"""


from math import sin, pi
from RIS_MIK import MIK


class GiMIK(MIK):
    '''
    Grain boundary-Interstitial-Modified Inverse Kirkendall (GiMIK) is used to
    describe the model which incorporates grain boundary structure and
    in terstitial binding effect.
    '''
    def __init__(self,
                 limit_condition=2,
                 theta=3.5,
                 R1=1,
                 R2=2,
                 RF=10,
                 DOSES=[80],
                 ECOHB=-4.21):
        '''
        Input parameters to modify boundary conditions on the GiMIK moedel
        theta: tilt angle (1° = PI/180)
        '''
        super(GiMIK, self).__init__(R1=R1,
                                    R2=R2,
                                    RF=RF,
                                    DOSES=DOSES,
                                    ECOHB=ECOHB)
        # 1: without sink strength
        # 2: with sink strength everywhere
        # 3: with sink strength everywhere + sink strength at GB
        self.limit_condition = limit_condition
        self.theta = theta
        # a parameter in unit if the inverse sauqre of the length (/m2)
        self.A = 1e20
        # sink strength
        self.SGB = self.A * sin(self.theta*(pi/180)/2)
        # Internal sink strength for vacancy
        # eletron irradiation 1e14 /m2; neutron irradiation 1e15 /m2
        self.Sv = 0#1e17
        # Internal sink strength for interstitial
        self.Si = 1.1*self.Sv

    def calculate_dCv_dCi(self):
        '''
        increment of fractional concentration for vancancy and interstitial
        '''
        NAT = self.NAT
        BIASV = self.BIASV
        BIASI = self.BIASI
        self.CVDOT[0] = 0.0
        self.CIDOT[0] = 0.0
        Gv = self.DISPV[1:] - BIASV*self.VACSINK[1:]*self.CV[1:] +\
            self.VACSOUR[1:]
        Gi = self.DISPI[1:] - BIASI*self.INTSINK[1:]*self.CI[1:]
        self.CVDOT[1:] = -self.DIVJV[1:]/NAT + Gv -\
                         self.RECOMB[1:]*self.CV[1:]*self.CI[1:]
        self.CIDOT[1:] = -self.DIVJI[1:]/NAT + Gi -\
                         self.RECOMB[1:]*self.CV[1:]*self.CI[1:]
        if self.limit_condition == 1:
            self.CVDOT[1:] = -self.DIVJV[1:]/NAT + self.DISPV[1:] -\
                             self.RECOMB[1:]*self.CV[1:]*self.CI[1:]
            self.CIDOT[1:] = -self.DIVJI[1:]/NAT + self.DISPI[1:] -\
                             self.RECOMB[1:]*self.CV[1:]*self.CI[1:]
        if self.limit_condition == 3:
            self.CVDOT[1] -= self.SGB*self.DIFV[1]*(self.CV[1]-self.CVTHER[1])
            self.CIDOT[1] -= self.SGB*self.DIFI[1]*self.CI[1]



if __name__ == '__main__':
    '''
    gimik = GiMIK(theta=60, R1=4, R2=18, RF=2018, DOSES=[140.0],
                  limit_condition=3)
    gimik.solve()
    gimik.plot_result(is_save=True, with_test=True)

    gimik.convolution(is_plot=True)
    gimik.plot_DCr()
    '''
    sols = []
    for lc in [2, 3]:
        gimik = GiMIK(theta=3.5, R1=4, R2=18, RF=218, DOSES=[140],
                      limit_condition=lc)
        sol = gimik.solve()
        print(gimik.convolution(is_plot=False))
        sols.append(sol.y[gimik.NSTEP:2*gimik.NSTEP, 0])
    gimik.compare_lc(sols)