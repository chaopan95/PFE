#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 16:56:04 2020

@author: panchao
"""


import os
import unittest
import numpy as np
from RIS_MIK import MIK


class MIK_Test(unittest.TestCase):
    """
    test for figure
    test for comparing solution with ancient one
    current solution with reference
    test for each function
    blobal test for physique
    unit test for numerique programmination
    """
    def test_initial_values(self):
        '''
        Test all initial values of MIK model
        '''
        mik = MIK()

        EPS=0.1000E-08
        self.assertEqual(EPS, mik.EPS)

        DISPRT=0.140000000E-05
        self.assertEqual(DISPRT, mik.DISPRT)

        ETAV=1.00000000
        self.assertEqual(ETAV, mik.ETAV)

        ETAI=1.00000000
        self.assertEqual(ETAI, mik.ETAI)

        TEMP=360.00000000
        self.assertEqual(TEMP, mik.TEMPC)

        CB=0.2100000
        self.assertEqual(CB, mik.CONCB)

        CC=0.0900000
        self.assertEqual(CC, mik.CONCC)

        DISL=0.10000000E+15
        self.assertEqual(DISL, mik.DISL)

        NAT=0.91000000E+29
        self.assertEqual(NAT, mik.NAT)

        LAMBDA=0.35000000E-09
        self.assertEqual(LAMBDA, mik.LAMBDA)

        FAV=0.785
        self.assertEqual(FAV, mik.FAV)

        FBV=0.668
        self.assertEqual(FBV, mik.FBV)

        FCV=0.872
        self.assertEqual(FCV, mik.FCV)

        FI=0.660
        self.assertEqual(FI, mik.FI)

        WAV=1.80000000
        self.assertEqual(WAV, mik.WAV)

        WBV=3.20000000
        self.assertEqual(WBV, mik.WBV)

        WCV=1.00000000
        self.assertEqual(WCV, mik.WCV)

        WAI=1.00000000
        self.assertEqual(WAI, mik.WAI)

        WBI=1.00000000
        self.assertEqual(WBI, mik.WBI)

        WCI=1.00000000
        self.assertEqual(WCI, mik.WCI)

        ECOHA=-4.28000000
        self.assertEqual(ECOHA, mik.ECOHA)

        ECOHB=-4.21000000
        self.assertEqual(ECOHB, mik.ECOHB)

        ECOHC=-4.44000000
        self.assertEqual(ECOHC, mik.ECOHC)

        EMIA=0.90000000
        self.assertEqual(EMIA, mik.EMIA)

        EMIB=0.90000000
        self.assertEqual(EMIB, mik.EMIB)

        EMIC=0.90000000
        self.assertEqual(EMIC, mik.EMIC)

        SV=1.00000000
        self.assertEqual(SV, mik.SV)

        EMA=1.28000000
        self.assertEqual(EMA, mik.EMA)

        EMB=0.97000000
        self.assertEqual(EMB, mik.EMB)

        EMC=1.04000000
        self.assertEqual(EMC, mik.EMC)

        EFA=1.40000000
        self.assertEqual(EFA, mik.EFA)

        EFB=1.60000000
        self.assertEqual(EFB, mik.EFB)

        EFC=1.79000000
        self.assertEqual(EFC, mik.EFC)

        EFGB=1.40000000
        self.assertEqual(EFGB, mik.EFGB)

        EORDAB=0.00300000
        self.assertEqual(EORDAB, mik.EORDAB)

        EORDAC=-0.00100000
        self.assertEqual(EORDAC, mik.EORDAC)

        EORDBC=0.00500000
        self.assertEqual(EORDBC, mik.EORDBC)

        NUOV=0.1500E+14
        self.assertEqual(NUOV, mik.NUOV)

        NUOI=0.1500E+13
        self.assertEqual(NUOI, mik.NUOI)

        AL=1.00000000
        self.assertEqual(AL, mik.AL)

        Z=12.00
        self.assertEqual(Z, mik.Z)

        BIASV=1.00
        self.assertEqual(BIASV, mik.BIASV)

        BIASI=1.00
        self.assertEqual(BIASI, mik.BIASI)

        TFRAC=[1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000]
        self.assertEqual(TFRAC, mik.TFRAC[:-1].tolist())

        CAFRAC=[1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000]
        self.assertEqual(CAFRAC, mik.CAFRAC.tolist())

        CBFRAC=[1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000]
        self.assertEqual(CBFRAC, mik.CBFRAC.tolist())

        CCFRAC=[1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 1.0000]

        self.assertEqual(CCFRAC, mik.CCFRAC.tolist())
        DFRAC=[1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000]
        self.assertEqual(DFRAC, mik.DFRAC.tolist())

        SFRAC=[1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
               1.0000, 1.0000]
        self.assertEqual(SFRAC, mik.SFRAC.tolist())
        
        PREVA = 21195000000000.000
        self.assertEqual(PREVA, mik.PREVA)
        
        NUIA = 67576.686846171171
        self.assertEqual(NUIA, mik.NUIA[0])
        
        TKT = 5.4545611747016665E-002
        self.assertEqual(TKT, mik.TKT[0])

    def test_fex(self):
        '''
        After one temporal iteration, we have an increment for Fe, Cr, Ni,
        vancy and interstitial
        '''
        Y = np.loadtxt('test/perks.Y')[:, 1:]

        NA = np.loadtxt('test/perks.NA')[:, 1:]
        NB = np.loadtxt('test/perks.NB')[:, 1:]
        NC = np.loadtxt('test/perks.NC')[:, 1:]
        NV = np.loadtxt('test/perks.NV')[:, 1:]
        NI = np.loadtxt('test/perks.NI')[:, 1:]
        EA = np.loadtxt('test/perks.EA')[:, 1:]
        
        NUVA = np.loadtxt('test/perks.NUVA')[:, 1:]
        NUIA = np.loadtxt('test/perks.NUIA')[:, 1:]

        DAV = np.loadtxt('test/perks.DAV')[:, 1:]
        DBV = np.loadtxt('test/perks.DBV')[:, 1:]
        DCV = np.loadtxt('test/perks.DCV')[:, 1:]

        RECA = np.loadtxt('test/perks.RECA')[:, 1:]

        DIFV = np.loadtxt('test/perks.DIFV')[:, 1:]
        DIFI = np.loadtxt('test/perks.DIFI')[:, 1:]

        GRADCA = np.loadtxt('test/perks.GRADCA')[:, 1:]

        DA = np.loadtxt('test/perks.DA')[:, 1:]

        JA = np.loadtxt('test/perks.JA')[:, 1:]

        DIVJA = np.loadtxt('test/perks.DIVJA')[:, 1:]

        YDOT = np.loadtxt('test/perks.YDOT')[:, 1:]

        mik = MIK()
        for i in range(1):  # Y.shape[0]
            mik.fex([], Y[i, :])
            self.assertEqual(NA[i, :].tolist(), mik.NA.tolist())
            self.assertEqual(NB[i, :].tolist(), mik.NB.tolist())
            self.assertEqual(NC[i, :].tolist(), mik.NC.tolist())
            self.assertEqual(NV[i, :].tolist(), mik.NV.tolist())
            self.assertEqual(NI[i, :].tolist(), mik.NI.tolist())
            self.assertEqual(EA[i, :].tolist(), mik.EA.tolist())
            self.assertEqual(NUVA[i, :].tolist(), mik.NUVA.tolist())
            self.assertEqual(NUIA[i, :].tolist(), mik.NUIA.tolist())
            self.assertEqual(DAV[i, :].tolist(), mik.DAV.tolist())
            self.assertEqual(DBV[i, :].tolist(), mik.DBV.tolist())
            self.assertEqual(DCV[i, :].tolist(), mik.DCV.tolist())
            self.assertEqual(RECA[i, :].tolist(), mik.RECA.tolist())
            self.assertEqual(DIFV[i, :].tolist(), mik.DIFV.tolist())
            
            np.allclose(DIFI[i, :], mik.DIFI)
            
            self.assertEqual(GRADCA[i, :].tolist(), mik.GRADCA.tolist())
            self.assertEqual(DA[i, :].tolist(), mik.DA.tolist())
            self.assertEqual(JA[i, :].tolist(), mik.JA.tolist())
            self.assertEqual(DIVJA[i, :].tolist(), mik.DIVJA.tolist())
            self.assertEqual(YDOT[i, :].tolist(), mik.YDOT.tolist())

    def test_DAV_RECA_with_NUVA(self):
        Y = np.loadtxt('test/perks.Y')[:, 1:]
        NUVA = np.loadtxt('test/perks.NUVA')[:, 1:]
        NUIA = np.loadtxt('test/perks.NUIA')[:, 1:]
        RECA = np.loadtxt('test/perks.RECA')[:, 1:]
        mik = MIK()
        for i in range(Y.shape[0]):
            mik.fex([], Y[i, :])
            reca = np.zeros(mik.NSTEP)
            reca[:-1] = (NUVA[i, :-1]+NUIA[i, :-1])*mik.Z
            reca[-1] = reca[-2]
            for j in range(1, mik.NSTEP-1):
                reca[j] = 0.5*(reca[j]+reca[j-1])
            self.assertEqual(RECA[i, :].tolist(), reca.tolist())

    def test_NUVA(self):
        '''
        '''
        Y = np.loadtxt('test/perks.Y')[:, 1:]
        NUVA = np.loadtxt('test/perks.NUVA')[:, 1:]
        mik = MIK()
        for i in range(Y.shape[0]):
            mik.fex([], Y[i, :])
            np.allclose(NUVA[i, :], mik.NUVA)

    def test_solution(self):
        '''
        test final result at different DOSE
        '''
        DOSES = [float(file.split('.txt')[0])
                 for file in os.listdir('test/DOSE/')]
        mik = MIK(DOSES=sorted(DOSES))
        sol = mik.solve()
        y = sol.y
        Ys = np.zeros(y.shape)
        for i in range(len(DOSES)):
            Y = np.loadtxt("test/DOSE/{}.txt".format(DOSES[i]))
            Ys[:, i] = np.reshape(Y[:, 1:], (250,))
        np.allclose(y, Ys)
        mik.plot_result(with_test=True)



if __name__ == '__main__':
    unittest.main()
