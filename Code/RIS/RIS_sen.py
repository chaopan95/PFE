#!/usr/bin/env python
# -*- coding: utf-8 -*-
# coding = utf-8

# https://openturns.github.io/openturns/master/examples/meta_modeling/
# functional_chaos.html

import os
import openturns as ot
import numpy as np
from shutil import rmtree
from scipy.stats import norm
from RIS_MIK import MIK


# FAV, FBV, FCV, FI: Jump correlation factors ???
# DOSE = 80 a verifier
# remplacer ; avec \t
paras_ris_range = {
    "DISPRT": [1e-6, 1e-5],
    "ETAV": [0.8, 1],
    "TEMPC": [350, 370],
    "DISL": [1e12, 1e14],
    "WAV": [0.5, 4],
    "WBV": [0.5, 4],
    "WCV": [0.5, 4],
    "ECOHA": [1.1*-4.28, 0.9*-4.28],
    "ECOHB": [1.1*-4.21, 0.9*-4.21],
    "ECOHC": [1.1*-4.44, 0.9*-4.44],
    "SV": [0.5, 3],
    "EMA": [0.9*1.28, 1.1*1.28],
    "EMB": [0.9*0.97, 1.1*0.97],
    "EMC": [0.9*1.04, 1.1*1.04],
    "EFA": [0.8*1.4, 1.2*1.4],
    "EFB": [0.8*1.6, 1.2*1.6],
    "EFC": [0.8*1.79, 1.2*1.79],
    "EORDAB": [0*0.003, 2*0.003],
    "EORDAC": [2*-0.001, 0*-0.001],
    "EORDBC": [0*0.005, 2*0.005],
    "AL": [0.8, 1],
    "BIASV": [0.9*1.0, 1.1*1.0],
    "BIASI": [0.9*1.0, 1.1*1.0],
    'FAV': [0.5, 0.9],
    'FBV': [0.5, 0.9],
    'FCV': [0.5, 0.9],
    'FWTM': [1.8, 2.2],
    'sigma': [1.5, 2.5]
    }
# vaier le profil du Cr (CONCA = 0.7)
# 2 parametres pour produire un profile comme TNES_304


def GenerateSample(vectX, N, method="MC"):
    """
    """
    if method == "MC":
        X = vectX.getNumericalSample(N)
    elif method == "QMC":
        dim_in = vectX.getDimension()
        # Uniform quasi-random sample over the unit hypercube
        mySobolSeq = ot.SobolSequence(dim_in)
        U = mySobolSeq.generate(N)
        # Isoprobabilistic transform for each marginal
        #
        # Lengthy calculations: compute quantiles using a double loop
        #distrib = vectX.getDistribution()
        #X = empty((N,dim_in))
        #for i in range(dim_in):
          #margin_i = distrib.getMarginal(i)
          #for n in range(N):
    	#u = array(U)[n,i]
    	#Xi = margin_i.computeQuantile(u)
    	#X[n,i] = array(Xi)
        #
        # Alternative: use the inverse isoprobabilistic transform from
        # OpenTURNS object "myDistribution"
        # Xi is normally distributed (N(0,1))
        Xi = norm.ppf(U)
        myDistribution = vectX.getDistribution()
        transfo_inv = myDistribution.getInverseIsoProbabilisticTransformation()
        X = transfo_inv(ot.Sample(Xi))
    return X

def GenerateExperiencePlan(paras_range=paras_ris_range, n_sample=50000):
    """
    """
    for k,v in paras_range.items():
        assert v[0] < v[1], 'ERROR of v[0]>=v[1] for ranges {}'.format(k)

    # le nombre de paramètre doit être proportionné au dimension de ranges
    len_actual_parameters = len(paras_range)
    # Specify the input random vector.
    # instance de classe densite de proba a definir
    myCollection = ot.DistributionCollection(len_actual_parameters)
    
    for i, p_name in enumerate(paras_range):
        distribution = ot.Uniform(paras_range[p_name][0],
                                  paras_range[p_name][1])
        myCollection[i] = ot.Distribution(distribution)
    myDistribution = ot.ComposedDistribution(myCollection)
    vectX = ot.RandomVector(ot.Distribution(myDistribution))
    
    # Sample the input random vector.
    Xsample = GenerateSample(vectX, n_sample, method='QMC')
    xsample = np.array(Xsample)
    return xsample

def GenerateFolder(input_sets, n_folder=200):
    """
    """
    folder = './sensibility/'
    if os.path.exists(folder):
        rmtree(folder)
    os.makedirs(folder)
    n_sample = input_sets.shape[0]
    n_line = int(n_sample/n_folder)

    for idx in range(n_folder):
        sub_folder = folder + '{}/'.format(idx+1)
        os.makedirs(sub_folder)
        np.savetxt(sub_folder + 'input.txt',
                   input_sets[idx*n_line: (idx+1)*n_line, :])

def generate_Cr_profile(sigma=2, height=0.08, baseline=0.21):
    mu = 0
    gaussian = lambda x: 1/(sigma * np.sqrt(2 * np.pi)) *\
                np.exp( - (x - mu)**2 / (2 * sigma**2) )
    x = np.linspace(-20, 20, 101)
    y = gaussian(x)
    y = y * (height/y.max()) + baseline
    #plt.plot(x, y)
    profile = np.zeros((51, 2))
    profile[:, 0] = x[50:]
    profile[:, 1] = y[50:]*100
    return profile

def sensibility_analyse_RIS():
    folder = 'sensibility/'
    for subfolder in os.listdir(folder):
        input_file_path = folder + subfolder + '/input.txt'
        output_file_path = folder + subfolder + '/output.txt'
        paras = np.loadtxt(input_file_path)
        res = np.zeros((paras.shape[0], 2))
        for i in range(paras.shape[0]):
            DISPRT, ETAV, TEMPC, DISL, WAV, WBV, WCV, ECOHA, ECOHB, ECOHC,\
            SV, EMA, EMB, EMC, EFA, EFB, EFC, EORDAB, EORDAC, EORDBC, AL,\
            BIASV, BIASI, FAV, FBV, FCV, FWTM, sigma = paras[i, :]
            Cr_profile = generate_Cr_profile(sigma=sigma)
            mik = MIK(RF=518, DOSES=[80],
                      DISPRT=DISPRT, ETAV=ETAV, TEMPC=TEMPC, DISL=DISL,
                      WAV=WAV, WBV=WBV, WCV=WCV, ECOHA=ECOHA, ECOHB=ECOHB,
                      ECOHC=ECOHC, SV=SV, EMA=EMA, EMB=EMB, EMC=EMC, EFA=EFA,
                      EFB=EFB, EFC=EFC, EORDAB=EORDAB, EORDAC=EORDAC,
                      EORDBC=EORDBC, AL=AL, BIASV=BIASV, BIASI=BIASI,
                      FAV=FAV, FBV=FBV, FCV=FCV, Cr_profile=Cr_profile)
            try:
                mik.solve()
            except:
                print('subfolder: {}, line: {}'.format(subfolder, i))
                continue
            CB, CC = mik.convolution(FWTM=FWTM)
            # Write dozn concentration of Cr and Ni at grain boundary
            res[i, 0] = CB
            res[i, 1] = CC
        np.savetxt(output_file_path, res)
        return


if __name__ == '__main__':
    input_sets = GenerateExperiencePlan(paras_ris_range)
    GenerateFolder(input_sets)
    #sensibility_analyse_RIS()
    print(input_sets.min(axis=0))
    
