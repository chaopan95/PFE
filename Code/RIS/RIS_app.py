#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 16:41:30 2020

@author: panchao
"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from RIS_MIK import MIK


def read_from_csv(filename):
    """
    """
    with open(filename) as f:
        contents = [line.replace(',', '.') for line in f.readlines()]
    with open(filename, 'w') as f:
        f.writelines(contents)
    df = pd.read_csv(filename, sep=';')
    values = df.values
    datasets = {}
    for i in range(df.shape[0]):
        datasets[i] = values[i, 1:]
    return datasets

def plot(Fe_Fortran, Cr_Fortran, Ni_Fortran, Fe_Python, Cr_Python, Ni_Python):
    """
    """
    metric_Fe = np.linspace(min(Fe_Fortran.min(), Fe_Python.min())*0.9,
                            max(Fe_Fortran.max(), Fe_Python.max())*1.1,
                            10)*100
    metric_Cr = np.linspace(min(Cr_Fortran.min(), Cr_Python.min())*0.9,
                            max(Cr_Fortran.max(), Cr_Python.max())*1.1,
                            10)*100
    metric_Ni = np.linspace(min(Ni_Fortran.min(), Ni_Python.min())*0.9,
                            max(Ni_Fortran.max(), Ni_Python.max())*1.1,
                            10)*100
    metrics = [metric_Fe, metric_Cr, metric_Ni]

    substances = ['Fe', 'Cr', 'Ni']
    sols = [[Fe_Fortran*100, Fe_Python*100],
            [Cr_Fortran*100, Cr_Python*100],
            [Ni_Fortran*100, Ni_Python*100]]

    fig, axs  = plt.subplots(1, 3, figsize=(20, 5))
    
    for metric, ax, substance, sol in zip(metrics, axs, substances, sols):
        ax.plot(metric, metric, color='black', linewidth=1, linestyle='solid')
        Fortran, Python = sol
        ax.scatter(Fortran, Python, s=10, color='red')
        ax.grid()
        ax.set_title('Python vs Fortran')
        ax.set_xlabel('{}(%) prévu par Fortran'.format(substance))
        ax.set_ylabel('{}(%) prévu par Python'.format(substance))

    plt.show()

def python_vs_fortran():
    """
    """
    datasets = read_from_csv('analyse/data_set.csv')
    for key in datasets.keys():
        CONCA, CONCB, CONCC, R1, R2, RF, N1, N2, N3, ECOHA, ECOHB, ECOHC,\
        EORDAB, EORDAC, EORDBC, EFA, EFB, EFC, EMA, EMB, EMC,\
        EMIA, EMIB, EMIC, WAV, WBV, WCV, FAV, FBV, FCV, FI, WAI, WBI, WCI,\
        SV, EFGB, TEMPC, DOSE = datasets[key]
        DOSE = round(DOSE, 5)
        ris = MIK(CONCA, CONCB, CONCC,
                  R1, R2, RF,
                  int(N1), int(N2), int(N3),
                  ECOHA, ECOHB, ECOHC,
                  EORDAB, EORDAC, EORDBC,
                  EFA, EFB, EFC,
                  EMA, EMB, EMC,
                  EMIA, EMIB, EMIC,
                  WAV, WBV, WCV,
                  FAV, FBV, FCV, FI,
                  WAI, WBI, WCI,
                  SV,
                  EFGB,
                  TEMPC,
                  DOSE)
        ris.solve()
        Fe_Python, Cr_Python, Ni_Python = ris.CA, ris.CB, ris.CC
        Y = np.loadtxt("test/DOSE/{}.txt".format(DOSE))
        Fe_Fortran, Cr_Fortran, Ni_Fortran = Y[:, 1], Y[:, 2], Y[:, 3]
        plot(Fe_Fortran, Cr_Fortran, Ni_Fortran,
             Fe_Python, Cr_Python, Ni_Python)

def python_vs_observation():
    """
    """

def comparaison_python_fortran():
    doses = []
    for txt in os.listdir('test/DOSE/'):
        if 'txt' in txt:
            dose = txt.split('.txt')[0]
            doses.append(float(dose))

    idx_distance = [0, 1, 4, 22, 30, 39]
    columns = ['parameters_set', 'temperature', 'dose', 'Cr_initial',
               'Cr_0nm_py', 'Cr_0.25nm_py', 'Cr_1nm_py', 'Cr_10.75nm_py',
               'Cr_117.75nm_py', 'Cr_1017.75nm_py',
               'Cr_0nm_f', 'Cr_0.25nm_f', 'Cr_1nm_f', 'Cr_10.75nm_f',
               'Cr_117.75nm_f', 'Cr_1017.75nm_f']
    df = pd.DataFrame(columns=columns)
    idx = 0

    for dose in sorted(doses):
        print('****************DOSE={} Start****************'.format(dose))
        ris = MIK(DOSE=dose)
        ris.solve()
        Fe_Python, Cr_Python, Ni_Python = ris.CA, ris.CB, ris.CC
        Y = np.loadtxt("test/DOSE/{}.txt".format(dose))
        Fe_Fortran, Cr_Fortran, Ni_Fortran = Y[:, 1], Y[:, 2], Y[:, 3]
        plot(Fe_Fortran, Cr_Fortran, Ni_Fortran,
             Fe_Python, Cr_Python, Ni_Python)
        print('****************DOSE={} End****************'.format(dose))
        print()

        df.loc[idx, 'parameters_set'] = idx
        df.loc[idx, 'temperature'] = ris.TEMPC
        df.loc[idx, 'dose'] = dose
        df.loc[idx, 'Cr_initial'] = ris.CONCB
        df.loc[idx, ['Cr_0nm_py', 'Cr_0.25nm_py', 'Cr_1nm_py',
                     'Cr_10.75nm_py', 'Cr_117.75nm_py',
                     'Cr_1017.75nm_py']] = ris.CB[idx_distance]
        df.loc[idx, ['Cr_0nm_f', 'Cr_0.25nm_f', 'Cr_1nm_f',
                     'Cr_10.75nm_f', 'Cr_117.75nm_f',
                     'Cr_1017.75nm_f']] = Cr_Fortran[idx_distance]
        idx += 1
    df.to_csv('analyse/python_vs_fortran_Cr.csv', sep=';', index=False)

def plot_python_vs_fortran_Cr():
    """
    '0nm', '0.25nm', '1nm', '10.75nm', '117.75nm', '1017.75nm'
    """
    df = pd.read_csv('analyse/python_vs_fortran_Cr.csv', sep=';')
    for idx in range(df.shape[0]):
        dose = round(df.loc[idx, 'dose'], 5)
        Cr_Fortran = df.loc[idx, ['Cr_0nm_f', 'Cr_0.25nm_f', 'Cr_1nm_f',
                                  'Cr_10.75nm_f', 'Cr_117.75nm_f',
                                  'Cr_1017.75nm_f']].values
        Cr_Python = df.loc[idx, ['Cr_0nm_py', 'Cr_0.25nm_py', 'Cr_1nm_py',
                                 'Cr_10.75nm_py', 'Cr_117.75nm_py',
                                 'Cr_1017.75nm_py']].values
        metric_Cr = np.linspace(min(Cr_Fortran.min(), Cr_Python.min())*0.9,
                                max(Cr_Fortran.max(), Cr_Python.max())*1.1,
                                10)*100
        plt.figure()
        plt.plot(metric_Cr, metric_Cr, color='black',
                 linewidth=1, linestyle='solid')
        plt.scatter(Cr_Fortran*100, Cr_Python*100, s=10, color='red')
        plt.grid()
        plt.title('Python vs Fortran quand DOSE = {}'.format(dose))
        plt.xlabel('Cr(%) prévu par Fortran')
        plt.ylabel('Cr(%) prévu par Python')
        plt.show()

def plot_in_linear_scale():
    ris = MIK()
    ris.solve()
    ris.plot(is_in_linear_scale=True)

def give_Cr_profile(DOSES=[0.1, 0.5, 1, 2], file_name='TNES_304.csv'):
    """
    This function takes a deterministic profile of Cr and set Fe = 70%,
    Ni = 30% - Cr as initial values.
    Besides, this function receives a group of dose in order to represent a
    final diagram of each substance under these doses.
    """
    Cr_profile = pd.read_csv(file_name, sep=';').values
    Cr_profile_right = Cr_profile[np.nonzero(Cr_profile[:, 0]>=0)]
    Cr_profile_left = Cr_profile[np.nonzero(Cr_profile[:, 0]<=0)]
    Cr_profile_left[:, 0] = -Cr_profile_left[::-1, 0]
    Cr_profile_left[:, 1] = Cr_profile_left[::-1, 1]

    ris_right = MIK(RF=58, DOSES=DOSES, Cr_profile=Cr_profile_right)
    ris_left = MIK(RF=58, DOSES=DOSES, Cr_profile=Cr_profile_left)

    xvalue = np.concatenate((-ris_left.XVALUE[::-1],
                             ris_right.XVALUE[1:]), axis=0)*1e9

    sol_r = ris_right.solve()
    sol_l = ris_left.solve()
    NSTEP = ris_left.NSTEP

    CA_l = sol_l.y[:NSTEP, :].T
    CA_r = sol_r.y[:NSTEP, :].T
    CA = np.concatenate((CA_l[:, ::-1], CA_r[:, 1:]), axis=1)

    CB_l = sol_l.y[NSTEP:2*NSTEP, :].T
    CB_r = sol_r.y[NSTEP:2*NSTEP, :].T
    CB = np.concatenate((CB_l[:, ::-1], CB_r[:, 1:]), axis=1)

    CC_l = sol_l.y[2*NSTEP:3*NSTEP, :].T
    CC_r = sol_r.y[2*NSTEP:3*NSTEP, :].T
    CC = np.concatenate((CC_l[:, ::-1], CC_r[:, 1:]), axis=1)

    fig, axs  = plt.subplots(1, 3, figsize=(20, 5))
    substances = ['Fe', 'Cr', 'Ni']
    sols = [CA, CB, CC]
    colors = ['red', 'green', 'blue']
    for ax, substance, sol, color in zip(axs, substances, sols, colors):
        for i in range(len(DOSES)):
            DOSE = DOSES[i]
            ax.plot(xvalue, sol[i], linewidth=1.0, linestyle='solid',
                    label='DOSE={}'.format(DOSE), color=color)
        ax.legend()
        title = substance
        ax.set_title(title)
        ax.set_xlabel('Distance to grain boundary (nm)')
        ax.set_ylabel('Fractional concentration')
    if not os.path.exists('images/'):
        os.makedirs('images/')
    plt.savefig('images/Initial State Effect.jpg')
    plt.show()
    return sols, xvalue

def compare_initial_profil_effect():
    sols1, xvalue = give_Cr_profile(DOSES=[2], file_name='Default.csv')
    sols2, xvalue = give_Cr_profile(DOSES=[2], file_name='TNES_304.csv')
    _, axs = plt.subplots(1, 3, figsize=(20, 5))
    substances = ['Fe', 'Cr', 'Ni']
    colors = ['red', 'green', 'blue']
    for ax, substance, sol1, sol2, color in zip(axs, substances, sols1,
                                                sols2, colors):
        ax.plot(xvalue, sol1[-1], linewidth=1.0, linestyle='solid',
                label='Default')
        ax.plot(xvalue, sol2[-1], linewidth=1.0, linestyle='solid',
                label='TNES_304')
        ax.legend()
        title = substance
        ax.set_title(title)
        ax.set_xlabel('Distance to grain boundary (nm)')
        ax.set_ylabel('Fractional concentration')
    if not os.path.exists('images/'):
        os.makedirs('images/')
    plt.savefig('images/Initial State Effect.jpg')
    plt.show()

def evolutionCrWithDose():
    fluence = 1
    DOSES = np.linspace(0, fluence, 3000000)
    mik = MIK(DOSES=DOSES, CONCB=0.16, CONCC=0.14, RF=218)
    sol = mik.solve()
    Crs = sol.y[mik.NSTEP, :]
    plt.plot(DOSES, Crs)
    plt.xlabel('DOSE (dpa)')
    plt.ylabel('Fractional concentration of Cr')
    plt.title('Predicted Cr at grain boundary')
    plt.grid()
    plt.show()

def initial_Cr_profil(file_name='TNES_304.csv'):
    Cr_profile = pd.read_csv(file_name, sep=';').values
    _, axs = plt.subplots(1, 3, figsize=(20, 5))
    x = Cr_profile[:, 0]
    Cr = np.array(Cr_profile[:, 1]) / 100
    Ni = 0.3 - Cr
    Fe = 0.7*np.ones(Cr.shape)
    concentrations = [Fe, Cr, Ni]
    colors = ['red', 'green', 'blue']
    substances = ['Fe', 'Cr', 'Ni']
    for ax, color, concentration, substance in zip(axs, colors, concentrations,
                                                   substances):
        ax.plot(x, concentration, color=color, linestyle='solid')
        ax.set_title(substance)
        ax.set_xlabel('Distance to grain boundary (nm)')
        ax.set_ylabel('Fractional concentration')
    plt.savefig('images/Initial Cr profile.jpg')
    plt.show()


if __name__ == "__main__":
    #comparaison_python_fortran()
    #plot_in_linear_scale()
    give_Cr_profile(DOSES=[140], file_name='TNES_304.csv')
    #compare_initial_profil_effect()
    #evolutionCrWithDose()
    initial_Cr_profil()
