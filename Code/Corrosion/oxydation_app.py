#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 14:51:50 2020

@author: panchao
"""


import os
import numpy as np
import matplotlib.pyplot as plt


def read_output(filename):
    with open(filename, 'r+') as f:
        data = [line.replace(', ', '\t') for line in f.readlines()]
    headers = data[0].strip().split('\t')
    with open('temp.txt', 'w+') as f:
        f.writelines(data[1:])
    data = np.loadtxt('temp.txt')
    os.remove("temp.txt")
    return headers, data

def plot(x, y, modelname, label, xlabel, ylabel, save_file_name):
    plt.plot(x, y, label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.title(modelname)
    if not os.path.exists('images/'):
        os.makedirs('images/')
    plt.savefig('images/{}.jpg'.format(save_file_name))
    plt.show()

def show_Seyeux_2010():
    headers, data = read_output(
        '/Outputs/c_solver_corrosion_evolution_test_1.output')
    plot(data[:, 0]/3600, data[:, 1]*1e9, modelname='Seyeux_2010',
         label='calculated thickness', xlabel='time, t[h]',
         ylabel='thickness, x[nm]',
         save_file_name='Seyeux_2010_calculated thickness')

def show_Leistner_2012():
    headers, data = read_output(
        '/Outputs/c_solver_corrosion_evolution_test_4.output')
    plot(data[:, 0]/3600, data[:, 1]*1e9, modelname='Leistner_2012',
         label='calculated thickness', xlabel='time, t[h]',
         ylabel='thickness, x[nm]',
         save_file_name='Leistner_2012_calculated thickness')
    plot(data[1:, 0]/3600, data[1:, 3]/1e13, modelname='Leistner_2012',
         label='J_v0', xlabel='time, t[h]',
         ylabel='species flux, J[1e13 1/m2/s]',
         save_file_name='Leistner_2012_species flux')

def show_Voyshnis_2014():
    headers, data = read_output(
        'Outputs/c_solver_corrosion_evolution_test_8.output')
    plot(data[:, 0]/3600, data[:, 1]*1e9, modelname='Voyshnis_2014',
         label='calculated thickness', xlabel='time, t[h]',
         ylabel='thickness, x[nm]',
         save_file_name='Voyshnis_2014_calculated thickness')
    plot(data[1:, 0]/3600, data[1:, 3]/1e13, modelname='Voyshnis_2014',
         label='J_v0', xlabel='time, t[h]',
         ylabel='species flux, J[1e13 1/m2/s]',
         save_file_name='Voyshnis_2014_J_v0')
    plot(data[1:, 0]/3600, data[1:, 4]/1e13, modelname='Voyshnis_2014',
         label='J_H', xlabel='time, t[h]',
         ylabel='species flux, J[1e13 1/m2/s]',
         save_file_name='Voyshnis_2014_J_H')

if __name__ == '__main__':
    #show_Seyeux_2010()
    #show_Leistner_2012()
    show_Voyshnis_2014()