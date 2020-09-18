#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 15:37:19 2020

@author: panchao
"""
import os
import numpy as np
import unittest


def read_result(filename):
    with open(filename, 'r') as f:
        content = [line.replace(', ', '\t') for line in f.readlines()]

    header = content[0]
    with open('temp.txt', 'w') as f:
        f.writelines(content[1:])
    
    data = np.loadtxt('temp.txt')
    os.remove('temp.txt')
    return header, data


class CorrisonTest(unittest.TestCase):
    def test_result_123(self):
        for i in range(1, 4):
            result_file =\
                'Outputs/c_solver_corrosion_evolution_test_{}.output'.format(i)
            reference_file =\
                'refs/c_solver_corrosion_evolution_test_{}.output'.format(i)
            _, res_data = read_result(result_file)
            _, ref_data = read_result(reference_file)
    
            #err = abs(res_data[1:, :] - ref_data[1:, :]) / ref_data[1:, :]
            #print(err.max())
            assert np.allclose(res_data, ref_data), 'fail to pass the test'

    def test_result_4567(self):
        for i in range(4, 8):
            result_file =\
                'Outputs/c_solver_corrosion_evolution_test_{}.output'.format(i)
            reference_file =\
                'refs/c_solver_corrosion_evolution_test_{}.output'.format(i)
            _, res_data = read_result(result_file)
            _, ref_data = read_result(reference_file)

            #err = abs(res_data[1:, :-1] - ref_data[1:, :-1])/ref_data[1:, :-1]
            #print(err.max())
            assert np.allclose(res_data, ref_data), 'fail to pass the test'

    def test_result_8(self):
        result_file = 'Outputs/c_solver_corrosion_evolution_test_8.output'
        reference_file = 'refs/c_solver_corrosion_evolution_test_8.output'
        _, res_data = read_result(result_file)
        _, ref_data = read_result(reference_file)

        err = abs(res_data[1:, :] - ref_data[1:, :]) / ref_data[1:, :]
        print('\nMax error rate = %.2e'%err.max())
        assert np.allclose(res_data, ref_data), 'fail to pass the test'


if __name__ == '__main__':
    unittest.main()