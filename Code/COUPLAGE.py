# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:32:32 2020

@author: C23880
"""

import numpy as np
import matplotlib.pyplot as plt
from time import time
from matplotlib.animation import FuncAnimation, PillowWriter
from RIS.RIS_MIK import MIK
from Corrosion.oxydation import Oxydation


def calculate_oxydation(Cr, Ni, T=360, t_final=1e6):
    oxy = Oxydation(x_Fe=1-Cr-Ni, x_Cr=Cr, x_Ni=Ni, temperature_in_K=T)
    t = np.zeros(2)
    y = np.zeros(2)
    J_v0 = np.zeros(2)
    J_MCr = np.zeros(2)
    J_ICr = np.zeros(2)
    J_H = np.zeros(2)
    # Initial oxydation film
    y[0] = 1e-12
    x = np.zeros(int(t_final+1))
    for ti in range(int(t_final+1)):
        x[ti] = y[0]
        # Calculate an increment of thickness with 1 second
        t, y, J_v0, J_MCr, J_ICr, J_H =\
                        oxy.time_integration_RK2_Voyshnis_2014(1, 1, t, y,
                                                               J_v0, J_MCr,
                                                               J_ICr, J_H)
        t[0] = t[1]
        y[0] = y[1]
        J_v0[0] = J_v0[1]
        J_MCr[0] = J_MCr[1]
        J_ICr[0] = J_ICr[1]
        J_H[0] = J_H[1]
    return x

def couple_RIS_OXY(Cr, Ni, T=360, t_final=1e6, DISPRT=1.4e-6):
    '''Predict oxydation thickness with RIS'''
    assert 0 < (Cr + Ni) < 1, 'Initial of Cr and Ni should in [0, 1]'
    t_start = time()
    # Initialise MIK model
    mik = MIK(CONCA=(1-Cr-Ni), CONCB=Cr, CONCC=Ni, DISPRT=DISPRT, DOSES=None)
    # Initialise Oxydation model
    oxy = Oxydation(x_Fe=1-Cr-Ni, x_Cr=Cr, x_Ni=Ni, temperature_in_K=T)
    
    # Get all solutions of MIK
    # In this couple model, Cr is determined onlu by MIK.
    time_interval = np.linspace(0, t_final, int(t_final+1))
    sol = mik.solve(time_interval=time_interval).y
    # Fractional concentrations of Fe, Cr, Ni at grain boundary
    XFe = sol[0, :]
    XCr = sol[mik.NSTEP, :]
    XNi = sol[2*mik.NSTEP, :]
    # Temporary variables at Oxydation model: time, thickness,
    # oxygen vacancies flux, Cr cations flux, Cr interstitial flux
    # flux of all species
    t = np.zeros(2)
    y = np.zeros(2)
    J_v0 = np.zeros(2)
    J_MCr = np.zeros(2)
    J_ICr = np.zeros(2)
    J_H = np.zeros(2)
    # Initial oxydation film
    y[0] = 1e-12
    # All oxydation thickness
    x = np.zeros(int(t_final+1))
    for ti in range(int(t_final+1)):
        x[ti] = y[0]
        # Update concentration of Fe, Cr, Ni dominated by MIK
        oxy.x_Fe = XFe[ti]
        oxy.x_Cr = XCr[ti]
        oxy.x_Ni = XNi[ti]
        # Calculate an increment of thickness with 1 second
        t, y, J_v0, J_MCr, J_ICr, J_H =\
                        oxy.time_integration_RK2_Voyshnis_2014(1, 1, t, y,
                                                               J_v0, J_MCr,
                                                               J_ICr, J_H)
        t[0] = t[1]
        y[0] = y[1]
        J_v0[0] = J_v0[1]
        J_MCr[0] = J_MCr[1]
        J_ICr[0] = J_ICr[1]
        J_H[0] = J_H[1]
    # Time in hour
    h = time_interval / 3600
    t_end = time()
    print(round(t_end-t_start, 2))
    return h, x, XCr

def plot_couplage_result(h, x, XCr, x0):
    fig, ax1 = plt.subplots()
    color='tab:red'
    ax1.set_xlabel('Time (h)')
    ax1.set_ylabel('Thickness (nm)', color=color)
    ax1.plot(h, x*1e9, color=color, label='Film growth with RIS')
    ax1.plot(h, x0*1e9, color='black', label='Film growth without RIS')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.legend(bbox_to_anchor=(0.95, 0.6), loc='upper right')
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('Fractional concentration of Cr', color=color)
    ax2.plot(h, XCr, color=color, label='Cr at GB')
    ax2.legend(bbox_to_anchor=(0.95, 0.1), loc='lower right')
    ax2.tick_params(axis='y', labelcolor=color)
    plt.title('Oxydation with RIS')
    plt.savefig('couplage.jpg')
    plt.show()


def film_growth_animation1(h, x, XCr, x0):
    ''''''
    fig, ax1 = plt.subplots()
    fig.set_tight_layout(True)
    color1='tab:red'
    ax1.set_xlabel('Time (h)')
    ax1.set_ylabel('Thickness (nm)', color=color1)
    ax1.plot(h, x0*1e9, color='white')
    ax1.plot(h[0], x[0]*1e9, color=color1, label='Film growth with RIS')
    ax1.plot(h[0], x0[0]*1e9, color='black', label='Film growth without RIS')
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.legend(bbox_to_anchor=(0.95, 0.6), loc='upper right')

    ax2 = ax1.twinx()
    color2 = 'tab:blue'
    ax2.set_ylabel('Fractional concentration of Cr', color=color2)
    ax2.plot(h, XCr, color='white')
    ax2.plot(h[0], XCr[0], color=color2, label='Cr at GB')
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.legend(bbox_to_anchor=(0.95, 0.1), loc='lower right')

    plt.title('Oxydation with RIS')

    def update(split):
        ''''''
        if split > 0:
            ax1.plot(h[split-1: split+1], x[split-1: split+1]*1e9,
                     color=color1, label='Film growth with RIS')
            ax1.plot(h[split-1: split+1], x0[split-1: split+1]*1e9,
                     color='black', label='Film growth without RIS')
            ax2.plot(h[split-1: split+1], XCr[split-1: split+1],
                     color=color2, label='Cr at GB')
        plt.title('Time = {}h'.format(round(h[split], 2)))
        return ax1, ax2

    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    anim = FuncAnimation(fig, update, frames=np.arange(0, len(h)))
    writergif = PillowWriter(fps=80) 
    anim.save('couplage1.gif', writer=writergif)
    plt.show()

def film_growth_animation2(h, x, XCr, step):
    ''''''
    minVal = XCr[-1]
    maxVal = XCr[0]
    y = np.linspace(minVal, maxVal, 256*4)
    fig, ax = plt.subplots(frameon=False)
    for side in ['left', 'right', 'top', 'bottom']:
        ax.spines[side].set_visible(False)
    ax.plot(1, minVal)
    colors = []
    for i, j in zip(range(255), range(1, 256)):
        color = [0, i/255, 1]
        ax.plot([0, 0], [y[i], y[j]], color=color, linewidth=20)
        colors.append(color)
    for i, j in zip(range(255), range(1, 256)):
        color = [0, 1, 1-i/255]
        ax.plot([0, 0], [y[i+256], y[j+256]], color=color,
                linewidth=20)
        colors.append(color)
    for i, j in zip(range(255), range(1, 256)):
        color = [i/255, 1, 0]
        ax.plot([0, 0], [y[i+256*2], y[j+256*2]], color=color,
                linewidth=20)
        colors.append(color)
    for i, j in zip(range(255), range(1, 256)):
        color = [1, 1-i/255, 0]
        ax.plot([0, 0], [y[i+256*3], y[j+256*3]], color=color,
                linewidth=20)
        colors.append(color)
    colors = colors[::-1]
    ax.get_xaxis().set_visible(False)
    nCr = len(XCr)
    nColor = len(colors)

    xPos1 = 0.4
    xPos2 = 2
    yPos = (minVal+maxVal)/2
    length = int(x[-1]*1e9) + 1
    ax.text(xPos1, yPos*0.95, '0 nm')
    ax.text(xPos2, yPos*0.95, '{} nm'.format(length))
    def update(i):
        e = xPos1+x[i]*1e9/length*(xPos2-xPos1)
        ax.plot([e, xPos2], [yPos, yPos],
                color=colors[int(i*nColor/nCr)],
                linewidth=10)
        ax.plot([xPos1, e], [yPos, yPos], color='black', linewidth=10)
        plt.title('Time = {}h'.format(round(h[i], 2)))
        return ax
    
    anim = FuncAnimation(fig, update, frames=np.arange(0, nCr)[::step])
    writergif = PillowWriter(fps=80) 
    anim.save('couplage2.gif', writer=writergif)
    plt.show()

def run():
    x0 = calculate_oxydation(Cr=0.21, Ni=0.09, t_final=1e4)
    h, x, XCr = couple_RIS_OXY(Cr=0.21, Ni=0.09, t_final=1e4)
    plot_couplage_result(h, x, XCr, x0)
    step = 100
    film_growth_animation1(h[::step], x[::step], XCr[::step], x0[::step])
    film_growth_animation2(h, x, XCr, step)
    print('finished')


if __name__ == '__main__':
    run()
