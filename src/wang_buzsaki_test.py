#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 17:16:25 2019

@author: alafuzof
"""
import brian2.only as br2
from brian2.units import *
from wang_buzsaki import WangBuzsakiCell
import matplotlib.pyplot as plt
import numpy as np

def plot_baseline(simulation_time=3000*ms):
    wb_cell = WangBuzsakiCell()
    wb_group = wb_cell.generate_neuron_group(n_cells=1, dt=0.01*ms)
    wb_group.I_inj = 0.0*uA
    
    net = br2.Network(wb_group)
    
    mon = br2.StateMonitor(wb_group, True, record=True)
    net.add(mon)
    
    net.run(simulation_time)
    
    fig, ax = plt.subplots(1,1)
    ax.plot(mon.t/ms, mon.V[0]/mV, label='V')
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Membrane voltage (mV)')

    return ax

def plot_1uA_injection(simulation_time=100*ms):
    wb_cell = WangBuzsakiCell()
    wb_group = wb_cell.generate_neuron_group(n_cells=1, dt=0.01*ms)
    wb_group.I_inj = 1.0*uA
    
    net = br2.Network(wb_group)
    
    mon = br2.StateMonitor(wb_group, True, record=True)
    net.add(mon)
    
    net.run(simulation_time)
    
    fig, ax = plt.subplots(1,1)
    ax.plot(mon.t/ms, mon.V[0]/mV, label='V')
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Membrane voltage (mV)')

    return ax

def reproduce_fig1a():
    wb_cell = WangBuzsakiCell()
    wb_group = wb_cell.generate_neuron_group(n_cells=1000,  dt=0.01*ms)
    wb_group.I_inj = 'i*(20.0/N)*uA'
    
    net = br2.Network(wb_group)
    
    mon = br2.SpikeMonitor(wb_group, record=False)
    net.add(mon)
    
    net.run(1000*ms)
    
    fig, ax = plt.subplots(1,1)
    x = np.linspace(0,20.0,1000)
    y = mon.count/(1000*ms)
    ax.plot(x,y)
    ax.set_xlabel('Injection current (µA/cm2)')
    ax.set_ylabel('Firing frequency (Hz)')
    
    return ax, mon 

def main():
    import os
    os.makedirs('../figures/', exist_ok=True)
    
    ax = plot_baseline()
    ax.figure.savefig('../figures/test_wang_buzsaki_baseline.png')
    
    ax = plot_1uA_injection()
    ax.figure.savefig('../figures/test_wang_buzsaki_1uA_injection.png')
    
    ax, mon = reproduce_fig1a()
    ax.figure.savefig('../figures/test_wang_buzsaki_fig1a.png')
    

if __name__ == '__main__':
    main()
    