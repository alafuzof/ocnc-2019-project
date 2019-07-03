#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 16:14:19 2019

@author: alafuzof
"""
import brian2.only as br2
from brian2.units import *
from pinsky_rinzel import PinskyRinzelCell
import matplotlib.pyplot as plt

def simulate_baseline(time=10*second):
    prc = PinskyRinzelCell()
    prc_group = prc.generate_neuron_group(n_cells=1)
    net = br2.Network(prc_group)
    mon = br2.StateMonitor(prc_group, True, record=True)
    net.add(mon)
    net.run(time, report='text')
    
    return prc, prc_group, net, mon

def simulate_fig2a(time=2*second, dt=None):
    prc = PinskyRinzelCell()
    prc.somatic_parameters['Is_soma'] = 0.75*uA/cm2
    prc_group = prc.generate_neuron_group(n_cells=1, dt=None)
    net = br2.Network(prc_group)
    mon = br2.StateMonitor(prc_group, True, record=True)
    net.add(mon)
    net.run(time, report='text')
    
    return prc, prc_group, net, mon

def plot_fig2_style(mon):
    fig, axx = plt.subplots(1, 2, sharex=True, constrained_layout=True, figsize=(10,4))
    axx[0].plot(mon.t/second, mon.Vs[0]/mV, label='Vs')
    axx[0].plot(mon.t/second, mon.Vd[0]/mV, label='Vd')
    axx[0].set_xlabel('Time (s)')
    axx[0].set_ylabel('Voltage (mV)')
    axx[0].legend()
    
    axx[1].plot(mon.t/second, mon.Cad[0]/400, label='Cad')
    axx[1].plot(mon.t/second, mon.qd[0], label='qd')
    axx[1].legend()
    
    return axx[0]

def plot_fig3_style(mon):
    fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(8,6))
    ax.plot(mon.t/ms, mon.Vs[0]/mV, label='Vs')
    ax.plot(mon.t/ms, mon.Vd[0]/mV, label='Vd')
    ax.plot(mon.t/ms, mon.Cad[0]/20-60, label='Cad')
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Voltage (mV)')
    #ax.set_xlim(0, 50)
    ax.legend()
    
    return ax

def main():
    prc, prc_group, net, mon = simulate_baseline()
    ax = plot_fig2_style(mon)
    ax.figure.savefig('./figures/pinsky_rinzel_baseline.png')
    
    prc, prc_group, net, mon = simulate_fig2a()
    ax = plot_fig2_style(mon)
    ax.figure.savefig('./figures/pinsky_rinzel_fig2a.png')
    
    prc, prc_group, net, mon = simulate_fig2a(time=0.2*second, dt=10*us)
    ax = plot_fig3_style(mon)
    ax.set_xlim(100,150)
    ax.figure.savefig('./figures/pinsky_rinzel_fig3.png')
