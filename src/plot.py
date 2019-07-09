#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 23:59:09 2019

@author: alafuzof
"""
from brian2.units import *
import matplotlib.pyplot as plt

def plot_spike_raster(spike_monitor, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, constrained_layout=True)
    
    ax.scatter(spike_monitor.t/second, spike_monitor.i)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Cell')
    
    return ax

def plot_average(t, y, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, constrained_layout=True)
        
    ax.plot(t/second, y.mean(axis=0))
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Population mean')
    
    return ax
        
