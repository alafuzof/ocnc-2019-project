#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 17:19:47 2019

@author: alafuzof
"""

import brian2.only as br2
from brian2.units import *
from util import convert_specific_units

class WangBuzsakiCell():
    def __init__(self):
        self.parameters = {
            #
            'membrane_area': 1*cm2,
            # Membrane capacitance
            'Cm': 1*uF/cm2,
            # Maximum conductances
            'g_Na':  35.0*msiemens/cm2,
            'g_K':    9.0*msiemens/cm2,
            'g_leak': 0.1*msiemens/cm2,
            'g_syn':  0.1*msiemens/cm2,
            # Reversal potentials
            'E_Na':    55*mV,
            'E_K':    -90*mV,
            'E_leak': -65*mV,
            'E_syn':  -75*mV,
            # Synaptic decay time constant
            'tau_syn': 10*ms,
            # Phi is used to modulate the speed of activation and inactivation kinetics
            'phi': 5.0}
        
        self.equations = '''
        # Currents: sodium + delayed-rectifier-potassium + leak + synaptic + injected
        dV/dt = (-I_Na-I_K-I_leak-I_syn+I_inj)/Cm : volt
        I_Na   = g_Na*m_inf**3*h*(V-E_Na): amp
        I_K    = g_K*n**4*(V-E_K): amp
        I_leak = g_leak*(V-E_leak): amp
        I_syn  = g_syn*s*(V-E_syn): amp
        I_inj: amp
        # Instant activation (m_inf) and dynamic inactivation (h) of sodium current
        m_inf = alpha_m/(alpha_m+beta_m) : 1
        alpha_m = -0.1/mV*(V+35*mV)/(exp(-0.1/mV*(V+35*mV))-1)/ms : Hz
        beta_m = 4*exp(-(V+60*mV)/(18*mV))/ms : Hz
        dh/dt = phi*(alpha_h*(1-h)-beta_h*h) : 1
        alpha_h = 0.07*exp(-(V+58*mV)/(20*mV))/ms : Hz
        beta_h = 1./(exp(-0.1/mV*(V+28*mV))+1)/ms : Hz
        # Activation of delayed-rectifier potassium current
        dn/dt = phi*(alpha_n*(1-n)-beta_n*n) : 1
        alpha_n = -0.01/mV*(V+34*mV)/(exp(-0.1/mV*(V+34*mV))-1)/ms : Hz
        beta_n = 0.125*exp(-(V+44*mV)/(80*mV))/ms : Hz
        # Synaptic activation. Note the part that depends on pre-synaptic 
        # voltage is segregated into the synaptic model. tau_syn is the 1/beta
        ds/dt = s_tot - s/tau_syn : 1
        s_tot: Hz
        '''
        
        self.spike_threshold = -52*mV
        
        
        self.synaptic_parameters = {
            'alpha': 12*ms**-1,
            'theta_syn': 0*mV}
        self.synaptic_equations = '''
        ds_syn/dt = alpha*F*(1.0-s_post) : 1
        F = 1.0/(1.0+exp(-(V_pre/mV-theta_syn/mV)/2)): 1
        s_tot_post = s_syn : 1
        '''
        
    def generate_neuron_group(self, n_cells, dt=None):
        if dt is None:
            dt = br2.defaultclock.dt
            
        parameters = convert_specific_units(self.parameters, self.parameters['membrane_area'])
        equations = self.equations[:]
        
        group = br2.NeuronGroup(n_cells, equations,
                                threshold=f'V >= {self.spike_threshold/mV}*mV',
                                refractory=f'V >= {self.spike_threshold/mV}*mV',
                                namespace=parameters, 
                                method='exponential_euler', 
                                dt=dt)
        
        #group.V = -70*mV
        group.V = -70*mV
        group.h = 1.0
        #group.n = 0.001
        #group.s_tot = 0.0/ms
        #group.s = 0.0
        
        return group
    
    def generate_synapses(self, P, Q, dt=None):
        if dt is None:
            dt = br2.defaultclock.dt
        
        parameters = self.synaptic_parameters
        equations = self.synaptic_equations
        
        synapses = br2.Synapses(P, Q, model=equations, namespace=parameters, method='rk4', dt=dt)
        
        return synapses
    
    def tmp(self):
        self.parameters = {
            #
            #'membrane_area': 1*cm2,
            # Membrane capacitance
            'Cm': 1*uF,
            # Maximum conductances
            'g_Na':  35.0*msiemens,
            'g_K':    9.0*msiemens,
            'g_leak': 0.1*msiemens,
            'g_syn':  0.1*msiemens,
            # Reversal potentials
            'E_Na':    55*mV,
            'E_K':    -90*mV,
            'E_leak': -65*mV,
            'E_syn':  -75*mV,
            # Synaptic decay time constant
            'tau_syn': 10*ms,
            # Phi is used to modulate the speed of activation and inactivation kinetics
            'phi': 5.0
            }
        