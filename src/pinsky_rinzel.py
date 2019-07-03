#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 11:29:42 2019

@author: alafuzof
"""

import brian2.only as br2
from brian2.units import *
from util import convert_specific_units

class PinskyRinzelCell:
    def __init__(self):
        self.general_parameters = {
            # Cell membrane area (Pinsky & Rinzel don't give a specific value)
            'membrane_area': 1*um2,
            # Soma proportion
            'pp': 0.5,
            # Coupling conductance
            'gc': 2.1*msiemens/cm2,
            # Membrane capacitance
            'Cm': 3.0*uF/cm2,
            # Reversal potentials
            'E_leak': -60*mV,
            'E_Na':    60*mV,
            'E_K':    -75*mV,
            'E_Ca':    80*mV,
            'E_NMDA':   0*mV, # This is just Vsyn in the original model
            'E_AMPA':   0*mV} # -- " -- 
        
        self.somatic_parameters = {
            # Somatic injection current, used to provide a stable resting membrane voltage
            'Is_soma': -0.5*uA/cm2,
            # Maximal conductances
            'gs_leak':  0.1*msiemens/cm2, 
            'gs_Na':   30.0*msiemens/cm2,
            'gs_Kdr':  15.0*msiemens/cm2}
        self.somatic_equations = '''
        # Currents: leak + sodium + delayed-rectifier-potassium + soma-dendrite-coupling + somatic injection
        dVs/dt = (-Is_leak - Is_Na - Is_Kdr + Is_couple/pp + Is_soma/pp)/Cm : volt
        Is_leak   = gs_leak*(Vs-E_leak) : amp
        Is_Na     = gs_Na*Minfs*Minfs*hs*(Vs-E_Na) : amp
        Is_Kdr    = gs_Kdr*ns*(Vs-E_K): amp
        Is_couple = gc*(Vd-Vs): amp
        # Instant activation (Minfs) and dynamic inactivation (hs) of sodium current
        Minfs   = alphams/(alphams+betams) : 1
        alphams = 0.32*(-46.9-Vs/mV)/(exp((-46.9-Vs/mV)/4.0)-1.0)/ms : Hz
        betams  = 0.28*(Vs/mV+19.9)/(exp((Vs/mV+19.9)/5.0)-1.0)/ms : Hz
        dhs/dt  = alphahs-(alphahs+betahs)*hs : 1
        alphahs = 0.128*exp((-43.0-Vs/mV)/18.0)/ms : Hz
        betahs  = 4.0/(1.0+exp((-20.0-Vs/mV)/5.0))/ms : Hz
        # Activation of delayed rectifier potassium current
        dns/dt  = alphans-(alphans+betans)*ns : 1
        alphans = 0.016*(-24.9-Vs/mV)/(exp((-24.9-Vs/mV)/5.0)-1.0)/ms : Hz
        betans  = 0.25*exp(-1.0-0.025*Vs/mV)/ms : Hz
        '''
        
        self.dendritic_parameters = {
                # Dendritic injection current
                'Id_dendrite': 0.0*nA,
                # Maximal conductances / conductance densities
                'gd_leak':  0.1*msiemens/cm2, 
                'gd_Ca':   10.0*msiemens/cm2,
                'gd_Kahp':  0.8*msiemens/cm2,
                'gd_KCa':  15.0*msiemens/cm2,
                'gd_AMPA':  0.0*msiemens/cm2, # Default model doesn't simulate synapses
                'gd_NMDA':  0.0*msiemens/cm2, # --"--
                # Time constants
                'tau_AMPA': 2*ms,
                'tau_NMDA': 150*ms,
                # Calcium scaler (note correct value calculated before network is initialized)
                'calcium_scaler': 1.0,
                # NMDA saturation value
                'Smax': 125.0}
        self.dendritic_equations = '''
        # Currents: leak + calcium + after-hyperpolarization-potassium + 
        # calcium-activated-potassium + soma-dendrite-coupling + AMPA + NMDA
        dVd/dt = (-Id_leak - Id_Ca - Id_Kahp - Id_KCa - Id_syn/(1-pp) + Id_couple/(1-pp) + Id_dendrite/(1-pp))/Cm : volt
        Id_leak   = gd_leak*(Vd-E_leak): amp
        Id_Ca     = gd_Ca*sd*sd*(Vd-E_Ca) : amp
        Id_Kahp   = gd_Kahp*qd*(Vd-E_K): amp
        Id_KCa    = gd_KCa*cd*chid*(Vd-E_K): amp
        Id_couple = gc*(Vs-Vd): amp
        Id_syn    = Id_AMPA + Id_NMDA: amp
        Id_AMPA   = gd_AMPA*clip(syn_AMPA,0,7000)*(Vd-E_AMPA): amp
        Id_NMDA   = gd_NMDA*clip(Si,0.0,Smax)*(Vd-E_NMDA)/(1+0.28*exp(-0.062*(Vd/mV))) : amp
        # Activation of calcium current
        dsd/dt = alphasd-(alphasd+betasd)*sd : 1
        alphasd = 1.6/(1.0+exp(-0.072*(Vd/mV-5.0)))/ms : Hz
        betasd  = 0.02*(Vd/mV+8.9)/(exp((Vd/mV+8.9)/5.0)-1.0)/ms : Hz
        dCad/dt = -0.13*Id_Ca*calcium_scaler/uamp/ms-0.075*Cad/ms : 1
        # Activation of after-hyperpolarization-potassium current
        dqd/dt = alphaqd-(alphaqd+betaqd)*qd : 1
        alphaqd = clip(0.00002*Cad,0,0.01)/ms : Hz
        betaqd  = 0.001/ms : Hz
        # Activation (cd) and saturation (chid) of calcium-activated-potassium channel
        dcd/dt  = alphacd-(alphacd+betacd)*cd : 1
        alphacd = (int(Vd/mV<=-10)*exp((Vd/mV+50.0)/11-(Vd/mV+53.5)/27)/18.975+(Vd/mV>-10)*2.0*exp((-53.5-Vd/mV)/27.0))/ms  : Hz
        betacd  = ((Vd/mV<=-10)*(2.0*exp((-53.5-Vd/mV)/27.0)-alphacd*ms)+(Vd/mV>-10)*0)/ms : Hz
        chid    = clip(Cad/250.0,0,1.0) : 1
        # AMPA activation
        dsyn_AMPA/dt = -syn_AMPA/tau_AMPA: 1
        # NMDA current
        dSi/dt = -Si/tau_NMDA : 1
        '''
        
        self.other_equations = ''
        
    def generate_neuron_group(self, n_cells, dt=None):
        parameters = {}
        parameters.update(self.general_parameters)
        parameters.update(self.somatic_parameters)
        parameters.update(self.dendritic_parameters)
        parameters = convert_specific_units(parameters, parameters['membrane_area'])
        parameters['calcium_scaler'] = 1*cm2/parameters['membrane_area']
        
        
        equations = '\n'.join([self.somatic_equations, 
                               self.dendritic_equations, 
                               self.other_equations])
    
        if dt is None:
            dt = br2.defaultclock.dt
    
        group = br2.NeuronGroup(n_cells, equations, method='rk4', namespace=parameters, dt=dt)
        group.Vs = -64.6*mV # Somatic membrane voltage
        group.Vd = -64.5*mV # Dendritic membrane voltage
        group.hs = 0.999 # Somatic Na inactivation
        group.ns = 0.001 # Somatic K activation
        group.sd = 0.009 # Dendritic Ca activation
        group.cd = 0.007 # Dendritic KCa activation
        group.qd = 0.010 # Dendritic Kahp activation
        group.Cad = 0.20 # Dendritic Ca "concentration"
        group.Si = 0.0 # Dendritic NMDA opening
        group.syn_AMPA = 0.0 # Dendritic AMPA opening
        
        return group
    
        