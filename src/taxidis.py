#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 16:21:09 2019

@author: alafuzof
"""
import brian2.only as br2
from brian2.units import *
from pinsky_rinzel import PinskyRinzelCell
from util import convert_specific_units, replace_equation

class TaxidisPinskyRinzelCell(PinskyRinzelCell):
    def __init__(self):
        super().__init__()
        
        # Cell membrane area is specified in the Taxidis et al. model
        # ionic currents are specific, so this only affects the synaptic 
        # currents
        self.general_parameters['membrane_area'] = 50_000*um2
        
        # Taxidis et al. don't specify what reversal potential they use for 
        # GABA_A channels ("V_syn" is indicated for both synapses). Let's 
        # assume the correct ion.
        self.general_parameters['E_Cl'] = -75*mV
        
        # Remove E_leak, gs_leak and gd_leak, since we set them per cell below
        del self.general_parameters['E_leak']
        del self.somatic_parameters['gs_leak']
        del self.dendritic_parameters['gd_leak']
        
        # Add conductance and time constant for GABA synapses. Note that
        # Taxidis et al. don't really use a "maximum" conductance, so the value
        # below is mainly there to give the correct scale
        self.somatic_parameters['gs_GABA_A']  = 1.0*nsiemens
        self.somatic_parameters['tau_GABA_A'] = 7.0*ms
        
        # Also set AMPA conductance to the scaling value 
        self.dendritic_parameters['gd_AMPA'] = 1.0*nsiemens
        
        # Add GABA_A synaptic current
        self.somatic_equations = replace_equation(self.somatic_equations, 'dVs/dt', 'dVs/dt = (-Is_leak - Is_Na - Is_Kdr + Is_couple/pp + Is_GABA_A/pp + Is_soma/pp)/Cm : volt')
        self.somatic_equations += '''
        # Cell specific leak reversal for heterogeneity
        E_leak : volt
        # Cell specific maximum leak conductance for heterogeneity
        gs_leak: siemens
        # Note that in Gunn et al. GABA_A currents use E_K!
        Is_GABA_A = gs_GABA_A*clip(syn_GABA_A,0,7000)*(Vs-E_Cl): amp
        # Inhibitory synapses (GABA_A)
        dsyn_GABA_A/dt=-syn_GABA_A/tau_GABA_A: 1
        '''
        
        # We set the same leakage conductance in the dendrite as in the soma
        self.dendritic_equations += '''
        gd_leak = gs_leak : siemens
        '''
    
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
        # Initial conditions are Gaussian distributed around the "steady-state"
        # with 10% SD
        group.Vs = '(0.1*randn()*64.6 - 64.6)*mV' # Somatic membrane voltage
        group.Vd = '(0.1*randn()*64.5 - 64.5)*mV' # Dendritic membrane voltage
        # FIXME: Check if it's okay for activation variables to be outside [0,1] range
        group.hs = '(0.1*randn()*0.999 + 0.999)' # Somatic Na inactivation
        group.ns = '(0.1*randn()*0.001 + 0.001)' # Somatic K activation
        group.sd = '(0.1*randn()*0.009 + 0.009)' # Dendritic Ca activation
        group.cd = '(0.1*randn()*0.007 + 0.007)' # Dendritic KCa activation
        group.qd = '(0.1*randn()*0.010 + 0.010)' # Dendritic Kahp activation
        group.Cad = '(0.1*randn()*0.20 + 0.20)' # Dendritic Ca "concentration"
        
        # Taxidis model has no NMDA channels, so we set their activaty to 0
        group.Si = 0.0 # Dendritic NMDA activity
        # AMPA and GABA activity start at zero
        group.syn_AMPA = 0.0 # Dendritic AMPA activity
        group.syn_GABA_A = 0.0 # Somatic GABA activity
        
        # Maximum leak conductance and reversal are Gaussian distributed to 
        # increase heterogeneity
        group.gs_leak = '(0.1*randn()*0.1*msiemens*membrane_area/cm2 + 0.1*msiemens*membrane_area/cm2)'
        group.E_leak = '(0.1*randn()*60 - 60)*mV'
        
        return group
        
def generate_ca1_pyramidal_group(n_cells=1000, dt=None):
    tpr_cell = TaxidisPinskyRinzelCell()
    
    # Taxidis et al. reduce the calcium conductance for CA1 cells "to 
    # simulate better the in vitro characteristic firing properties of such 
    # cells, where the typical CA3 intrinsic bursting is replaced by tonic 
    # firing with frequency accommodation 
    tpr_cell.dendritic_parameters['gd_Ca'] = 7*msiemens/cm2
    
    return tpr_cell.generate_neuron_group(n_cells, dt)
    
def generate_ca3_pyramidal_group(n_cells=1000, dt=None):
    tpr_cell = TaxidisPinskyRinzelCell()
    
    return tpr_cell.generate_neuron_group(n_cells, dt)