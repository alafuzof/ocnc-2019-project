#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 16:21:09 2019

@author: alafuzof
"""
from brian2.units import *
from .pinsky_rinzel_cell import PinskyRinzelCell

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
        
        # Remove E_leak, since we set it per cell below
        del self.general_parameters['E_leak']
        
        # Add conductance and time constant for GABA synapses. Note that
        # Taxidis et al. don't really use a "maximum" conductance, so the value
        # below is mainly there to make the units work out
        self.somatic_parameters['gs_GABA']  = 1.0*nsiemens
        self.somatic_parameters['tau_GABA'] = 7.0*ms
        
        # Add GABA_A synaptic current
        self.somatic_equations = replace_equation(self.somatic_equations, 'dVs/dt', 'dVs/dt = (-Is_leak - Is_Na - Is_Kdr + Is_couple/pp + Is_GABA_A/pp + Is_soma/pp)/Cm : volt')
        self.somatic_equations += '''
        # Cell specific leak reversal for heterogeneity
        E_leak : volt
        # Cell specific maximum leak conductance for heterogeneity
        gs_leak: siemens
        # Note that in Gunn et al. GABA_A currents use E_K!
        Is_GABA   = gs_GABA*clip(GABA_syn,0,7000)*(Vs-E_Cl): amp
        # Inhibitory synapses (GABA_A)
        dGABA_syn/dt=-GABA_syn/tau_GABA: 1
        '''
        
        # Taxidis et al. reduce the calcium conductance " to simulate better 
        # the in vitro characteristic firing properties of such cells, where 
        # the typical CA3 intrinsic bursting is replaced by tonic firing with 
        # frequency accommodation 
        self.dendritic_parameters['gd_Ca'] = 7*msiemens/cm2
        
        # We set the same 
        self.dendritic_equations += '''
        gd_leak = gs_leak : siemens
        '''
        
    