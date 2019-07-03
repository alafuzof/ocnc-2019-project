#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 16:45:33 2019

@author: alafuzof
"""


class GunnCell(TaxidisCell):
    def __init__(self):
        super().__init__()
        
        # The slightly higher value corresponds to 3.5*msiemens/cm2 with 50000 um2 cells
        self.general_parameters['gc'] = 1.75*usiemens, 
        
        {'Is_soma': -0.3*nA} # Note this is the Gunn et al. value
        
        # This is used for counting arriving spikes at each cell!?
        '''
        counter: 1
        '''