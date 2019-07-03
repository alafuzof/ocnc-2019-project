#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 11:28:24 2019

@author: alafuzof
"""
from brian2.units import *
#from brian2.units.fundamentalunits import Quantity

def replace_equation(equations, match, new_equation):
    split_eq = equations.split('\n')
    n_match = sum(0 if match not in x else 1 for x in split_eq) 
    if n_match != 1:
        raise ValueError(f'Only one equation should match {match}, now there are {n_matches} matches!')
    split_eq = [x if match not in x else new_equation for x in split_eq]
    return '\n'.join(split_eq)

def convert_specific_units(parameters, membrane_area):
    res = parameters.copy()
    for key,value in res.items():
        if isinstance(value, Quantity) and value.dimensions in [(uF/cm2).dim, (msiemens/cm2).dim, (uA/cm2).dim]:
            #print('Old value:', value)
            res[key] = value*membrane_area
            #print('New value:', res[key])
    return res