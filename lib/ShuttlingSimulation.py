#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#########################################################################
#                                                                       #
#               Simulation of the relaxometry experiments               #
#                                                                       #
#########################################################################
import numpy as np
from scipy import linalg
from random import uniform
import time
import copy

import FitFunctions as FitF

from _RelaxMat import RelaxMat as RM


def make_field_list(self, increment):
    """
    make a single field list for all the experiments:
       1) considers the highest height the sample is shuttled to;
       2) compute the field at every positions separated by the increment;
       3) create the delay-time list for each experiments;
       4) use the slowest speed to compute these delays.

    Parameters
    ----------
    increment : TYPE: float
        DESCRIPTION: distance increment to compute the propagator at during the
        simulations.

    Returns
    -------
    field_list : TYPE: array
        DESCRIPTION: field list from the center of the magnet up to the highest
        height the sample is shuttled to.
    delays : TYPE: dictionnary
        DESCRIPTION: delays used to compute the propagators when simulating the
        experiments.

    """
    #get all the positions along the trajectory
    all_heights = [self.set_up[exp]['height'] for exp in self.set_up.keys()]
    max_height = max(all_heights)
    
    positions =  np.arange(increment, max_height, increment)
    for h in all_heights:
        positions = np.append(positions, h)
    positions = np.sort(positions)
    
    #compute field at each positions
    field_list = [FitF.Calc_B0(p, self.B0_cal_coeff) for p in positions]
    
    #get the delay between sucessive positions
    delays = {'up': {}, 'down': {}}
    if self.shuttling_type == "Constant Speed":
        for exp in self.set_up.keys():
            speed_up = self.set_up[exp]['height'] / self.set_up[exp]['SLF']
            speed_down = self.set_up[exp]['height'] / self.set_up[exp]['SHF']
            
            delays['up'][exp] = [positions[0] / speed_up]
            delays['down'][exp] = [positions[0] / speed_down]
            for counter, p in enumerate(positions[1:]):
                if p > self.set_up[exp]['height']:
                    break
                else:
                    delays['up'][exp].append( (positions[counter+1] - positions[counter]) / speed_up )
                    delays['down'][exp].append( (positions[counter+1] - positions[counter]) / speed_down )
    else:
        for exp in self.set_up.keys():
            acceleration_up = 4.0*self.set_up[exp]['height'] / self.set_up[exp]['SLF']**2
            acceleration_down = 4.0*self.set_up[exp]['height'] / self.set_up[exp]['SHF']**2

            delays['up'][exp] = []
            delays['down'][exp] = []
            for p in positions:
                if p > self.set_up[exp]['height']:
                    break
                else:
                    if p < self.set_up[exp]['height'] / 2.:
                        dt_up = np.sqrt(2. * p / acceleration_up)
                        dt_down = np.sqrt(2. * p / acceleration_up)
                    else:
                        dt_up = np.sqrt(2. / acceleration_up * (4.*self.set_up[exp]['height'] / self.set_up[exp]['SLF'] -
                                                          p - self.set_up[exp]['height']))
                        dt_down = np.sqrt(2. / acceleration_down * (4.*self.set_up[exp]['height'] / self.set_up[exp]['SHF'] -
                                                          p - self.set_up[exp]['height']))
                delays['up'][exp].append(dt_up)
                delays['down'][exp].append(dt_down)
        
    return field_list, delays


def Propagator_Diagonalization(relax_mat, delay):
    """
    calculated the propagator by diagonalizing the relaxation matrix

    Parameters
    ----------
    relax_mat : TYPE: array
        DESCRIPTION: relaxation matrix.
    delay : TYPE: float
        DESCRIPTION: delay for the evolution.

    Returns
    -------
    P : TYPE: array
        DESCRIPTION: propagator.

    """
    relax_mat = np.asarray(relax_mat)
    P = linalg.expm(-delay * relax_mat)
    
    return P


def Propagator_Exponential(relax_mat, delay):
    """
    calculated the propagator by computing the exponential of the relaxation matrix

    Parameters
    ----------
    relax_mat : TYPE: array
        DESCRIPTION: relaxation matrix.
    delay : TYPE: float
        DESCRIPTION: delay for the evolution.

    Returns
    -------
    P : TYPE: array
        DESCRIPTION: propagator.

    """
    relax_mat = np.asarray(relax_mat)
    eigenvalue, eigenvector = np.linalg.eig(relax_mat)
    
    P = eigenvector @ np.diag(np.exp(-delay * eigenvalue)) @ np.linalg.inv(eigenvector)
    
    return P


def Calc_Propagator(exp, param, tauc, other_input, Static_MagField, relax_mat_shuttle, delays_shuttle, set_up, B0_low_fields, calculation_method):
    """
    calculates the propagator for a given experiment

    Parameters
    ----------
    exp : TYPE: str
        DESCRIPTION: experiment number.
    param : TYPE: array
        DESCRIPTION: MCMC parameters.
    tauc : TYPE: float
        DESCRIPTION: value of tau_c.
    other_input : TYPE: dictionnary
        DESCRIPTION: other inputs to calculate the relaxation matrix.
    Static_MagField : TYPE: float
        DESCRIPTION: static magnetic field of the relaxometry spectrometer.
    relax_mat_shuttle : TYPE: array
        DESCRIPTION: list of relaxation matrices along the sample trajectory.
    delays_shuttle : TYPE: dictionnary
        DESCRIPTION: list of delays between each position in field_list.
    set_up : TYPE: dictionnary
        DESCRIPTION: experimental setup of the relaxometry experiments.
    B0_low_fields : TYPE: dictionnary
        DESCRIPTION: low fields at which the relaxometry decay is recorded.
    calculation_method : TYPE: function
        DESCRIPTION: optimized approach to compute the propagators.

    Returns
    -------
    Propagators : TYPE: dictionnary
        DESCRIPTION: list of propagators during the relaxometry experiment.

    """
    Propagators = {}
    #High Field
    relaxmat_HF = RM(Static_MagField, param, tauc, other_input)
    Propagators['HF_1'] = calculation_method(relaxmat_HF, set_up[exp]['WTHF'])
    Propagators['HF_2'] = calculation_method(relaxmat_HF, set_up[exp]['d25'])
    
    #Low Field
    relaxmat_LF = RM(B0_low_fields[exp], param, tauc, other_input)
    Propagators['LF'] = {}
    for vc in set_up[exp]['vc']:
        Propagators['LF'][vc] = calculation_method(relaxmat_LF, set_up[exp]['LF_times'][vc])
        
    #shuttling
    prop_shuttling_up = []
    for count, delay in enumerate(delays_shuttle['up'][exp]):
        prop_shuttling_up.append(calculation_method(relax_mat_shuttle[count], delay))
    prop_shuttling_down = []
    for count, delay in enumerate(delays_shuttle['down'][exp]):
        prop_shuttling_down.append(calculation_method(relax_mat_shuttle[count], delay))

    Propagators['LF->HF'] = prop_shuttling_down[0]
    for p in prop_shuttling_down[1:]:
        Propagators['LF->HF'] = Propagators['LF->HF'] @ p
    Propagators['HF->LF'] = prop_shuttling_up[0]
    for p in prop_shuttling_up[1:]:
        Propagators['HF->LF'] = p @ Propagators['HF->LF']
        
    return Propagators
        

def Expected_Values(param, tauc, other_input, Static_MagField, field_list, delays_shuttle, set_up, B0_low_fields, PosAuto, calculation_method):
    """
    calculates the expected value at the end of the relaxometry experiments

    Parameters
    ----------
    param : TYPE: array
        DESCRIPTION: MCMC parameters.
    tauc : TYPE: float
        DESCRIPTION: value of tau_c.
    other_input : TYPE: dictionnary
        DESCRIPTION: other inputs to calculate the relaxation matrix.
    Static_MagField : TYPE: float
        DESCRIPTION: static magnetic field of the relaxometry spectrometer.
    field_list : TYPE: array
        DESCRIPTION: list of fields along the sample trajectory.
    delays_shuttle : TYPE: dictionnary
        DESCRIPTION: list of delays between each position in field_list.
    set_up : TYPE: dictionnary
        DESCRIPTION: experimental setup of the relaxometry experiments.
    B0_low_fields : TYPE: dictionnary
        DESCRIPTION: low fields at which the relaxometry decay is recorded.
    PosAuto : TYPE: int
        DESCRIPTION: position of the rate of interest in the relaxation matrix.
    calculation_method : TYPE: function
        DESCRIPTION: optimized approach to compute the propagators.

    Returns
    -------
    expected_values : TYPE: dictionnary
        DESCRIPTION: list of expectation values during the relaxometry experiment.

    """
    #Shuttle periods
    relax_mat_shuttle = []
    for b0_shuttle in field_list:
        relax_mat_shuttle.append(RM(b0_shuttle, param, tauc, other_input))
    
    expected_values = {}
    for exp in set_up.keys():
        expected_values[exp] = {}
        Propagators = Calc_Propagator(exp, param, tauc, other_input, Static_MagField, relax_mat_shuttle, delays_shuttle, set_up, B0_low_fields, calculation_method)
        
        for vc in set_up[exp]['vc']:
            total_propagrator = Propagators['HF_2'] @ Propagators['LF->HF'] @ Propagators['LF'][vc] @ Propagators['HF->LF'] @ Propagators['HF_1']
            expected_values[exp][vc] = total_propagrator[PosAuto][PosAuto]
           
    return expected_values


def Expected_Values_For_Plot(param, tauc, other_input, Static_MagField, field_list, delays_shuttle, set_up, B0_low_fields, PosAuto, lf_times, exp, calculation_method):
    """
    calculates the expected value at the end of the relaxometry experiments,
    only for figures

    Parameters
    ----------
    param : TYPE: array
        DESCRIPTION: MCMC parameters.
    tauc : TYPE: float
        DESCRIPTION: value of tau_c.
    other_input : TYPE: dictionnary
        DESCRIPTION: other inputs to calculate the relaxation matrix.
    Static_MagField : TYPE: float
        DESCRIPTION: static magnetic field of the relaxometry spectrometer.
    field_list : TYPE: array
        DESCRIPTION: list of fields along the sample trajectory.
    delays_shuttle : TYPE: dictionnary
        DESCRIPTION: list of delays between each position in field_list.
    set_up : TYPE: dictionnary
        DESCRIPTION: experimental setup of the relaxometry experiments.
    B0_low_fields : TYPE: dictionnary
        DESCRIPTION: low fields at which the relaxometry decay is recorded.
    PosAuto : TYPE: int
        DESCRIPTION: position of the rate of interest in the relaxation matrix.
    lf_times : TYPE: array
        DESCRIPTION: relaxation delays.
    exp : TYPE: str
        DESCRIPTION: experiment number.
    calculation_method : TYPE: function
        DESCRIPTION: optimized approach to compute the propagators.

    Returns
    -------
    TYPE: array
        DESCRIPTION: expected value for the decay of experiment exp.

    """
    set_up_new = copy.deepcopy(set_up)
    set_up_new[exp]['vc'] = lf_times
    set_up_new[exp]['LF_times'] = {t: t for t in lf_times}

    #Shuttle periods
    relax_mat_shuttle = []
    for b0_shuttle in field_list:
        relax_mat_shuttle.append(RM(b0_shuttle, param, tauc, other_input))
        
    Propagators = Calc_Propagator(exp, param, tauc, other_input, Static_MagField, relax_mat_shuttle, delays_shuttle, set_up_new, B0_low_fields, calculation_method)
    
    expected_values = []
    for vc in lf_times:
        total_propagrator = Propagators['HF_2'] @ Propagators['LF->HF'] @ Propagators['LF'][vc] @ Propagators['HF->LF'] @ Propagators['HF_1']
        expected_values.append(total_propagrator[PosAuto][PosAuto])
           
    return np.asarray(expected_values)


def optimize_shuttling_increment(self, PosAuto):
    """
    optimizes the shuttling increment distance

    Parameters
    ----------
    PosAuto : TYPE: int
        DESCRIPTION: position of the rate of interest in the relaxation matrix.

    Returns
    -------
    TYPE: float
        DESCRIPTION: distance increment.

    """
    print()
    print("Optimizing shuttling increment...")
    #initialize
    lowest_field_idx = list(self.B0_low_field.values()).index(min(list(self.B0_low_field.values())))
    exp_lowest_field = list(self.B0_low_field.keys())[lowest_field_idx]
    
    RandomParam = [[] for i in range(10)]
    for i in range(10):
        for name in self.bonds_starting_point.keys():
            for P in range(len(self.bonds_starting_point[name])):
                RandomParam[i].append(uniform(self.bonds_starting_point[name][P][0], self.bonds_starting_point[name][P][1]))
                
    Increment = 1e-3    #starting increment in meter
    field_list, delays = make_field_list(self, Increment)
    
    sub_setup = {exp_lowest_field: self.set_up[exp_lowest_field]}
    sub_B0LF = {exp_lowest_field: self.B0_low_field[exp_lowest_field]}
    
    ExpVal_Init = []
    aa = list(self.other_inputs.keys())[0]
    for param in RandomParam:
        ExpVal_Init_Full = Expected_Values(param, self.TauC, self.other_inputs[aa], self.Static_MagField, field_list, delays, sub_setup, sub_B0LF, PosAuto, Propagator_Diagonalization)
        ExpVal_Init.append(ExpVal_Init_Full[exp_lowest_field][max(sub_setup[exp_lowest_field]['vc'])])
    print(f" Initial value : {Increment} m")
    
    #optimize
    count = 0
    while True:
        count += 1
        Increment = round(Increment + 1e-3, 3)
        field_list, delays = make_field_list(self, Increment)
        
        for counter, param in enumerate(RandomParam):
            ExpVal_Full = Expected_Values(param, self.TauC, self.other_inputs[aa], self.Static_MagField, field_list, delays, sub_setup, sub_B0LF, PosAuto, Propagator_Diagonalization)
            ExpVal = ExpVal_Full[exp_lowest_field][max(sub_setup[exp_lowest_field]['vc'])]
        
            if abs(ExpVal-ExpVal_Init[counter])/ExpVal_Init[counter] > 0.01:
                print(f"Final used increment: {round(Increment - 1e-3, 3)} m")
                return round(Increment - 1e-3, 3)
            
        print(f" Updated value {count} : {Increment} m")
        
        
def optimize_calc_propagator(self):
    """
    sets the calculation method for matrix exponentials
    
    """
    print()
    print("Choosing propagator calculation method...")
    
    aa = list(self.other_inputs.keys())[0]
    RandomParam = []
    for name in self.bonds_starting_point.keys():
        for P in range(len(self.bonds_starting_point[name])):
            RandomParam.append(uniform(self.bonds_starting_point[name][P][0], self.bonds_starting_point[name][P][1]))
    relax_mat = RM(5.0, RandomParam, self.TauC, self.other_inputs[aa])
    dt = 1e-3

    print(" Method 1: matrix exponential")
    start = time.time()
    for i in range(10000):
        _ = Propagator_Exponential(relax_mat, dt)
    end = time.time()
    Duration_exp = end-start
    print(f"  time for 10,000 iterations: {round(Duration_exp, 2)} s")
    
    print(" Method 2: matrix diagonalization")
    start = time.time()
    for i in range(10000):
        _ = Propagator_Diagonalization(relax_mat, dt)
    end = time.time()
    Duration_diag = end-start
    print(f"  time for 10,000 iterations: {round(Duration_diag, 2)} s")
    
    if min(Duration_exp, Duration_diag) == Duration_exp:
        self.PropFunction = Propagator_Exponential
        print("Choosing calculation method done. Method 1 chosen.")
    else:
        self.PropFunction = Propagator_Diagonalization
        print("Choosing calculation method done. Method 2 chosen.")
