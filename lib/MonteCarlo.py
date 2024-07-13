#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################################
#                                                    #
#               Monte Carlo simulation               #
#                                                    #
######################################################
from random import uniform
from multiprocessing import Pool

import Parameters as ParamFile
import FigOutputs as FigOut
import ShuttlingSimulation as ShSim

import numpy as np
import emcee
from math import exp, log
import sys


def scaling_factor(simulated_intensity, relaxometry_data, AA):
    """
    computes the scaling factors such that the simulated intensity decays matches the experimental intensities

    Parameters
    ----------
    simulated_intensity : TYPE: dictionnary
        DESCRIPTION: simulated intensity decays for one residue.
    intensity_relaxometry : TYPE: dictionnary
        DESCRIPTION: experimental intensities and errors.
    AA : TYPE: str
        DESCRIPTION: residue.

    Returns
    -------
    scaling : TYPE: dictionnary
        DESCRIPTION: scaling factor for each experiments and one particular residue.

    """
    scaling = {}
    
    for exp_num in relaxometry_data.keys():
        Scaling_Num, Scaling_Den = 0.0, 0.0
        for vc in relaxometry_data[exp_num][AA].keys():
            if relaxometry_data[exp_num][AA][vc]['data'] != "NA":
                Scaling_Num += relaxometry_data[exp_num][AA][vc]['data'] * simulated_intensity[exp_num][vc] / (relaxometry_data[exp_num][AA][vc]['error'])**2
                Scaling_Den += (simulated_intensity[exp_num][vc] / relaxometry_data[exp_num][AA][vc]['error'])**2
        scaling[exp_num] = Scaling_Num / Scaling_Den
            
    return scaling


def lnlike_high_field(X, AA, relax_data, av_scaling, std_scaling, tauc, other_inputs):
    """
    lnlike function for high-field data. See the emcee doc for more info.

    Parameters
    ----------
    X : TYPE: array
        DESCRIPTION: values of the free parameters of the MCMC.
    AA : TYPE: str
        DESCRIPTION: residue.
    relax_data : TYPE: dictionary
        DESCRIPTION: rates and errors, shifted to have an average and std of 0 and 1.
    av_scaling : TYPE: dictionary
        DESCRIPTION: averages for the scaling for high-field data only.
    std_scaling : TYPE: dictionary
        DESCRIPTION: standard deviations for the scaling for high-field data only.
    tauc : TYPE: float
        DESCRIPTION: value of correlation time tc.
    other_inputs : TYPE: dictionary
        DESCRIPTION: other inputs.

    Returns
    -------
    distCalc : TYPE: float
        DESCRIPTION: value of the lnlike function for high-field data only.

    """
    relaxation_function = ParamFile.ImportFunc()
    
    distCalc = 0.0
    expLnF = exp(2.0*X[-1])
    
    for Func in relax_data.keys():
        for B in relax_data[Func][AA].keys():
            if 'data' in relax_data[Func][AA][B].keys() and relax_data[Func][AA][B]['data'] != 'NA':
                exp_data = relax_data[Func][AA][B]['data']
                exp_err = relax_data[Func][AA][B]['error']
                model = (relaxation_function[Func](B, X[:-1], tauc, other_inputs[AA]) - av_scaling[Func][AA]) / std_scaling[Func][AA]
                
                inv_sigma2 = exp_err**2 + model**2 * expLnF
                distCalc += -1.0/2.0 * ((exp_data - model)**2 / inv_sigma2 + log(2.0*np.pi*inv_sigma2))
                
    return distCalc


def lnlike_relaxometry(X, AA, intensity_relaxometry, intensity_relaxometry_shifted, av_scaling, std_scaling,
           tauc, other_inputs, mag_field, field_list, Delays, set_up, B0_low_fields, prop_function):
    """
    lnlike function for relaxometry data. See the emcee doc for more info.

    Parameters
    ----------
    X : TYPE: array
        DESCRIPTION: values of the free parameters of the MCMC.
    AA : TYPE: str
        DESCRIPTION: residue.
    intensity_relaxometry : TYPE: dictionary
        DESCRIPTION: intensities and errors.
    intensity_relaxometry : TYPE: dictionary
        DESCRIPTION: intensities and errors, shifted to have an average and std of 0 and 1.
    av_scaling : TYPE: dictionary
        DESCRIPTION: averages for the scaling of the relaxometry data only.
    std_scaling : TYPE: dictionary
        DESCRIPTION: standard deviations for the scaling of the relaxometry data only..
    tauc : TYPE: float
        DESCRIPTION: value of correlation time tc.
    other_inputs : TYPE: dictionary
        DESCRIPTION: other inputs.
    mag_field : TYPE: float
        DESCRIPTION: static magnetic field where relaometry experiments are recorded.
    field_list : TYPE: dictionary
        DESCRIPTION: contains the fields along the sample trajectory.
    Delays : TYPE: dictionary
        DESCRIPTION: contains the delays between successive positions in the relaxometry experiments.
    set_up : TYPE: dictionary
        DESCRIPTION: relaxometry experimental setup.
    B0_low_fields : TYPE: dictionary
        DESCRIPTION: low fields.
    prop_function : TYPE: function
        DESCRIPTION: optimized function to compute the propagator.

    Returns
    -------
    distCalc : TYPE: float
        DESCRIPTION: value of the lnlike function for relaxometry data only.

    """
    distCalc = 0.0
    expLnF = exp(2.0*X[-1])

    simulated_intensities = ShSim.Expected_Values(X[:-1], tauc, other_inputs[AA], mag_field, field_list, Delays, set_up, B0_low_fields, ParamFile.PositionAuto, prop_function)
    intensity_scaling = scaling_factor(simulated_intensities, intensity_relaxometry, AA)
    
    for exp_num in set_up.keys():
        N_vc, dist = 0., 0.
        for vc in set_up[exp_num]['vc']:
            if intensity_relaxometry[exp_num][AA][vc]['data'] != 'NA':
                N_vc += 1.
                exp_data = intensity_relaxometry_shifted[exp_num][AA][vc]['data']
                exp_err = intensity_relaxometry_shifted[exp_num][AA][vc]['error']
                model = (intensity_scaling[exp_num] * simulated_intensities[exp_num][vc] - av_scaling[exp_num][AA]) / std_scaling[exp_num][AA]
            
                inv_sigma2 = exp_err**2 + model**2 * expLnF
                dist += -1.0/2.0 * ((exp_data - model)**2 / inv_sigma2 + log(2.0*np.pi*inv_sigma2))
                    
        distCalc += dist / N_vc
                
    return distCalc

    
def lnprob(X, AA, intensity_relaxometry, intensity_relaxometry_shifted, relax_data, av_scaling, std_scaling,
           tauc, other_inputs, mag_field, field_list, Delays, set_up, B0_low_fields, prop_function):
    """
    lnprb function to evaluate the validity of the MCMC step. See the emcee documentations for more info.

    Parameters
    ----------
    X : TYPE: array
        DESCRIPTION: values of the free parameters of the MCMC.
    AA : TYPE: str
        DESCRIPTION: residue.
    intensity_relaxometry : TYPE: dictionary
        DESCRIPTION: intensities and errors.
    intensity_relaxometry_shifted : TYPE: dictionary
        DESCRIPTION: intensities and errors, shifted to have an average and std of 0 and 1.
    relax_data : TYPE: dictionary
        DESCRIPTION: rates and errors, shifted to have an average and std of 0 and 1.
    av_scaling : TYPE: dictionary
        DESCRIPTION: averages for the scaling.
    std_scaling : TYPE: dictionary
        DESCRIPTION: standard deviations for the scaling.
    tauc : TYPE: float
        DESCRIPTION: value of correlation time tc.
    other_inputs : TYPE: dictionary
        DESCRIPTION: other inputs.
    mag_field : TYPE: float
        DESCRIPTION: static magnetic field where relaometry experiments are recorded.
    field_list : TYPE: dictionary
        DESCRIPTION: contains the fields along the sample trajectory.
    Delays : TYPE: dictionary
        DESCRIPTION: contains the delays between successive positions in the relaxometry experiments.
    set_up : TYPE: dictionary
        DESCRIPTION: relaxometry experimental setup.
    B0_low_fields : TYPE: dictionary
        DESCRIPTION: low fields.
    prop_function : TYPE: function
        DESCRIPTION: optimized function to compute the propagator.

    Returns
    -------
    TYPE: float
        DESCRIPTION: value of the lnlike function.

    """
    lp = ParamFile.Cons(X)
    if lp != 0.0:
        return -np.inf
    else:
        return lnlike_relaxometry(X, AA, intensity_relaxometry, intensity_relaxometry_shifted, av_scaling['intensities'], std_scaling['intensities'],
                          tauc, other_inputs, mag_field, field_list, Delays, set_up, B0_low_fields, prop_function) +\
                lnlike_high_field(X, AA, relax_data, av_scaling['high field'], std_scaling['high field'], tauc, other_inputs)


def number_walker_check(number_walker_b, nParam):
    """
    checks the number of chains in the MCMC, and updates if necessary. See the emcee doc for more info.

    Parameters
    ----------
    number_walker_b : TYPE: int
        DESCRIPTION: number of chains as given by user.
    nParam : TYPE: int
        DESCRIPTION: number of parameters in the MCMC.

    Returns
    -------
    number_walker : TYPE: int
        DESCRIPTION: new (potentially identical) number of chains.

    """
    if number_walker_b < 2*nParam:
        number_walker = 2*nParam
        print("")
        print("WARNING")
        print(f"The number of chain is changed to *{number_walker}* because it must be at least 2-times higher than the number of parameters")
        print("")
    else:
        if number_walker_b % 2 == 0:
            number_walker = number_walker_b
        else:
            number_walker = number_walker_b + 1
            print("")
            print("WARNING")
            print(f"The number of chain is changed to *{number_walker}* because it must be even")
            print("")
        
    return number_walker


def initialize_mcmc(number_walker, number_param, bounds):
    """
    initializes the MCMC according to bounds provided by the user. Checks are performed, exists if they fail.

    Parameters
    ----------
    number_walker : TYPE: int
        DESCRIPTION: number of chains.
    number_param : TYPE: int
        DESCRIPTION: number of parameters in the MCMC.
    bounds : TYPE: dictionnary
        DESCRIPTION: bounds provided by the user.

    Returns
    -------
    init_chain : TYPE: array
        DESCRIPTION: initial parameter values for the MCMC.

    """
    init_chain = np.empty(shape=(number_walker, number_param), dtype=float)
    
    for w in range(number_walker):
        count_p = 0
        for name in bounds.keys():
            for i in range(len(bounds[name])):
                init_chain[w][count_p] = uniform(bounds[name][i][0], bounds[name][i][1])
                count_p += 1
        if ParamFile.Cons(init_chain[w]) != 0.:
            print("")
            print("ERROR")
            print("The given bounds to initialize the MCMC do not match the constrain function.")
            print("")
            sys.exit()
            
    return init_chain
    
    
def Markov_Chain_Monte_Carlo(self, AA):
    """
    MCMC function. See the emcee doc for more info

    Parameters
    ----------
    AA : TYPE: str
        DESCRIPTION: residue.

    Returns
    -------
    Results : TYPE: dictionnary
        DESCRIPTION: contains the mean and positive and negative error for each parameter.
    MAF : TYPE: float
        DESCRIPTION: mean acceptance fraction of the MCMC.
    full_samples : TYPE: array
        DESCRIPTION: the whole MCMC trajectory.

    """
    Labels = []
    for name in ParamFile.Names.keys():
        for p in ParamFile.Names[name]:
            Labels.append(p)
    Labels.append("f")

    cutoff = min(int(0.2*self.number_mcmc_steps), 500)
    number_param = len(ParamFile.Names['OrderParam']) + len(ParamFile.Names['CorrTimes']) + len(ParamFile.Names['others']) + 1
    init_chain = initialize_mcmc(self.number_walker, number_param, self.bonds_starting_point)
    
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(self.number_walker, number_param, lnprob,
                                        args=(AA, self.intensities, self.intensities_scaled, self.HF_data_scaled, self.average, self.std,
                                              self.TauC, self.other_inputs, self.Static_MagField,
                                              self.shuttling_fields, self.shuttling_delays, self.set_up, self.B0_low_field, self.PropFunction),
                                        pool = pool)
        try:
            sampler.run_mcmc(init_chain, self.number_mcmc_steps, progress=True);
        except:
            try:
                width = 30
                for i, result in enumerate(sampler.sample(init_chain, iterations=self.number_mcmc_steps)):
                    n = int((width+1) * float(i) / self.number_mcmc_steps)
                    sys.stdout.write("\r[{0}{1}]".format('#' * n, ' ' * (width - n)))
                sys.stdout.write("\n")
            except:
                sampler.run_mcmc(init_chain, self.number_mcmc_steps);
    
    samples = sampler.chain[:, cutoff:, :].reshape((-1, number_param))
    full_samples = sampler.chain[:, :, :].reshape((-1, number_param))
            
    Mean = np.percentile(samples, 50, axis=0)
    Pos_Error = np.percentile(samples, 84, axis=0) - Mean
    Neg_Error = Mean - np.percentile(samples, 16, axis=0)
    Results = {'Mean': Mean, '+': Pos_Error, '-': Neg_Error}
            
    MAF = np.mean(sampler.acceptance_fraction)
            
    print("")
    print('Param\tMean\t+error\t-error')
    for P in range(number_param):
        print(f"{Labels[P]}\t{'{:.3e}'.format(Mean[P])}\t{'{:.3e}'.format(Pos_Error[P])}\t{'{:.3e}'.format(Neg_Error[P])}")
    print("")
    print(f"Mean Acceptance Fraction: {round(MAF, 4)}")

    FigOut.figure_trajectory(self.pdf_trajectories, f'{self.dir_fit_output}/Trajectories', sampler, Mean, Labels, cutoff, AA)
    FigOut.figure_correlation_plot(self.pdf_correlations, f'{self.dir_fit_output}/Correlations', samples, Mean, Labels, AA)
            
    return Results, MAF, full_samples