#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#########################################################
#                                                       #
#               Functions used during fit               #
#                                                       #
#########################################################
import math
import numpy as np
from scipy.optimize import curve_fit

import Parameters as ParamFile
RelaxFunc = ParamFile.ImportFunc()


#################################################### Exponential decay ####################################################
def exp(x, a, b):
    return b*np.exp(-a*x)


#################################################### Fit exp decay ####################################################
def fit_intensity_decay(time, data):
    param, cov = curve_fit(exp, time, data)
    
    rate = param[0]
    pre_exp = param[1]
    
    err_mat = np.sqrt(np.diag(cov))
    err = err_mat[0]
    pre_exp_err = err_mat[1]
    
    return rate, err, pre_exp, pre_exp_err


#################################################### Scaling factor Intensities ####################################################
# def ScalingFactor(SimulatedIntensity, IntensityRelaxometry, ErrIntensity):       #X is a vector contaning the parameters of the spectral density funciton
#     scaling = np.asarray([[0.0] for B in IntensityRelaxometry])
#     for B in range(len(IntensityRelaxometry)):
#         ScalingNum = 0.0
#         ScalingDen = 0.0
#         for VC in range(len(IntensityRelaxometry[B])):
#             if IntensityRelaxometry[B][VC][1] != "NA":
#                 ScalingNum +=  IntensityRelaxometry[B][VC][1]*SimulatedIntensity[B][VC]/(IntensityRelaxometry[B][VC][2])**2
#                 ScalingDen += (SimulatedIntensity[B][VC]/IntensityRelaxometry[B][VC][2])**2
#         scaling[B] = ScalingNum/ScalingDen
            
#     return scaling


#################################################### Chi2 HF rates ####################################################
def Chi2HF(X, HF_data, HF_err, tauc, OtherInput, AA):       #X is a vector contaning the parameters of the spectral density funciton
    Chi2Calc = 0.0
    
    for Func in HF_data.keys():
        for B0 in HF_data[Func].keys():
            try:
                diff = RelaxFunc[Func](B0, X, tauc, OtherInput) - HF_data[Func][B0][AA]
                err = HF_err[Func][B0][AA]
                Chi2Calc += (diff/err)**2
            except:
                pass
            
    return Chi2Calc

    
#################################################### Chi2 Intensities ####################################################
def Chi2I(SimulatedIntensity, scaling, Intensity, Err_Intensity, AA):
    Chi2Calc = 0.0

    for exp in Intensity.keys():          #experiment number dimension
        for vc in Intensity[exp][AA].keys():  #relaxation delay dimension
            if Intensity[exp][AA][vc] != "NA":
                diff = scaling[exp]*SimulatedIntensity[exp][vc] - Intensity[exp][AA][vc]
                err = Err_Intensity[exp][AA][vc]
                Chi2Calc += (diff/err)**2
            
    return Chi2Calc
            

#################################################### Chi2 Total ####################################################
def Chi2TOT(X, SimulatedIntensity, scaling, Intensity, Err_Intensity, HF_data, HF_err, tauc, OtherInput, AA):
    chi2I = Chi2I(SimulatedIntensity, scaling, Intensity, Err_Intensity, AA)
    chi2HF = Chi2HF(X, HF_data, HF_err, tauc, OtherInput, AA)
    
    return chi2I + chi2HF


#################################################### Back-calculated R1 at LF ####################################################
def CalcR1_LF(mcmc_param, tauc, other_input, B0_LF):
    back_calc_rates = {}
    
    for exp in B0_LF.keys():
        back_calc_rates[exp] = RelaxFunc['R1'](B0_LF[exp], mcmc_param, tauc, other_input)
        
    return back_calc_rates


#################################################### fit R1 at LF ####################################################
def FitR1_LF(intensities, AA):
    R1_fitted = {}
    
    for exp in intensities.keys():
        int_fit = []
        delays = []
        for vc in intensities[exp][AA].keys():
            if intensities[exp][AA][vc] != 'NA':
                delays.append(vc)
                int_fit.append(intensities[exp][AA][vc])
                
        r1, a, b, c = fit_intensity_decay(delays, int_fit)
        R1_fitted[exp] = r1
        
    return R1_fitted


#################################################### Fit of B0 ####################################################
def B0Fit(x, LowerCoefs, MiddleCoefs, HigherCoefs):
    Field = 0.0
    if x <= 0.47:
        for i in range(0, len(HigherCoefs)):
            Field += HigherCoefs[i] * math.pow(x, len(HigherCoefs) - i - 1)
    else:
        if x <= 0.6:
            for i in range(0, len(MiddleCoefs)):
                Field += MiddleCoefs[i] * math.pow(x, len(MiddleCoefs) - i - 1)
        else:
            for i in range(0, len(LowerCoefs)):
                Field += LowerCoefs[i] * math.pow(x, len(LowerCoefs) - i - 1)
    return Field
