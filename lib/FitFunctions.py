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
relaxation_function = ParamFile.ImportFunc()


def exp(x, a, b):
    """
    exponential decay with a pre-exponential factor

    Parameters
    ----------
    x : TYPE: float
        DESCRIPTION: entry variable.
    a : TYPE: float
        DESCRIPTION: exponential decay rate.
    b : TYPE: float
        DESCRIPTION: pre-exponential factor.

    Returns
    -------
    TYPE: float
        DESCRIPTION: b*e^{-a*x}.

    """
    return b*np.exp(-a*x)


def residuals(variable, data, function, *args):
    """
    computes residuals for a given fitted function.

    Parameters
    ----------
    variable : TYPE: array
        DESCRIPTION: variable over which residuals are computed (typically field).
    data : TYPE: array
        DESCRIPTION: experimental data.
    function : TYPE: python function
        DESCRIPTION: function to back-calculate experimental data.
    *args : TYPE
        DESCRIPTION: list of arguments for function.

    Returns
    -------
    residuals : TYPE: array
        DESCRIPTION: array of residuals.

    """
    residuals = np.empty(shape=(len(variable)), dtype=float)
    
    for c, val in enumerate(variable):
        residuals[c] = data[c] - function(val, *args)
        
    return residuals


def fit_intensity_decay(time, data):
    """
    exponential fitting of the experimental intensity decays. Uses curve_fit.

    Parameters
    ----------
    time : TYPE: array
        DESCRIPTION: relaxation delays.
    data : TYPE: array
        DESCRIPTION: intensities.

    Returns
    -------
    rate : TYPE: float
        DESCRIPTION: decay rate.
    rate_err : TYPE: float
        DESCRIPTION: error on the decay rate.
    pre_exp : TYPE: float
        DESCRIPTION: pre-exponential factor.
    pre_exp_err : TYPE: float
        DESCRIPTION: error on the pre-exponential factor.

    """
    param, cov = curve_fit(exp, time, data)
    
    rate = param[0]
    pre_exp = param[1]
    
    err_matrix = np.sqrt(np.diag(cov))
    rate_err = err_matrix[0]
    pre_exp_err = err_matrix[1]
    
    return rate, rate_err, pre_exp, pre_exp_err


def Chi2_HF(X, HF_data, tauc, other_input, AA):
    """
    Calculates the chi2 of the high-field rates

    Parameters
    ----------
    X : TYPE: array
        DESCRIPTION: parameters determined during MCMC.
    HF_data : TYPE: dictionnary
        DESCRIPTION: high-field data.
    tauc : TYPE: float
        DESCRIPTION: value of tau_c.
    other_input : TYPE: dictionnary
        DESCRIPTION: other inputs of the relaxation rates functions.
    AA : TYPE: str
        DESCRIPTION: residue.

    Returns
    -------
    Chi2Calc : TYPE: float
        DESCRIPTION: value of the chi2 for high-field data.

    """
    Chi2Calc = 0.0
    
    for Func in HF_data.keys():
        for B0 in HF_data[Func].keys():
            try:
                diff = relaxation_function[Func](B0, X, tauc, other_input[AA]) - HF_data[Func][B0]['data'][AA]
                err = HF_data[Func][B0]['error'][AA]
                Chi2Calc += (diff / err)**2
            except:
                pass
            
    return Chi2Calc

    
def Chi2_Intensities(simulated_intensity, scaling, experimental_intensity, AA):
    """
    Calculates the chi2 of the relaxometry intensity decays

    Parameters
    ----------
    simulated_intensity : TYPE: dictionnary
        DESCRIPTION: simulated intensities.
    scaling : TYPE: dictionnary
        DESCRIPTION: scaling factors for the simulated intensities.
    experimental_intensity : TYPE: dictionnary
        DESCRIPTION: experimental intensity decays.
    AA : TYPE: str
        DESCRIPTION: residue.

    Returns
    -------
    Chi2Calc : TYPE: float
        DESCRIPTION: value of the chi2 for relaxometry intensity decays.

    """
    Chi2Calc = 0.0

    for exp in experimental_intensity.keys():
        for vc in experimental_intensity[exp][AA].keys():
            if experimental_intensity[exp][AA][vc]['data'] != "NA":
                diff = scaling[AA][exp]*simulated_intensity[AA][exp][vc] - experimental_intensity[exp][AA][vc]['data']
                err = experimental_intensity[exp][AA][vc]['error']
                Chi2Calc += (diff / err)**2
            
    return Chi2Calc
            

def Chi2_TOT(X, simulated_intensity, scaling, experimental_intensity, HF_data, tauc, other_input, AA):
    """
    Calculates the chi2 after MCMC.

    Parameters
    ----------
    X : TYPE: array
        DESCRIPTION: parameters determined during MCMC.
    SimulatedIntensity : TYPE: dictionnary
        DESCRIPTION: simulated intensities.
    scaling : TYPE: dictionnary
        DESCRIPTION: scaling factors for the simulated intensities.
    experimental_intensity : TYPE: dictionnary
        DESCRIPTION: experimental intensity decays.
    HF_data : TYPE: dictionnary
        DESCRIPTION: high-field data.
    tauc : TYPE: float
        DESCRIPTION: value of tau_c.
    other_input : TYPE: dictionnary
        DESCRIPTION: other inputs of the relaxation rates functions.
    AA : TYPE: str
        DESCRIPTION: residue.

    Returns
    -------
    TYPE: float
        DESCRIPTION: value of the chi2.

    """
    chi2I = Chi2_Intensities(simulated_intensity, scaling, experimental_intensity, AA)
    chi2HF = Chi2_HF(X, HF_data, tauc, other_input, AA)
    
    return chi2I + chi2HF


def Calc_R1_LF(mcmc_param, tauc, other_input, B0_LF):
    """
    Computes the R1 at the fields where relaxometry decays are recorded.

    Parameters
    ----------
    mcmc_param : TYPE: array
        DESCRIPTION: parameters determined during MCMC.
    tauc : TYPE: float
        DESCRIPTION: value of tau_c.
    other_input : TYPE: dictionnary
        DESCRIPTION: other inputs of the relaxation rates functions.
    B0_LF : TYPE: dictionnary
        DESCRIPTION: fields where relaxometry decays are recorded.

    Returns
    -------
    back_calc_rates : TYPE: dictionnary
        DESCRIPTION: low-field R1.

    """
    back_calc_rates = {}
    for exp in B0_LF.keys():
        back_calc_rates[exp] = relaxation_function['R1'](B0_LF[exp], mcmc_param, tauc, other_input)
        
    return back_calc_rates


def Fit_R1_LF(intensities, AA):
    """
    fit the experimental intensity decays to a single-exponential decay
    to get an apparent R1

    Parameters
    ----------
    intensities : TYPE: dictionnary
        DESCRIPTION: intensities decays.
    AA : TYPE: str
        DESCRIPTION: residue.

    Returns
    -------
    R1_fitted : TYPE: dictionnary
        DESCRIPTION: fitted R1 using single-exponential decay fit.

    """
    R1_fitted = {}
    for exp in intensities.keys():
        int_for_fit = []
        delays = []
        for vc in intensities[exp][AA].keys():
            if intensities[exp][AA][vc]['data'] != 'NA':
                delays.append(vc)
                int_for_fit.append(intensities[exp][AA][vc]['data'])
                
        r1, _, _, _ = fit_intensity_decay(delays, int_for_fit)
        R1_fitted[exp] = r1
        
    return R1_fitted


def detect_magnetic_tunnel(field_cal):
    """
    detects a potential magnetic tunnel

    Parameters
    ----------
    field_cal : TYPE: dict
        DESCRIPTION: matches height with measured magnetic field.

    Returns
    -------
    tunnel_position : TYPE: float
        DESCRIPTION: position at which the tunnel starts.
    tunnel_field : TYPE: float
        DESCRIPTION: field of the tunnel.

    """
    h_0, b_0 = list(field_cal.items())[1]
    for h, b in list(field_cal.items())[2:]:
        if b == b_0 and b < 4.0:
            tunnel_position = h_0
            tunnel_field = b_0
            return tunnel_position, tunnel_field
        h_0, b_0 = h, b
        
    tunnel_position = max(list(field_cal.keys()))
    tunnel_field = field_cal[tunnel_position]
    return tunnel_position, tunnel_field


def calibrate_B0(field_cal):
    """
    fitting of the field within the magnet.

    Parameters
    ----------
    field_cal : TYPE: dictionnary
        DESCRIPTION: matches height with measured magnetic field.

    Returns
    -------
    list_coeff : TYPE: list
        DESCRIPTION: contains coefficients of the polynomial fits.
    tunnel_position: TYPE: float
        DESCRIPTION: position at which the tunnel starts.
    tunnel_field : TYPE: float
        DESCRIPTION: field of the tunnel.

    """
    tunnel_position, tunnel_field = detect_magnetic_tunnel(field_cal)
    
    field_cal_high, field_cal_low = {}, {}
    for h in field_cal.keys():
        if h < tunnel_position:
            if h < 0.4:
                field_cal_high[h] = field_cal[h]
            else:
                field_cal_low[h] = field_cal[h]

    HigherCoefs = np.polyfit(list(field_cal_high.keys()), list(field_cal_high.values()), 15)
    if len(field_cal_low) != 0:
        LowerCoefs = np.polyfit(list(field_cal_low.keys()), list(field_cal_low.values()), 15)
    else:
        LowerCoefs = [0.0]
    
    list_coeff = {'high_fields': HigherCoefs,
                  'low_fields': LowerCoefs}

    return list_coeff, tunnel_position, tunnel_field
    
    
def Calc_B0(height, list_coeffs, tunnel_position, tunnel_field):
    """
    calculates B0 based on calibration results

    Parameters
    ----------
    height : TYPE: float
        DESCRIPTION: height in the spectometer at which the field is calculated.
    list_coeffs : TYPE: list
        DESCRIPTION: results from the field calibration.
    tunnel_position: TYPE: float
        DESCRIPTION: position at which the tunnel starts.
    tunnel_field : TYPE: float
        DESCRIPTION: field of the tunnel.

    Returns
    -------
    field : TYPE: float
        DESCRIPTION: field at height.

    """
    if height >= tunnel_position:
        return tunnel_field
    
    field = 0.
    if height < 0.4:
        deg = len(list_coeffs['high_fields'])
        for c, coeff in enumerate(list_coeffs['high_fields']):
            field += coeff * math.pow(height, deg - c - 1)
        return field
    
    deg = len(list_coeffs['low_fields'])
    for c, coeff in enumerate(list_coeffs['low_fields']):
        field += coeff * math.pow(height, deg - c - 1)
    return field


def Get_B0_low_field(setup, list_coeffs, tunnel_position, tunnel_field):
    """
    Calculates the low field where relaxation occurs

    Parameters
    ----------
    setup : TYPE: dictionnary
        DESCRIPTION: setup for the relaxometry experiments.
    list_coeffs : TYPE: dictionnary
        DESCRIPTION: results from the field calibration.
    tunnel_position: TYPE: float
        DESCRIPTION: position at which the tunnel starts.
    tunnel_field : TYPE: float
        DESCRIPTION: field of the tunnel.

    Returns
    -------
    B0_low_field : TYPE: dictionnary
        DESCRIPTION: contains the low fields where relaxation rates are recorded.

    """
    B0_low_field = {}
    for exp in setup.keys():
        if 'field' in setup[exp].keys():
            B0_low_field[exp] = setup[exp]['field']
        else:
            B0_low_field[exp] = Calc_B0(setup[exp]['height'], list_coeffs, tunnel_position, tunnel_field)
            
    return B0_low_field