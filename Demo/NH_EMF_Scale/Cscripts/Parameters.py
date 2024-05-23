import numpy as np

#################################################### Defined relaxation rates ####################################################

RelaxationRates = ["R1", "R2", "Sigma"]

#################################################### Variables ####################################################

PositionAuto = 0

Names = {'OrderParam': ["Sf2", "Ss2"], 'CorrTimes': ["tau_f", "tau_s"], 'others': []}

def ImportFunc():
    from Rates import _R1calculation
    from Rates import _R2calculation
    from Rates import _Sigmacalculation
    
    rates_func = {}
    rates_func['R1'] = _R1calculation.R1calculation
    rates_func['R2'] = _R2calculation.R2calculation
    rates_func['Sigma'] = _Sigmacalculation.Sigmacalculation
    
    return rates_func


def Cons(X):
    Sf2, Ss2, tau_f, tau_s, lnf = X
    if 0.0 < Sf2 < 1 and 0.0 < Ss2 < 1 and 1e-13 < tau_f < tau_s < 1e-8 and tau_f < 5e-11 and -10. < lnf < 1.:
        return 0.0
    else:
        return -np.inf


