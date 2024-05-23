import numpy as np

#################################################### Defined relaxation rates ####################################################

RelaxationRates = ["R1", "R2", "Sigma", "EtaZ", "EtaXY"]

#################################################### Variables ####################################################

PositionAuto = 0

Names = {'OrderParam': ['p1'], 'CorrTimes': ['Drot','lk21'], 'others': ['CSAfactor']}

def ImportFunc():
    from Rates import _R1calculation
    from Rates import _R2calculation
    from Rates import _Sigmacalculation
    from Rates import _EtaZcalculation
    from Rates import _EtaXYcalculation
    
    rates_func = {}
    rates_func['R1'] = _R1calculation.R1calculation
    rates_func['R2'] = _R2calculation.R2calculation
    rates_func['Sigma'] = _Sigmacalculation.Sigmacalculation
    rates_func['EtaZ'] = _EtaZcalculation.EtaZcalculation
    rates_func['EtaXY'] = _EtaXYcalculation.EtaXYcalculation
    
    return rates_func


def Cons(X):
    p1, Drot, lk21, CSAfactor, lnf = X
    if 1e10 < Drot < 1e11 and 0.75 < p1 < 1.0 and 8. < lk21 < 10. and 0. < CSAfactor < 2. and -10. < lnf < 1.:
        return 0.0
    else:
        return -np.inf


