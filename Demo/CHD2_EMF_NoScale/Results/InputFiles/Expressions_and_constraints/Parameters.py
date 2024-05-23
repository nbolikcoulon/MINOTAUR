import numpy as np

#################################################### Defined relaxation rates ####################################################

RelaxationRates = ["R2", "R1", "Sigma", "R1H"]

#################################################### Variables ####################################################

PositionAuto = 0
pA = 1.0
Anisotropy1 = 'NO'

Names = {'OrderParam': ['Sf2','Ss2'], 'CorrTimes': ['tm','tf','ts'], 'others': ['CSApC']}

def ImportFunc():
    import _R2calculation
    import _R1calculation
    import _Sigmacalculation
    import _R1Hcalculation
    
    rates_func = {}
    rates_func['R1'] = _R1calculation.R1calculation
    rates_func['R2'] = _R2calculation.R2calculation
    rates_func['Sigma'] = _Sigmacalculation.Sigmacalculation
    rates_func['R1H'] = _R1Hcalculation.R1Hcalculation
    
    return rates_func





def Cons(X):
    Sf2, Ss2, tm, tf, ts, CSApC, lnf = X
    if 0.0 < Sf2 < 1.0 and 0.0 < Ss2 < 1.0 and 0.0 < tm < tf < ts < 30.0e-9 and 0.0 < CSApC < 35.0 and -10. < lnf < 1.:
        return 0.0
    else:
        return -np.inf