import numpy as np

#################################################### Defined relaxation rates ####################################################

RelaxationRates = ["R2", "R1", "Sigma"]

#################################################### Variables ####################################################

PositionAuto = 0

Names = {'OrderParam': ['Sf2','Ss2'], 'CorrTimes': ['tm','tf','ts'], 'others': ['CSApC']}

def ImportFunc():
    import _R2calculation
    import _R1calculation
    import _Sigmacalculation
    
    
    return [_R2calculation.R2calculation, _R1calculation.R1calculation, _Sigmacalculation.Sigmacalculation]





def Cons(X):
    Sf2, Ss2, tm, tf, ts, CSApC, lnf = X
    if 0.0 < Sf2 < 1.0 and 0.0 < Ss2 < 1.0 and 1e-13 < tm < tf < ts < 2.0e-8 and 0.0 < CSApC < 35.0 and -10. < lnf < 1.:
        return 0.0
    else:
        return -np.inf