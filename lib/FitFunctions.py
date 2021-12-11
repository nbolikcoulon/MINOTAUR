#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#########################################################
#                                                       #
#               Functions used during fit               #
#                                                       #
#########################################################

import math
import numpy as np
import scipy.stats as st

import Parameters as ParamFile
RelaxFunc = ParamFile.ImportFunc()


#################################################### Exponential decay ####################################################
def exp(x, a, b):
    return b*np.exp(-a*x)



#################################################### Chi2 HF rates ####################################################
def Chi2HF(X, RelaxData, B0, tauc, OtherInput):       #X is a vector contaning the parameters of the spectral density funciton
    
    
    Chi2Calc = 0.0
    
    for Func in range(len(RelaxData)):
        for B in range(len(RelaxData[Func])):
            diff = RelaxFunc[Func](B0[Func][B], X, tauc, OtherInput)[0] - RelaxData[Func][B][1]
            err = RelaxData[Func][B][2]
            Chi2Calc += (diff/err)**2
            
    return Chi2Calc



#################################################### Scaling factor Intensities ####################################################
def ScalingFactor(SimulatedIntensity, IntensityRelaxometry):       #X is a vector contaning the parameters of the spectral density funciton
    scaling = np.asarray([[0.0] for B in IntensityRelaxometry])
    for B in range(len(IntensityRelaxometry)):
        ScalingNum = 0.0
        ScalingDen = 0.0
        for VC in range(len(IntensityRelaxometry[B])):
            if IntensityRelaxometry[B][VC][1] != "NA":
                ScalingNum +=  IntensityRelaxometry[B][VC][1]*SimulatedIntensity[B][VC]/(IntensityRelaxometry[B][VC][2])**2
                ScalingDen += (SimulatedIntensity[B][VC]/IntensityRelaxometry[B][VC][2])**2
        scaling[B] = ScalingNum/ScalingDen
            
    return scaling
            
    
    
#################################################### Chi2 Intensities ####################################################
def Chi2I(SimulatedIntensity, IntensityRelaxometry):       #X is a vector contaning the parameters of the spectral density funciton
    Chi2Calc = 0.0
    scaling = [[0.0] for B in IntensityRelaxometry]

    for B in range(len(IntensityRelaxometry)):
        ScalingNum = 0.0
        ScalingDet = 0.0
        for VC in range(len(IntensityRelaxometry[B])):
            if IntensityRelaxometry[B][VC][1] != "NA":
                ScalingNum += IntensityRelaxometry[B][VC][1]*SimulatedIntensity[B][VC]/(IntensityRelaxometry[B][VC][2])**2
                ScalingDet += (SimulatedIntensity[B][VC])**2/(IntensityRelaxometry[B][VC][2])**2
        scaling[B] = ScalingNum/ScalingDet

    for B in range(len(IntensityRelaxometry)):          #experiment number dimension
        for VC in range(len(IntensityRelaxometry[B])):  #relaxation delay dimension
            if IntensityRelaxometry[B][VC][1] != "NA":
                diff = scaling[B]*SimulatedIntensity[B][VC] - IntensityRelaxometry[B][VC][1]
                err = IntensityRelaxometry[B][VC][2]
                Chi2Calc += (diff/err)**2
            
    return Chi2Calc
            



#################################################### Chi2 Total ####################################################
def Chi2TOT(X, SimulatedIntensity, IntensityRelaxometry, RelaxData, B0, tauc, OtherInput):
    chi2I = Chi2I(SimulatedIntensity, IntensityRelaxometry)
    chi2HF = Chi2HF(X, RelaxData, B0, tauc, OtherInput)
    
    return chi2I + chi2HF




#################################################### Proton - Deuterium position fit ####################################################    
def Chi2_forProton(X, param, tc, B0, B02, OtherInputs, R, E, RelaxData, IntensityRelaxometry, PropFunction, MagField, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, HzFunc):
    
    Chi2calc = 0.0
    
    for field in range(len(B0)):
        if R[field] != "NA":
            calc = HzFunc(X, param, tc, B0[field][0], OtherInputs)[0]
            meas = float(R[field])
            err = float(E[field])
        
            diff = calc - meas
                
            Chi2calc += (diff/err)**2
            
    Inputs = OtherInputs.append(X[0])
    Inputs = OtherInputs.append(X[1])
            
    Chi2calc += Chi2HF(param, RelaxData, B02, tc, Inputs)
    
    SimulatedIntensity = PropFunction(param, tc, OtherInputs, MagField, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, ParamFile.PositionAuto)
    Chi2calc += Chi2I(SimulatedIntensity, IntensityRelaxometry)

    return Chi2calc




def FitFunction_forProton(param, tc, B0, B02, OtherInputs, R, E, RelaxData, IntensityRelaxometry, PropFunction, MagField, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, HzFunc):
    x0 = [-1.0e-10, -1.0e-10]
    scaling = [1.0e-13, 1e-13]
    
    
    Nround = 1
    Niter = 10000
    print("  Round ", str(Nround))
    while True:
        
        naccept = 0.0
        
        for i in range(Niter):
            x0_p = [x0[j] + st.norm(0, 0.5).rvs()*scaling[j] for j in range(2)]
            if x0_p[1] < 0.0 and np.sqrt(x0_p[0]**2+x0_p[1]**2) > 1.0e-10:
                rho = Chi2_forProton(x0_p, param, tc, B0, B02, OtherInputs, R, E, RelaxData, IntensityRelaxometry, PropFunction, MagField, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, HzFunc)/Chi2_forProton(x0, param, tc, B0, B02, OtherInputs, R, E, RelaxData, IntensityRelaxometry, PropFunction, MagField, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, HzFunc)
                if rho < 1.0:
                    naccept += 1.0
                    x0 = x0_p
                    
                    
            
                           
#Convergence evaluation
        #Nothing happened
        if naccept == 0.0:
            Chi2calc = Chi2_forProton(x0, param, tc, B0, B02, OtherInputs, R, E, RelaxData, IntensityRelaxometry, PropFunction, MagField, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, HzFunc)
            print(" ")
            print(x0, Chi2calc)
            print(" ")
            return x0
        else:
            if Nround == 1:
                minRedChi2ToCompare = Chi2_forProton(x0, param, tc, B0, B02, OtherInputs, R, E, RelaxData, IntensityRelaxometry, PropFunction, MagField, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, HzFunc)
                Nround = 2
                Niter = 1000
                print("  Round ", str(Nround))

            else:
                NewMinRedChi2 = Chi2_forProton(x0, param, tc, B0, B02, OtherInputs, R, E, RelaxData, IntensityRelaxometry, PropFunction, MagField, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, HzFunc)
            
                #Very small variation of Chi2
                if abs(NewMinRedChi2-minRedChi2ToCompare)/minRedChi2ToCompare < 0.001:
                    print(" ")
                    print(x0, NewMinRedChi2)
                    print(" ")
                    return x0
                else:
                    minRedChi2ToCompare = NewMinRedChi2
                    Nround += 1
                    print("  Round ", str(Nround), " New Chi2: ", minRedChi2ToCompare)





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