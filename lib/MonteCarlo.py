#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################################
#                                                    #
#               Monte Carlo simulation               #
#                                                    #
######################################################
from random import uniform
from multiprocessing import Pool

from lib import FigOutputs as FigOut

import numpy as np
import emcee
from math import exp, log


import Parameters as ParamFile




#################################################### Scaling factor Intensities ####################################################
def ScalingFactor(SimulatedIntensity, IntensityRelaxometry):
    scaling = np.asarray([[0.0] for B in IntensityRelaxometry])
    for B in range(len(IntensityRelaxometry)):
        ScalingNum = 0.0
        ScalingDen = 0.0
        for VC in range(len(IntensityRelaxometry[B])):
            if IntensityRelaxometry[B][VC][1] != "NA":
                ScalingNum += IntensityRelaxometry[B][VC][1]*SimulatedIntensity[B][VC]/(IntensityRelaxometry[B][VC][2])**2
                ScalingDen += (SimulatedIntensity[B][VC]/IntensityRelaxometry[B][VC][2])**2
        scaling[B] = ScalingNum/ScalingDen
            
    return scaling


#################################################### lnlike ####################################################
def lnlike(X, IntensityRelaxometry, RelaxData, B0HF, tauc, OtherInputs, magfield, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, PropFunction):       #X is a vector contaning the parameters of the spectral density funciton
    
    RelaxFunc = ParamFile.ImportFunc()
    
    distCalc = 0.0
    
    expLnF = exp(2.0*X[-1])
    
    for Func in range(len(RelaxData)):
        for B in range(len(B0HF[Func])):
            model_real = RelaxFunc[Func](B0HF[Func][B], X[:-1], tauc, OtherInputs)[0]
            model = model_real
            
            err = RelaxData[Func][B][2]

            inv_sigma2 = err**2 + model**2 * expLnF
            dist = -1.0/2.0 * ((RelaxData[Func][B][1] - model)**2 / inv_sigma2 + log(2.0*np.pi*inv_sigma2))
                
            distCalc += dist
            

    SimulatedInt = PropFunction(X[:-1], tauc, OtherInputs, magfield, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, ParamFile.PositionAuto)
    

    #scaling factor
    scaling = ScalingFactor(SimulatedInt, IntensityRelaxometry)
            
    for B in range(len(IntensityRelaxometry)):
        Nvc = len(list(filter(lambda x: x!='NA', np.asarray(IntensityRelaxometry[B])[:,1])))
        for VC in range(len(IntensityRelaxometry[B])):
            if IntensityRelaxometry[B][VC][1] != "NA":
                model_real = scaling[B]*SimulatedInt[B][VC]
                model = model_real
                
                err = IntensityRelaxometry[B][VC][2]
    
                inv_sigma2 = err**2 + model**2 * expLnF
                dist = -1.0/2.0 * ((IntensityRelaxometry[B][VC][1] - model)**2 / inv_sigma2 + log(2.0*np.pi*inv_sigma2))
                    
                distCalc += dist/Nvc
            
                
    return distCalc



    
    
#################################################### lnprob ####################################################
def lnprob(X, IntensityRelaxometry, RelaxData, B0HF, tauc, OtherInputs, magfield, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, PropFunction):
    lp = ParamFile.Cons(X)
    if lp != 0.0:
        return -np.inf
    else:
        return lp + lnlike(X, IntensityRelaxometry, RelaxData, B0HF, tauc, OtherInputs, magfield, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, PropFunction)
    


#################################################### nWalker check ####################################################
def nWalkerCheck(nWalkers_b, nParam):
    if nWalkers_b < 2*nParam:
        nWalkers = 2*nParam
        print("")
        print("WARNING")
        print("The number of chain is changed to *" + str(nWalkers) + "* because it must be at least 2-times higher than the number of parameters")
        print("")
    else:
        if nWalkers_b % 2 == 0:
            nWalkers = nWalkers_b
        else:
            nWalkers = nWalkers_b + 1
            print("")
            print("WARNING")
            print("The number of chain is changed to *" + str(nWalkers) + "* because it must be even")
            print("")
        
    return nWalkers


#################################################### MCMC ####################################################
def MarkovChainMonteCarlo(CorrDirName, TrajDirName, IntensityRelaxometry, RelaxData, B0HF, tauc, OtherInputs, magfield, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, nWalkers, nMcmc, Bonds, labels, AA, PropFunction, nParam):
    
    Labels = [[] for p in range(nParam)]
    for i in range(len(labels)):
        Labels[i] = labels[i]
    Labels[-1] = "f"

    Bonds.append([-10.0, 1.0])
    
    lim = min(int(0.2*nMcmc), 500)

    pos = []
    for w in range(nWalkers):
        toadd = []
        for i in range(nParam):
            val = uniform(Bonds[i][0], Bonds[i][1])
            toadd.append(val)
        pos.append(toadd)
    pos = np.array(pos)
    
    
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nWalkers, nParam, lnprob, args=(IntensityRelaxometry, RelaxData, B0HF, tauc, OtherInputs, magfield, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, PropFunction), pool = pool)

        sampler.run_mcmc(pos, nMcmc, progress=True);
    
    samples = sampler.chain[:, lim:, :].reshape((-1, nParam))
    FullSamples = sampler.chain[:, :, :].reshape((-1, nParam))
            
            
    Mean = np.percentile(samples, 50, axis=0)
    Sixteenth = np.percentile(samples, 16, axis=0)
    Heigthyfourth = np.percentile(samples, 84, axis=0)
            
    MAF = np.mean(sampler.acceptance_fraction)
            
    PosError = []
    NegError = []
    for p in range(nParam):
        PosError.append(Heigthyfourth[p] - Mean[p])
        NegError.append(Mean[p] - Sixteenth[p])
    
    print("")
    for P in range(nParam):
        print(Labels[P], Mean[P], PosError[P], NegError[P])
    print("")
    print("Mean Acceptance Fracttion: ", MAF)
    

    FigOut.FigTraj(TrajDirName, sampler, Mean, Labels, lim, AA)
    FigOut.FigCorr(CorrDirName, samples, Mean, Labels, AA)
            
    return [Mean, PosError, NegError], MAF, FullSamples

