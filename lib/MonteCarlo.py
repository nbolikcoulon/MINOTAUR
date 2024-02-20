#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################################
#                                                    #
#               Monte Carlo simulation               #
#                                                    #
######################################################
from random import uniform
from multiprocessing import Pool

# from lib import FigOutputs as FigOut
import FigOutputs as FigOut

import numpy as np
import emcee
from math import exp, log
import sys

import Parameters as ParamFile


#################################################### Scaling factor Intensities ####################################################
def ScalingFactor(SimulatedIntensity, IntensityRelaxometry, ErrRelaxometry, AA):
    scaling = {}
    
    for exp_num in IntensityRelaxometry.keys():
        ScalingNum = 0.0
        ScalingDen = 0.0
        for vc in IntensityRelaxometry[exp_num][AA].keys():
            if IntensityRelaxometry[exp_num][AA][vc] != "NA":
                ScalingNum += IntensityRelaxometry[exp_num][AA][vc]*SimulatedIntensity[exp_num][vc]/(ErrRelaxometry[exp_num][AA][vc])**2
                ScalingDen += (SimulatedIntensity[exp_num][vc]/ErrRelaxometry[exp_num][AA][vc])**2
        scaling[exp_num] = ScalingNum/ScalingDen
            
    return scaling


#################################################### lnlike ####################################################
def lnlike(X, AA, IntensityRelaxometry, ErrRelaxometry, RelaxData, RelaxErr, tauc, OtherInputs, magfield, Increment, FieldList, Delays, SetUp, B0LFields, PropFunction):       #X is a vector contaning the parameters of the spectral density funciton
    
    RelaxFunc = ParamFile.ImportFunc()
    
    distCalc = 0.0
    
    expLnF = exp(2.0*X[-1])
    
    for Func in RelaxData.keys():
        for B in RelaxData[Func].keys():
            try:
                exp_data = RelaxData[Func][B][AA]
                exp_err = RelaxErr[Func][B][AA]
                model = RelaxFunc[Func](B, X[:-1], tauc, OtherInputs[AA])
                
                inv_sigma2 = exp_err**2 + model**2 * expLnF
                dist = -1.0/2.0 * ((exp_data - model)**2 / inv_sigma2 + log(2.0*np.pi*inv_sigma2))

                distCalc += dist
            except:
                pass

    SimulatedInt = PropFunction(X[:-1], tauc, OtherInputs[AA], magfield, Increment, FieldList, Delays, SetUp, B0LFields, ParamFile.PositionAuto)
    scaling = ScalingFactor(SimulatedInt, IntensityRelaxometry, ErrRelaxometry, AA)
    
    for exp_num in SetUp.keys():
        N_vc = len(list(filter(lambda x: x!='NA', np.asarray(list(IntensityRelaxometry[exp_num].values())))))
        for vc in SetUp[exp_num]['vc']:
            if IntensityRelaxometry[exp_num][AA][vc] != 'NA':
                exp_data = IntensityRelaxometry[exp_num][AA][vc]
                exp_err = ErrRelaxometry[exp_num][AA][vc]
                model = scaling[exp_num]*SimulatedInt[exp_num][vc]
            
                inv_sigma2 = exp_err**2 + model**2 * expLnF
                dist = -1.0/2.0 * ((exp_data - model)**2 / inv_sigma2 + log(2.0*np.pi*inv_sigma2))
                    
                distCalc += dist / N_vc
                
    return distCalc

    
#################################################### lnprob ####################################################
def lnprob(X, AA, IntensityRelaxometry, ErrRelaxometry, RelaxData, RelaxErr, tauc, OtherInputs, magfield, Increment, FieldList, Delays, SetUp, B0LFields, PropFunction):
    lp = ParamFile.Cons(X)
    if lp != 0.0:
        return -np.inf
    else:
        return lp + lnlike(X, AA, IntensityRelaxometry, ErrRelaxometry, RelaxData, RelaxErr, tauc, OtherInputs, magfield, Increment, FieldList, Delays, SetUp, B0LFields, PropFunction)


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
def MarkovChainMonteCarlo(pdf_traj, pdf_corr, AA, nParam, CorrDirName, TrajDirName, tauc, magfield, nWalkers, nMcmc, Bonds, labels,
                          Increment, PropFunction, shuttling_fields, shuttling_delays, 
                          IntensityRelaxometry, ErrRelaxometry, RelaxData, RelaxErr, OtherInputs, set_up, B0LFields):
    
    Labels = []
    for l in labels:
        Labels.append(l)
    Labels.append("f")

    Bonds.append([-10.0, 1.0])
    
    cutoff = min(int(0.2*nMcmc), 500)

    pos = []
    for w in range(nWalkers):
        toadd = []
        for i in range(nParam):
            val = uniform(Bonds[i][0], Bonds[i][1])
            toadd.append(val)
        pos.append(toadd)
    pos = np.array(pos)
    
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nWalkers, nParam, lnprob,
                                        args=(AA, IntensityRelaxometry, ErrRelaxometry, RelaxData, RelaxErr, tauc, OtherInputs, magfield, Increment,
                                              shuttling_fields, shuttling_delays, set_up, B0LFields, PropFunction),
                                        pool = pool)
        try:
            sampler.run_mcmc(pos, nMcmc, progress=True);
        except:
            try:
                width = 30
                for i, result in enumerate(sampler.sample(pos, iterations=nMcmc)):
                    n = int((width+1) * float(i) / nMcmc)
                    sys.stdout.write("\r[{0}{1}]".format('#' * n, ' ' * (width - n)))
                sys.stdout.write("\n")
            except:
                sampler.run_mcmc(pos, nMcmc);
    
    samples = sampler.chain[:, cutoff:, :].reshape((-1, nParam))
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
        print(Labels[P], '\t', '{:.3e}'.format(Mean[P]), '\t', '{:.3e}'.format(PosError[P]), '\t', '{:.3e}'.format(NegError[P]))
    print("")
    print(f"Mean Acceptance Fracttion: {round(MAF, 4)}")
    

    FigOut.FigTraj(pdf_traj, TrajDirName, sampler, Mean, Labels, cutoff, AA)
    FigOut.FigCorr(pdf_corr, CorrDirName, samples, Mean, Labels, AA)
            
    return [Mean, PosError, NegError], MAF, FullSamples

