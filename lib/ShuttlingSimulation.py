#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################################################
#                                                                          #
#   Here, we simulate the shuttling and give the intensity as the output   #
#      for one amino acids, all magnetic fields and relaxation delays      #
#                                                                          #
############################################################################

import numpy as np
from scipy import linalg

import FitFunctions as FitF

from _RelaxMat import RelaxMat as RM



def FieldList(self, Increment):
    if self.AccelerationType == "Constant Speed":
        SpeedUp = []
        SpeedDown = []
        for i in range(len(self.ExperimentNumber)):
            SpeedUp.append(self.Height[i] / self.SLF[i])
            SpeedDown.append(self.Height[i] / self.SHF[i])
        PositionListUp = [np.arange(0, self.Height[i], Increment * SpeedUp[i]) for i in range(len(self.ExperimentNumber))]     #a vector containing positions at different increment for each height
        for i in range(len(self.ExperimentNumber)):
            PositionListUp[i] = np.append(PositionListUp[i], self.Height[i])
        FieldListUp = [[FitF.B0Fit(PositionListUp[i][j], self.LowerCoefs, self.MiddleCoefs, self.HigherCoefs) for j in range(len(PositionListUp[i]))] for i in range(len(PositionListUp))]       #a vector containing the value of the field for the corresponding height
            
        PositionListDown = [np.arange(0, self.Height[i], Increment * SpeedDown[i]) for i in range(len(self.ExperimentNumber))]     #a vector containing positions at different increment for each height
        for i in range(len(self.ExperimentNumber)):
            PositionListDown[i] = np.append(PositionListDown[i], self.Height[i])
        PositionListDown = [PositionListDown[i][::-1] for i in range(len(PositionListDown))]
        FieldListDown = [[FitF.B0Fit(PositionListDown[i][j], self.LowerCoefs, self.MiddleCoefs, self.HigherCoefs) for j in range(len(PositionListDown[i]))] for i in range(len(PositionListDown))]       #a vector containing the value of the field for the corresponding height
        
    else:
        Acceleration = []
        for i in range(len(self.ExperimentNumber)):
            Acceleration.append(4.0*self.Height[i] / (self.SLF[i] * self.SLF[i]))
        TimeListUp = [np.arange(0, self.SLF[i], Increment) for i in range(len(self.ExperimentNumber))]
        TimeListDown = [np.arange(0, self.SHF[i], Increment) for i in range(len(self.ExperimentNumber))]
        def getPosition(t, T, A, H):         #t for time, T for time to get to target point, A for acceleraiton, H target point
            if t <= T/2.0:
                h = 1.0 / 2.0 * A * t**2
            else:
                h = -1.0 / 2.0 * A * t**2 + 4.0*H/T * t - H
            return h
    
        PositionListUp = []
        PositionListDown = []
        for i in range(len(self.ExperimentNumber)):
            IntermediateListUp = []
            for Time in TimeListUp[i]:
                IntermediateListUp.append(getPosition(Time, self.SLF[i], Acceleration[i], self.Height[i]))
            IntermediateListDown = []
            for Time in TimeListDown[i]:
                IntermediateListDown.append(getPosition(Time, self.SHF[i], Acceleration[i], self.Height[i]))
            PositionListUp.append(IntermediateListUp)
            PositionListDown.append(IntermediateListDown)
    
        FieldListUp = [[FitF.B0Fit(PositionListUp[i][j], self.LowerCoefs, self.MiddleCoefs, self.HigherCoefs) for j in range(len(PositionListUp[i]))] for i in range(len(PositionListUp))]
        PositionListDown = [PositionListDown[i][::-1] for i in range(len(PositionListDown))]
        FieldListDown = [[FitF.B0Fit(PositionListDown[i][j], self.LowerCoefs, self.MiddleCoefs, self.HigherCoefs) for j in range(len(PositionListDown[i]))] for i in range(len(PositionListDown))]
        
    return FieldListUp, FieldListDown





def PropCalDiag(param, tauc, OtherInputs, magfield, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, PosAuto):
    P = [[] for i in range(5)]      #Propagator
    
#HF
    RMHF = np.array(RM(magfield, param, tauc, OtherInputs))[0]
    
    eigHF, eigvecHF = np.linalg.eig(RMHF)
    for wthf in range(len(WTHF)):
        P[0].append(eigvecHF @ np.diag(np.exp(-WTHF[wthf] * eigHF)) @ np.linalg.inv(eigvecHF))
        P[4].append(eigvecHF @ np.diag(np.exp(-d25[wthf] * eigHF)) @ np.linalg.inv(eigvecHF))
        
    
#Shuttle HF to LF
    RMSHLF = [[[] for k in FieldListUp[FieldUp]] for FieldUp in range(len(FieldListUp))]
    for FieldUp in range(len(FieldListUp)):
        for k in range(len(FieldListUp[FieldUp])):
            RMSHLF[FieldUp][k] = np.array(RM(FieldListUp[FieldUp][k], param, tauc, OtherInputs))[0]
    
    P[1] = [[] for i in FieldListUp]
    rise = [[[] for k in FieldListUp[i]] for i in range(len(FieldListUp))]
    for FieldUp in range(len(FieldListUp)):
        for k in range(len(FieldListUp[FieldUp])):
            Eig, Vec = np.linalg.eig(RMSHLF[FieldUp][k])
            rise[FieldUp][k] = Vec @ np.diag(np.exp(- Increment * Eig)) @ np.linalg.inv(Vec)
        prod = rise[FieldUp][0]
        for k in range(1, len(FieldListUp[FieldUp])):
            prod = prod @ rise[FieldUp][k]
        P[1][FieldUp] = prod
        
#WT LF
    P[2] = [[] for i in LFtimes]
    for wtlf in range(len(B0LFields)):
        RMLF = np.array(RM(B0LFields[wtlf], param, tauc, OtherInputs))[0]
        Eig, Vec = np.linalg.eig(RMLF)
        
        for vc in range(len(LFtimes[wtlf])):
            P[2][wtlf].append(Vec @ np.diag(np.exp(-LFtimes[wtlf][vc] * Eig)) @ np.linalg.inv(Vec))

#Shuttle LF to HF
    RMSHF = [[[] for k in FieldListDown[FieldDown]] for FieldDown in range(len(FieldListDown))]
    for FieldDown in range(len(FieldListDown)):
        for k in range(len(FieldListDown[FieldDown])):
            RMSHF[FieldDown][k] = np.array(RM(FieldListDown[FieldDown][k], param, tauc, OtherInputs))[0]
            
    P[3] = [[] for i in FieldListDown]
    down = [[[] for k in FieldListDown[i]] for i in range(len(FieldListDown))]
    for FieldDown in range(len(FieldListDown)):
        for k in range(len(FieldListDown[FieldDown])):
            Eig, Vec = np.linalg.eig(RMSHF[FieldDown][k])
            down[FieldDown][k] = Vec @ np.diag(np.exp(-Increment * Eig)) @ np.linalg.inv(Vec)
        prod = down[FieldDown][0]
        for k in range(1, len(FieldListDown[FieldDown])):
            prod = prod @ down[FieldDown][k]
        P[3][FieldDown] = prod

            
#Total Propagation
    TotalExpVal = [[] for i in ExperimentNumber]
                    
    for ExpNum in range(len(ExperimentNumber)):
        PropInt1 = P[4][ExpNum] @ P[3][ExpNum]
        PropInt2 = P[1][ExpNum] @ P[0][ExpNum]
        for vc in range(len(LFtimes[ExpNum])):
            TotalProp = PropInt1 @ P[2][ExpNum][vc] @ PropInt2
            TotalExpVal[ExpNum].append(TotalProp[PosAuto][PosAuto])
            
           
    return TotalExpVal




def PropCalExp(param, tauc, OtherInputs, magfield, Increment, FieldListUp, FieldListDown, ExperimentNumber, WTHF, d25, LFtimes, B0LFields, PosAuto):
    P = [[] for i in range(5)]      #Propagator
    
#HF
    RMHF = np.array(RM(magfield, param, tauc, OtherInputs))[0]
    for wthf in range(len(WTHF)):
        P[0].append(linalg.expm(-WTHF[wthf] * RMHF))
        P[4].append(linalg.expm(-d25[wthf] * RMHF))
        
    
#Shuttle HF to LF
    RMSHLF = [[[] for k in FieldListUp[FieldUp]] for FieldUp in range(len(FieldListUp))]
    for FieldUp in range(len(FieldListUp)):
        for k in range(len(FieldListUp[FieldUp])):
            RMSHLF[FieldUp][k] = np.array(RM(FieldListUp[FieldUp][k], param, tauc, OtherInputs))[0]
    
    P[1] = [[] for i in FieldListUp]
    rise = [[[] for k in FieldListUp[i]] for i in range(len(FieldListUp))]
    for FieldUp in range(len(FieldListUp)):
        for k in range(len(FieldListUp[FieldUp])):
            rise[FieldUp][k] = linalg.expm(- Increment * RMSHLF[FieldUp][k])
        prod = rise[FieldUp][0]
        for k in range(1, len(FieldListUp[FieldUp])):
            prod = prod @ rise[FieldUp][k]
        P[1][FieldUp] = prod
        
#WT LF
    P[2] = [[] for i in LFtimes]
    for wtlf in range(len(B0LFields)):
        RMLF = np.array(RM(B0LFields[wtlf], param, tauc, OtherInputs))[0]
        
        for vc in range(len(LFtimes[wtlf])):
            P[2][wtlf].append(linalg.expm(-LFtimes[wtlf][vc] * RMLF))

#Shuttle LF to HF
    RMSHF = [[[] for k in FieldListDown[FieldDown]] for FieldDown in range(len(FieldListDown))]
    for FieldDown in range(len(FieldListDown)):
        for k in range(len(FieldListDown[FieldDown])):
            RMSHF[FieldDown][k] = np.array(RM(FieldListDown[FieldDown][k], param, tauc, OtherInputs))[0]
            
    P[3] = [[] for i in FieldListDown]
    down = [[[] for k in FieldListDown[i]] for i in range(len(FieldListDown))]
    for FieldDown in range(len(FieldListDown)):
        for k in range(len(FieldListDown[FieldDown])):
            down[FieldDown][k] = linalg.expm(- Increment * RMSHF[FieldDown][k])
        prod = down[FieldDown][0]
        for k in range(1, len(FieldListDown[FieldDown])):
            prod = prod @ down[FieldDown][k]
        P[3][FieldDown] = prod

            
#Total Propagation
    TotalExpVal = [[] for i in ExperimentNumber]
                    
    for ExpNum in range(len(ExperimentNumber)):
        PropInt1 = P[4][ExpNum] @ P[3][ExpNum]
        PropInt2 = P[1][ExpNum] @ P[0][ExpNum]
        for vc in range(len(LFtimes[ExpNum])):
            TotalProp = PropInt1 @ P[2][ExpNum][vc] @ PropInt2
            TotalExpVal[ExpNum].append(TotalProp[PosAuto][PosAuto])
            
            
    return TotalExpVal




def optShuttling(self, ExpNum, ParamOpt, PosAuto):
    Increment = 1e-3
    FieldListUp_Init, FieldListDown_Init = FieldList(self, Increment)
    
    posMaxVC = self.VC[ExpNum].index(max(self.VC[ExpNum]))
    
    ExpVal_Init = []
    for i in range(len(ParamOpt)):
        ExpVal_InitF = PropCalDiag(ParamOpt[i], self.TauC, np.array(self.OtherInputs[0][1]), self.MagField, Increment, FieldListUp_Init, FieldListDown_Init, [self.ExperimentNumber[ExpNum]], [self.WTHF[ExpNum]], [self.d25[ExpNum]], [self.LFtimes[ExpNum]], [self.B0LFields[0][ExpNum]], PosAuto)
        ExpVal_Init.append(ExpVal_InitF[-1][posMaxVC])
    print(" Initial value : ", Increment, "s")
        
    Success = 0
    count = 0
    while Success == 0:
        count += 1
        Increment = Increment + 1e-3
        FieldListUp, FieldListDown = FieldList(self, Increment)
        
        ExpVal = []
        for i in range(len(ParamOpt)):
            ExpVal_F = PropCalDiag(ParamOpt[i], self.TauC, np.array(self.OtherInputs[0][1]), self.MagField, Increment, FieldListUp, FieldListDown, [self.ExperimentNumber[ExpNum]], [self.WTHF[ExpNum]], [self.d25[ExpNum]], [self.LFtimes[ExpNum]], [self.B0LFields[0][ExpNum]], PosAuto)
            ExpVal.append(ExpVal_F[-1][posMaxVC])
        
        for i in range(len(ParamOpt)):
            if abs(ExpVal[i]-ExpVal_Init[i])/ExpVal_Init[i] > 0.01:
                Success = 1
                break
        if Success ==0:
            print(" Updated value ", count, " : ", Increment, " s")
        else:
            Increment = Increment - 1e-3
                

    return Increment