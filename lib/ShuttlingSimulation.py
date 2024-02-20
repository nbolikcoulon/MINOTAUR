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


def FieldList(self, increment):
    ### 
    #   make a single field list for all the experiments:
    #   consider the highest height the sample is shuttled to;
    #   compute the field at every optimized increment;
    #   create the delay time list for each experiments 
    #   use the slowest speed to compute these delays
    ###
    all_heights = []
    for exp in self.set_up.keys():
        all_heights.append(self.set_up[exp]['height'])
    max_height = max(all_heights)
    
    positions =  np.arange(increment, max_height, increment)
    for exp in self.set_up.keys():
        positions = np.append(positions, self.set_up[exp]['height'])
    positions = np.sort(positions)
    
    field_list = [FitF.B0Fit(p, self.LowerCoefs, self.MiddleCoefs, self.HigherCoefs) for p in positions]
    
    delays = {}
    delays['up'] = {}
    delays['down'] = {}
    if self.Shuttling_type == "Constant Speed":
        for exp in self.set_up.keys():
            speed_up = self.set_up[exp]['height'] / self.set_up[exp]['SLF']
            speed_down = self.set_up[exp]['height'] / self.set_up[exp]['SHF']
            
            delays['up'][exp] = []
            delays['up'][exp].append(positions[0] / speed_up)
            delays['down'][exp] = []
            delays['down'][exp].append(positions[0] / speed_down)
            for counter, p in enumerate(positions[1:]):
                if p > self.set_up[exp]['height']:
                    break
                else:
                    delays['up'][exp].append( (positions[counter+1] - positions[counter]) / speed_up )
                    delays['down'][exp].append( (positions[counter+1] - positions[counter]) / speed_down )

    else:
        for exp in self.set_up.keys():
            acceleration_up = 4.0*self.set_up[exp]['height'] / self.set_up[exp]['SLF'] * self.set_up[exp]['SLF']
            acceleration_down = 4.0*self.set_up[exp]['height'] / self.set_up[exp]['SHF'] * self.set_up[exp]['SHF']

            delays['up'][exp] = []
            delays['down'][exp] = []
            for p in positions:
                if p > self.set_up[exp]['height']:
                    break
                else:
                    if p < self.set_up[exp]['height'] / 2.:
                        dt_up = np.sqrt(2. * p / acceleration_up)
                        dt_down = np.sqrt(2. * p / acceleration_up)
                    else:
                        dt_up = np.sqrt(2. / acceleration_up * (4.*self.set_up[exp]['height'] / max(self.set_up[exp]['SLF'], self.set_up[exp]['SHF']) -
                                                          p - self.set_up[exp]['height']))
                        dt_down = np.sqrt(2. / acceleration_down * (4.*self.set_up[exp]['height'] / max(self.set_up[exp]['SLF'], self.set_up[exp]['SHF']) -
                                                          p - self.set_up[exp]['height']))
                delays['up'][exp].append(dt_up)
                delays['down'][exp].append(dt_down)
        
    return field_list, delays


def PropCalDiag(param, tauc, OtherInputs, magfield, Increment, field_list, delays_shuttle, set_up, B0LFields, PosAuto):
    Propagators = {}
    
    #Shuttle periods
    shuttle_eigval = []
    shuttle_eigvec = []
    for b0_shuttle in range(len(field_list)):
        relax_mat = np.asarray(RM(field_list[b0_shuttle], param, tauc, OtherInputs))
        Eig, Vec = np.linalg.eig(relax_mat)
        shuttle_eigval.append(Eig)
        shuttle_eigvec.append(Vec)

            
        
    expected_values = {}
    for exp in set_up.keys():
        expected_values[exp] = {}
        #High Field
        relaxmat_HF = np.array(RM(magfield, param, tauc, OtherInputs))
        eigHF, eigvecHF = np.linalg.eig(relaxmat_HF)
        Propagators['HF_1'] = eigvecHF @ np.diag(np.exp(-set_up[exp]['WTHF'] * eigHF)) @ np.linalg.inv(eigvecHF)
        Propagators['HF_2'] = eigvecHF @ np.diag(np.exp(-set_up[exp]['d25'] * eigHF)) @ np.linalg.inv(eigvecHF)
        
        #Low Field
        relaxmat_LF = np.array(RM(B0LFields[exp], param, tauc, OtherInputs))
        eigLF, eigvecLF = np.linalg.eig(relaxmat_LF)
        Propagators['LF'] = {}
        for vc in set_up[exp]['vc']:
            Propagators['LF'][vc] = eigvecLF @ np.diag(np.exp(-set_up[exp]['LF_times'][vc] * eigLF)) @ np.linalg.inv(eigvecLF)
            
        #shuttling
        prop_shuttling_up = []
        for count, delay in enumerate(delays_shuttle['up'][exp]):
            prop_shuttling_up.append(shuttle_eigvec[count] @ np.diag(np.exp(-delay * shuttle_eigval[count])) @ np.linalg.inv(shuttle_eigvec[count]))
        prop_shuttling_down = []
        for count, delay in enumerate(delays_shuttle['down'][exp]):
            prop_shuttling_down.append(shuttle_eigvec[count] @ np.diag(np.exp(-delay * shuttle_eigval[count])) @ np.linalg.inv(shuttle_eigvec[count]))

        Propagators['LF->HF'] = prop_shuttling_down[0]
        for p in prop_shuttling_down[1:]:
            Propagators['LF->HF'] = Propagators['LF->HF'] @ p
        Propagators['HF->LF'] = prop_shuttling_up[0]
        for p in prop_shuttling_up[1:]:
            Propagators['HF->LF'] = p @ Propagators['HF->LF']
        
        #Total Propagation
        for vc in set_up[exp]['vc']:
            total_propagrator = Propagators['HF_2'] @ Propagators['LF->HF'] @ Propagators['LF'][vc] @ Propagators['HF->LF'] @ Propagators['HF_1']
            expected_values[exp][vc] = total_propagrator[PosAuto][PosAuto]
           
    return expected_values


def PropCalDiag_ForPlot(param, tauc, OtherInputs, magfield, Increment, field_list, delays_shuttle, set_up, B0LFields, PosAuto, lf_times, exp):
    Propagators = {}
    
    #Shuttle periods
    shuttle_eigval = []
    shuttle_eigvec = []
    for b0_shuttle in range(len(field_list)):
        relax_mat = np.asarray(RM(field_list[b0_shuttle], param, tauc, OtherInputs))
        Eig, Vec = np.linalg.eig(relax_mat)
        shuttle_eigval.append(Eig)
        shuttle_eigvec.append(Vec)
        
    expected_values = []
    #High Field
    relaxmat_HF = np.array(RM(magfield, param, tauc, OtherInputs))
    eigHF, eigvecHF = np.linalg.eig(relaxmat_HF)
    Propagators['HF_1'] = eigvecHF @ np.diag(np.exp(-set_up[exp]['WTHF'] * eigHF)) @ np.linalg.inv(eigvecHF)
    Propagators['HF_2'] = eigvecHF @ np.diag(np.exp(-set_up[exp]['d25'] * eigHF)) @ np.linalg.inv(eigvecHF)
    
    #Low Field
    relaxmat_LF = np.array(RM(B0LFields[exp], param, tauc, OtherInputs))
    eigLF, eigvecLF = np.linalg.eig(relaxmat_LF)
    Propagators['LF'] = {}
    for vc in lf_times:
        Propagators['LF'][vc] = eigvecLF @ np.diag(np.exp(-vc * eigLF)) @ np.linalg.inv(eigvecLF)
        
    #shuttling
    prop_shuttling_up = []
    for count, delay in enumerate(delays_shuttle['up'][exp]):
        prop_shuttling_up.append(shuttle_eigvec[count] @ np.diag(np.exp(-delay * shuttle_eigval[count])) @ np.linalg.inv(shuttle_eigvec[count]))
    prop_shuttling_down = []
    for count, delay in enumerate(delays_shuttle['down'][exp]):
        prop_shuttling_down.append(shuttle_eigvec[count] @ np.diag(np.exp(-delay * shuttle_eigval[count])) @ np.linalg.inv(shuttle_eigvec[count]))

    Propagators['LF->HF'] = prop_shuttling_down[0]
    for p in prop_shuttling_down[1:]:
        Propagators['LF->HF'] = Propagators['LF->HF'] @ p
    Propagators['HF->LF'] = prop_shuttling_up[0]
    for p in prop_shuttling_up[1:]:
        Propagators['HF->LF'] = p @ Propagators['HF->LF']

    #Total Propagation
    for vc in lf_times:
        total_propagrator = Propagators['HF_2'] @ Propagators['LF->HF'] @ Propagators['LF'][vc] @ Propagators['HF->LF'] @ Propagators['HF_1']
        expected_values.append(total_propagrator[PosAuto][PosAuto])
           
    return np.asarray(expected_values)


def PropCalExp(param, tauc, OtherInputs, magfield, Increment, field_list, delays_shuttle, set_up, B0LFields, PosAuto):
    Propagators = {}
    
    #Shuttle periods
    relax_mat = []
    for b0_shuttle in range(len(field_list)):
        relax_mat.append(np.asarray(RM(field_list[b0_shuttle], param, tauc, OtherInputs)))


    expected_values = {}
    for exp in set_up.keys():
        expected_values[exp] = {}
        #High Field
        relaxmat_HF = np.array(RM(magfield, param, tauc, OtherInputs))
        Propagators['HF_1'] = linalg.expm(-set_up[exp]['WTHF'] * relaxmat_HF)
        Propagators['HF_2'] = linalg.expm(-set_up[exp]['d25'] * relaxmat_HF)
        
        #Low Field
        relaxmat_LF = np.array(RM(B0LFields[exp], param, tauc, OtherInputs))
        Propagators['LF'] = {}
        for vc in set_up[exp]['vc']:
            Propagators['LF'][vc] = linalg.expm(-set_up[exp]['LF_times'][vc] * relaxmat_LF)
            
        #shuttling
        prop_shuttling_up = []
        for count, delay in enumerate(delays_shuttle['up'][exp]):
            prop_shuttling_up.append(linalg.expm(-delay * relax_mat[count]))
        prop_shuttling_down = []
        for count, delay in enumerate(delays_shuttle['down'][exp]):
            prop_shuttling_down.append(linalg.expm(-delay * relax_mat[count]))

        Propagators['LF->HF'] = prop_shuttling_down[0]
        for p in prop_shuttling_down[1:]:
            Propagators['LF->HF'] = Propagators['LF->HF'] @ p
        Propagators['HF->LF'] = prop_shuttling_up[0]
        for p in prop_shuttling_up[1:]:
            Propagators['HF->LF'] = p @ Propagators['HF->LF']
        
        #Total Propagation
        for vc in set_up[exp]['vc']:
            total_propagrator = Propagators['HF_2'] @ Propagators['LF->HF'] @ Propagators['LF'][vc] @ Propagators['HF->LF'] @ Propagators['HF_1']
            expected_values[exp][vc] = total_propagrator[PosAuto][PosAuto]
           
    return expected_values


def optShuttling(self, ExpNum, ParamOpt, PosAuto):
    Increment = 1e-3    #starting increment in meter
    
    field_list, delays = FieldList(self, Increment)
    
    sub_setup = {}
    sub_setup[ExpNum] = self.set_up[ExpNum]
    sub_B0LF = {}
    sub_B0LF[ExpNum] = self.B0LFields[ExpNum]
    
    ExpVal_Init = []
    aa = list(self.OtherInputs.keys())[0]
    for param in ParamOpt:
        ExpVal_Init_Full = PropCalDiag(param, self.TauC, self.OtherInputs[aa], self.MagField, Increment, field_list, delays, sub_setup, sub_B0LF, PosAuto)
        ExpVal_Init.append(ExpVal_Init_Full[ExpNum][max(sub_setup[ExpNum]['vc'])])
    print(f" Initial value : {Increment} m")
        
    count = 0
    while True:
        count += 1
        Increment += 1e-3
        field_list, delays = FieldList(self, Increment)
        
        for counter, param in enumerate(ParamOpt):
            ExpVal_Full = PropCalDiag(param, self.TauC, self.OtherInputs[aa], self.MagField, Increment, field_list, delays, sub_setup, sub_B0LF, PosAuto)
            ExpVal = ExpVal_Full[ExpNum][max(sub_setup[ExpNum]['vc'])]
        
            if abs(ExpVal-ExpVal_Init[counter])/ExpVal_Init[counter] > 0.01:
                return Increment - 1e-3
            
        print(f" Updated value {count} : {Increment} m")
        
        
