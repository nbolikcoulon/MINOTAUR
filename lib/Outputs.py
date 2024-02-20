#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##########################################################################
#                                                                        #
#                       Creates output text files:                       #
#       - parameters used for the analysis                               #
#       - dynamics parameters                                            #
#       - script to visualize the relaxation matrix and relaxation rates #
#       - PDB file                                                       #
#       - Low field R1                                                   #
#                                                                        #
##########################################################################
from MonteCarlo import ScalingFactor
import numpy as np
import os

            
def WritePreSavedParam(self, IcarusDir):
    SavedFile = IcarusDir + "/PreSaved_1_Parameters.txt"
    
    t = 0
    n = 1
    while t == 0:
        if not os.path.exists(SavedFile):
            t = 1
        else:
            SavedFile = IcarusDir + "/PreSaved_" + str(n) + "_Parameters.txt"
            n += 1
    defaultVal = open(SavedFile, 'w')
        
    ExperimentalSetUpPath = "Not provided"
    FieldCalibrationPath = "Not provided"
    IntRelaxoPath = "Not provided"
    InputPath = "Not provided"
    if hasattr(self, "ExperimentalSetUp"):
        ExperimentalSetUpPath = self.ExperimentalSetUp
    if hasattr(self, "FieldCalibration"):
        FieldCalibrationPath = self.FieldCalibration
    if hasattr(self, "Intrelax"):
        IntRelaxoPath = self.Intrelax
    if hasattr(self, "InputFile"):
        InputPath = self.InputFile
        
    Nchains = int(self.Nwalker.get())
    Nstep = int(self.Nmcmc.get())

    TauC = self.TauC.get()
                
    defaultVal.write("PRE SAVED\n")
        
    defaultVal.write(f"Global tumbling correlation time\t{TauC}\n")
    defaultVal.write(f"Shuttling type\t{self.Shuttling_TYPE.get()}\n")
    defaultVal.write(f"Number of chains\t{Nchains}\n")
    defaultVal.write(f"Number of MCMC steps\t{Nstep}\n")
    defaultVal.write(f"PDB ID\t{self.PDBid.get()}\n")
    defaultVal.write(f"Setup path\t{ExperimentalSetUpPath}\n")
    defaultVal.write(f"Field calibration path\t{FieldCalibrationPath}\n")
    defaultVal.write(f"Relaxometry data path\t{IntRelaxoPath}\n")
    defaultVal.write(f"Ohter Input path\t{InputPath}\n")
        
    for i in range(len(self.OP)):
        mins2 = self.MinS2[i].get()
        maxs2 = self.MaxS2[i].get()
        defaultVal.write(f'{self.OP[i]}\t{mins2}\t{maxs2}\t')
    
    defaultVal.write("\n")
    for i in range(len(self.CT)):
        mintau = self.MinTau[i].get()
        maxtau = self.MaxTau[i].get()
        defaultVal.write(f'{self.CT[i]}\t{mintau}\t{maxtau}\t')
        
    defaultVal.write("\n")
    for i in range(len(self.Others)):
        minOthers = self.MinOthers[i].get()
        maxOthers = self.MaxOthers[i].get()
        defaultVal.write(f'{self.Others[i]}\t{minOthers}\t{maxOthers}\t')
            
    for DataSet in range(len(self.RATES)):
        defaultVal.write("\n")
        defaultVal.write(f"HF data set {DataSet+1}\t{self.RATES[DataSet]}\t{self.RelaxationDataType[DataSet].get()}\t{self.Fields[DataSet].get()}")

    defaultVal.close()
        
        
def writeParam(self, loadFile, TauC, Shuttling_Type, Nmcmc, Nwalkers, checkPDB, PDB, ExperimentalSetUpPath, FieldCalibrationPath, IntRelaxoPath, InputPath):
    defaultVal = open(loadFile, 'w')
        
    defaultVal.write("UNCORRECT\n")
    defaultVal.write(f"Global tumbling correlation time\t{TauC}\n")
    defaultVal.write(f"Shuttling type\t{Shuttling_Type}\n")
    defaultVal.write(f"Number of Chains\t{Nwalkers}\n")
    defaultVal.write(f"Number of MCMC steps\t{Nmcmc}\n")
    defaultVal.write("PDB ID\t")
    if checkPDB:
        defaultVal.write(PDB)
    else:
        defaultVal.write('Not available')
    defaultVal.write("\n")
    defaultVal.write(f"Setup Path\t{ExperimentalSetUpPath}\n")
    defaultVal.write(f"Field calibration path\t{FieldCalibrationPath}\n")
    defaultVal.write(f"Relaxometry data path\t{IntRelaxoPath}\n")
    defaultVal.write(f"Other inputs path\t{InputPath}\n")
        
    bnds = []
    for i in range(len(self.OP)):
        mins2 = float(self.MinS2[i].get())
        maxs2 = float(self.MaxS2[i].get())
        bnds.append((mins2, maxs2))
        defaultVal.write(f'{self.OP[i]}\t{mins2}\t{maxs2}\t')
        
    defaultVal.write("\n")
    for i in range(len(self.CT)):
        mintau = float(self.MinTau[i].get())
        maxtau = float(self.MaxTau[i].get())
        bnds.append((mintau, maxtau))
        defaultVal.write(f'{self.CT[i]}\t{mintau}\t{maxtau}\t')
        
    defaultVal.write("\n")
    for i in range(len(self.Others)):
        minOthers = float(self.MinOthers[i].get())
        maxOthers = float(self.MaxOthers[i].get())
        bnds.append((minOthers, maxOthers))
        defaultVal.write(f'{self.Others[i]}\t{minOthers}\t{maxOthers}\t')
            
    for DataSet in range(len(self.RelaxationDataSet)):
        defaultVal.write("\n")
        defaultVal.write(f"HF data set {DataSet+1}\t{self.RATES[DataSet]}\t{self.RelaxationDataType[DataSet].get()}\t{self.Fields[DataSet].get()}")
        
        
    defaultVal = open(loadFile, 'w')
    defaultVal.write("CONFIRMED\n")
    defaultVal.close()
    
    return bnds

        
def WriteMCMCTraj(self, Traj, AA):
    trajFile = open(f'{self.directory_name}/FittingResults/Trajectories/Trajectory_Residue{AA}.txt', 'w')
    
    trajFile.write('Step')
    for param in self.TotParam:
        trajFile.write(f'\t{param}')
    trajFile.write('\tf')
        
    for count, step in enumerate(Traj):
        trajFile.write('\n{count}')
        for param in step:
            trajFile.write(f'\t{param}')
    trajFile.close()
    
    
def WriteMCMCParam(self, AAList):
    OutName = self.directory_name + "/FittingResults/MCMCparameters.txt"
    OutFile = open(OutName, 'w')
    OutFile.write("Residue\tMAF\tAverage scaling factor\tVariance scaling factor")
    for L in self.TotParam:
        OutFile.write(f"\t{L}\t+ error\t- error")
    OutFile.write("lnf\t+ error\t- error")
    for AA in AAList:
        Scaling = ScalingFactor(self.FinalSimulatedIntensities[AA], self.Intensities, self.Err_Int, AA)
        AverageScal = np.average(list(Scaling.values()))
        VarianceScal = np.std(list(Scaling.values()))
    
        OutFile.write(f"\n{AA}\t{self.Acceptance[AA]}\t{AverageScal}\t{VarianceScal}")

        for L in range(len(self.TotParam)+1):
            OutFile.write(f"\t{self.MCMCparam[AA][0][L]}\t{self.MCMCparam[AA][1][L]}\t{self.MCMCparam[AA][2][L]}")
        
    OutFile.close()
        
        
def WritePDB(self, AllChi2):
    
    if len(AllChi2) == 1:
        print('Only one residue provided: PDB coloring files will not be created')
        return
        
    dirPDB = f'{self.directory_name}/FittingResults/PDBFiles'
    os.makedirs(dirPDB)
    
    maxParam = {}
    minParam = {}
    
    labels = np.copy(self.TotParam)

    AllParam = {}
    for count, param in enumerate(labels):
        AllParam[param] = []
        for AA in self.MCMCparam.keys():
            AllParam[param].append(self.MCMCparam[AA][0][count])
        minParam[param], maxParam[param] = min(AllParam[param]), max(AllParam[param])
        
    labels = np.append(labels, 'Chi2')
    AllParam['Chi2'] = []
    for AA in self.MCMCparam.keys():
        AllParam['Chi2'].append(AllChi2[AA])
    minParam['Chi2'], maxParam['Chi2'] = min(list(AllChi2.values())), max(list(AllChi2.values()))

    NormalizedParam = {}
    for param in labels:
        NormalizedParam[param] = {}
        for count, AA in enumerate(list(self.MCMCparam.keys())):
            NormalizedParam[param][AA] = (AllParam[param][count] - minParam[param]) / (maxParam[param] - minParam[param])
                    
    for param in labels:
        PDBfile = open(f'{dirPDB}/{self.PDB}_{param}.pml', 'w')
                
        PDBfile.write(f"load {self.PDB}.pdb\n\n")
        PDBfile.write("hide\n\n")
        PDBfile.write("show cartoon, all\n\n")
                
        PDBfile.write("#Color the whole structure in grey as default color\n")
        PDBfile.write(f"color grey, {self.PDB}\n\n")
                
        PDBfile.write("#Define a color in the blue to red range for each residue we have data for\n")
        for AA in NormalizedParam[param].keys():
            PDBfile.write(f"set_color resi{AA} = [{NormalizedParam[param][AA]}, 0.0, {1.0 - NormalizedParam[param][AA]}]\n")

        PDBfile.write("\n#Coloring each residue according to its color\n")
        for AA in NormalizedParam[param].keys():
            PDBfile.write("color resi{AA}, resi {AA}\n")
                    
        PDBfile.close()
        
        
def WriteLFR1(self):
    frelax = open(f'{self.directory_name}/FittingResults/LowFieldRelaxationR1.txt', 'w')
    
    frelax.write("\t\tField\tIntensity scaling\tR1 back-calculated\tFitted R1\n")
            
    for AA in self.ScalingIntensities.keys():
        frelax.write(f"Residue {AA}\n")
                        
        for exp in self.set_up.keys():
            frelax.write(f"{self.B0LFields[exp]}\t{self.ScalingIntensities[AA][exp]}\t{self.R1LF_BackCalc[AA][exp]}\t{self.R1LF_Fitted[AA][exp]}\n")
            
    frelax.close()
