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

import FitFunctions as FitF
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
        
    defaultVal.write("Global tumbling correlation time\t")
    defaultVal.write(str(TauC))

    defaultVal.write("\n")
    defaultVal.write("Acceleration type\t")
    defaultVal.write(str(self.AccelerationTYPE.get()))
    defaultVal.write("\n")
    defaultVal.write("Number of chains\t")
    defaultVal.write(str(Nchains))
    defaultVal.write("\n")
    defaultVal.write("Number of MCMC steps\t")
    defaultVal.write(str(Nstep))
    defaultVal.write("\n")
    defaultVal.write("PDB ID\t")
    defaultVal.write(self.PDBid.get())
    defaultVal.write("\n")
    defaultVal.write("Setup path\t")
    defaultVal.write(str(ExperimentalSetUpPath))
    defaultVal.write("\n")
    defaultVal.write("Field calibration path\t")
    defaultVal.write(str(FieldCalibrationPath))
    defaultVal.write("\n")
    defaultVal.write("Relaxometry data path\t")
    defaultVal.write(str(IntRelaxoPath))
    defaultVal.write("\n")
    defaultVal.write("Ohter Input path\t")
    defaultVal.write(str(InputPath))
    defaultVal.write("\n")
        
    for i in range(len(self.OP)):
        mins2 = self.MinS2[i].get()
        maxs2 = self.MaxS2[i].get()
        
        defaultVal.write(self.OP[i])
        defaultVal.write("\t")
        defaultVal.write(str(mins2))
        defaultVal.write("\t")
        defaultVal.write(str(maxs2))
        defaultVal.write("\t")
    
    defaultVal.write("\n")
    for i in range(len(self.CT)):
        mintau = self.MinTau[i].get()
        maxtau = self.MaxTau[i].get()
        
        defaultVal.write(self.CT[i])
        defaultVal.write("\t")
        defaultVal.write(str(mintau))
        defaultVal.write("\t")
        defaultVal.write(str(maxtau))
        defaultVal.write("\t")
        
    defaultVal.write("\n")
    for i in range(len(self.Others)):
        minOthers = self.MinOthers[i].get()
        maxOthers = self.MaxOthers[i].get()
        
        defaultVal.write(self.Others[i])
        defaultVal.write("\t")
        defaultVal.write(str(minOthers))
        defaultVal.write("\t")
        defaultVal.write(str(maxOthers))
        defaultVal.write("\t")
            
    for DataSet in range(len(self.RATES)):
        defaultVal.write("\n")
        defaultVal.write("HF data set " + str(DataSet+1) + "\t")

        rate = self.RATES[DataSet]
        Type = self.RelaxationDataType[DataSet].get()
        Field = self.Fields[DataSet].get()

        defaultVal.write(str(rate))
        defaultVal.write("\t")
        defaultVal.write(str(Type))
        defaultVal.write("\t")
        defaultVal.write(str(Field))
        
    defaultVal.close()
        
        
def writeParam(self, loadFile, TauC, AccelerationType, Nmcmc, Nwalkers, checkPDB, PDB, ExperimentalSetUpPath, FieldCalibrationPath, IntRelaxoPath, InputPath):
    defaultVal = open(loadFile, 'w')
        
    defaultVal.write("UNCORRECT\n")
    defaultVal.write("Global tumbling correlation time\t")
    defaultVal.write(str(TauC))
    
    defaultVal.write("\n")
    defaultVal.write("Acceleration type\t")
    defaultVal.write(AccelerationType)
    defaultVal.write("\n")
    defaultVal.write("Number of Chains\t")
    defaultVal.write(str(Nwalkers))
    defaultVal.write("\n")
    defaultVal.write("Number of MCMC steps\t")
    defaultVal.write(str(Nmcmc))
    defaultVal.write("\n")
    defaultVal.write("PDB ID\t")
    if checkPDB:
        defaultVal.write(PDB)
    else:
        defaultVal.write('Not available')
    defaultVal.write("\n")
    defaultVal.write("Setup Path\t")
    defaultVal.write(str(ExperimentalSetUpPath))
    defaultVal.write("\n")
    defaultVal.write("Field calibration path\t")
    defaultVal.write(str(FieldCalibrationPath))
    defaultVal.write("\n")
    defaultVal.write("Relaxometry data path\t")
    defaultVal.write(str(IntRelaxoPath))
    defaultVal.write("\n")
    defaultVal.write("Other inputs path\t")
    defaultVal.write(str(InputPath))
    defaultVal.write("\n")
        
    bnds = []
    for i in range(len(self.OP)):
        mins2 = float(self.MinS2[i].get())
        maxs2 = float(self.MaxS2[i].get())
        bnds.append((mins2, maxs2))

        defaultVal.write(self.OP[i])
        defaultVal.write("\t")
        defaultVal.write(str(mins2))
        defaultVal.write("\t")
        defaultVal.write(str(maxs2))
        defaultVal.write("\t")
        
    defaultVal.write("\n")
    for i in range(len(self.CT)):
        mintau = float(self.MinTau[i].get())
        maxtau = float(self.MaxTau[i].get())
        bnds.append((mintau, maxtau))
        
        defaultVal.write(self.CT[i])
        defaultVal.write("\t")
        defaultVal.write(str(mintau))
        defaultVal.write("\t")
        defaultVal.write(str(maxtau))
        defaultVal.write("\t")
        
    defaultVal.write("\n")
    for i in range(len(self.Others)):
        minOthers = float(self.MinOthers[i].get())
        maxOthers = float(self.MaxOthers[i].get())
        bnds.append((minOthers, maxOthers))
        
        defaultVal.write(self.Others[i])
        defaultVal.write("\t")
        defaultVal.write(str(minOthers))
        defaultVal.write("\t")
        defaultVal.write(str(maxOthers))
        defaultVal.write("\t")
            
    for DataSet in range(len(self.RelaxationDataSet)):
        defaultVal.write("\n")
        defaultVal.write("HF data set " + str(DataSet+1) + "\t")
        defaultVal.write(str(self.RATES[DataSet]))
        defaultVal.write("\t")
        defaultVal.write(str(self.RelaxationDataType[DataSet].get()))
        defaultVal.write("\t")
        defaultVal.write(str(self.Fields[DataSet].get()))
        
        
    defaultVal = open(loadFile, 'w')
    defaultVal.write("CONFIRMED\n")
    defaultVal.close()
    
    return bnds



        
def WriteMCMCTraj(self, Traj, AA):
    OutName = self.directoryName + "/FittingResults/Trajectories/Trajectory_Residue" + str(AA) + ".txt"
        
    trajFile = open(OutName, 'w')
    
    trajFile.write("Step")
    for i in range(len(self.TotParam)):
        trajFile.write("\t")
        trajFile.write(self.TotParam[i])
    trajFile.write("\tf")
        
    for i in range(len(Traj)):
        trajFile.write("\n")
        trajFile.write(str(i))
        for P in range(len(Traj[i])):
            trajFile.write("\t")
            trajFile.write(str(Traj[i][P]))
    trajFile.close()
    
    
    
    
def WriteMCMCParam(self):
    Scaling = [FitF.ScalingFactor(self.FinalSimulatedIntensities[AA], self.Intensities[AA]) for AA in range(len(self.AAList))]
    AverageScal = [sum(Scaling[AA])/float(len(Scaling[AA])) for AA in range(len(self.AAList))]
    ScalingScore = [[Scaling[AA][i]**2 for i in range(len(Scaling[AA]))] for AA in range(len(self.AAList))]
    VarianceScal = [sum(ScalingScore[AA])/float(len(ScalingScore[AA])) - AverageScal[AA]**2 for AA in range(len(self.AAList))]
        
        
    OutName = self.directoryName + "/FittingResults/MCMCparameters.txt"
    OutFile = open(OutName, 'w')
    OutFile.write("Residue\tMAF\tAverage scaling factor\tVariance scaling factor")
    for L in range(len(self.TotParam)):
        OutFile.write("\t")
        OutFile.write(str(self.TotParam[L]))
        OutFile.write("\t")
        OutFile.write("+ error")
        OutFile.write("\t")
        OutFile.write("- error")
    OutFile.write("lnf\t+ error\t- error")
    for AA in range(len(self.AAList)):
        OutFile.write("\n")
        OutFile.write(str(self.AAList[AA]))
        OutFile.write("\t")
        OutFile.write(str(self.Acceptance[AA]))
        OutFile.write("\t")
        OutFile.write(str(AverageScal[AA]))
        OutFile.write("\t")
        OutFile.write(str(VarianceScal[AA]))
        for L in range(len(self.TotParam)):
            OutFile.write("\t")
            OutFile.write(str(self.MCMCparam[AA][0][L]))
            OutFile.write("\t")
            OutFile.write(str(self.MCMCparam[AA][1][L]))
            OutFile.write("\t")
            OutFile.write(str(self.MCMCparam[AA][2][L]))
        OutFile.write("\t")
        OutFile.write(str(self.MCMCparam[AA][0][-1]))
        OutFile.write("\t")
        OutFile.write(str(self.MCMCparam[AA][1][-1]))
        OutFile.write("\t")
        OutFile.write(str(self.MCMCparam[AA][2][-1]))
        
        
    OutFile.close()
    
    
        
        
def WritePDB(self, AllChi2):
    
    if len(self.AAList) == 1:
        print('Only one residue provided: PDB coloring files will not be created')
        
    else:
    
        dirPDB = self.directoryName + "/FittingResults/PDBFiles"
        os.makedirs(dirPDB)
        maxChi2 = max(AllChi2)
        minChi2 = min(AllChi2)
        maxParam = [[] for i in range(len(self.TotParam))]
        minParam = [[] for i in range(len(self.TotParam))]
        
        AllParam = [[[] for AA in self.AAList] for P in self.TotParam]
        for AA in range(len(self.AAList)):
            for P in range(len(self.TotParam)):
                AllParam[P][AA] = self.MCMCparam[AA][0][P]
        
        for param in range(len(AllParam)):
            maxParam[param] = max(AllParam[param])
            minParam[param] = min(AllParam[param])
                    
        NormalizedParam = [[] for AA in self.AAList]
        for AA in range(len(self.AAList)):
            NormalizedParam[AA].append((AllChi2[AA]-minChi2)/(maxChi2-minChi2))
            for param in range(len(AllParam)):
                NormalizedParam[AA].append((AllParam[param][AA]-minParam[param])/(maxParam[param]-minParam[param]))
                        
        for param in range(len(AllParam)+1):
            if param == 0:
                PDBfilename = dirPDB + "/" + self.PDB + "_Chi2.pml"
                PDBfile = open(PDBfilename, 'w')
            else:
                PDBfilename = dirPDB + "/" + self.PDB + "_" + str(self.TotParam[param-1]) + ".pml"
                PDBfile = open(PDBfilename, 'w')
                    
            PDBfile.write("load ")
            PDBfile.write(str(self.PDB))
            PDBfile.write(".pdb\n\nhide\n\nshow cartoon, all\n\n")
                    
            PDBfile.write("#Color the whole structure in grey as default color\n")
            PDBfile.write("color grey, ")
            PDBfile.write(str(self.PDB))
            PDBfile.write("\n\n")
                    
            PDBfile.write("#Define a color in the blue to red range for each residue we have data for")
            for AA in range(len(self.AAList)):
                PDBfile.write("\n")
                PDBfile.write("set_color resi")
                PDBfile.write(str(self.AAList[AA]))
                PDBfile.write(" = [")
                PDBfile.write(str(NormalizedParam[AA][param]))
                PDBfile.write(", 0.0, 1.0-")
                PDBfile.write(str(NormalizedParam[AA][param]))
                PDBfile.write("]")
    
            PDBfile.write("\n#Coloring each residue according to its color")
            for AA in range(len(self.AAList)):
                PDBfile.write("\ncolor resi")
                PDBfile.write(str(self.AAList[AA]))
                PDBfile.write(", resi ")
                PDBfile.write(str(self.AAList[AA]))
                        
        PDBfile.close()
            

def WriteLFR1(self):
    filenameRelaxSpeed = self.directoryName + "/FittingResults/LowFieldRelaxationR1.txt"
    frelax = open(filenameRelaxSpeed, 'w')
    
    frelax.write("\t\tField\tIntensity scaling\tR1 back-calculated\tFitted R1")
            
    for AA in range(len(self.AAList)):
        frelax.write("\n")
        Res = self.AAList[AA]
        frelax.write("Residue ")
        frelax.write(str(Res))
                        
        for Field in range(len(self.B0LFields[AA])):
            frelax.write("\n")
            frelax.write(str(self.B0LFields[AA][Field]))
            frelax.write("\t")
            frelax.write(str(self.ScalingIntensities[AA][Field]))
            frelax.write("\t")
            frelax.write(str(self.R1LFDataForCurve_BackCalc[AA][Field]))
            frelax.write("\t")
            frelax.write(str(self.R1LFDataForCurve_Fitted[AA][Field][1]))
            
    frelax.close()