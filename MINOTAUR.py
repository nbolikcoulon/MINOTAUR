#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tkinter as tk
from tkinter.filedialog import askopenfilename, askdirectory
from PIL import ImageTk, Image
from shutil import copyfile
from distutils.dir_util import copy_tree
from scipy.optimize import curve_fit
from random import uniform
from scipy import linalg

import numpy as np
import time
import sys
import os



CSFolder = str(sys.argv[1])
RatesFolder = CSFolder + "/Rates"

sys.path.append(os.path.abspath(CSFolder))
sys.path.append(os.path.abspath(RatesFolder))
import Parameters as ParamFile
import _RelaxMat

sys.path.append(os.path.abspath('lib'))
import FitFunctions as FitF
import ShuttlingSimulation as ShSim
import Outputs as Out
import FigOutputs as FigOut
import MonteCarlo as MCMC


#################################################### Graphical User Interface ####################################################
class GUI(tk.Frame):
    def __init__(self, parent):
        
        tk.Frame.__init__(self, parent)
        
        self.parent = parent
        
        self.RelaxFunc = ParamFile.ImportFunc()
        self.RATES = []
        
        self.parameters()

############ Main window        
    def parameters(self):
        
        
        def Quit(self):
            sys.exit(0)
            
        def SaveParam(self):
            CurrentDir = os.path.dirname(os.path.abspath(__file__))
            
            Out.WritePreSavedParam(self, CurrentDir)
        
        #Browsing functions
        def BrowseFieldCal(self):
            self.FieldCalibration = askopenfilename()
            self.FieldCal.config(text=os.path.basename(self.FieldCalibration))
            
        def BrowseExpSetUp(self):
            self.ExperimentalSetUp = askopenfilename()
            self.ExpSetUp.config(text=os.path.basename(self.ExperimentalSetUp))
            
        def BrowseRelaxomInt(self):
            self.Intrelax = askdirectory()
            self.Irelaxom.config(text=os.path.basename(self.Intrelax))
            
        def BrowseInputs(self):
            self.InputFile = askopenfilename()
            self.Input.config(text=os.path.basename(self.InputFile))
            
            
        def BrowseRates(self, n):
            n = n[0]
            if n+1 > len(self.RATES):
                self.RATES.append(askopenfilename())
            else:
                self.RATES[n] = askopenfilename()
                
            self.RelaxationDataSet[n].config(text=os.path.basename(self.RATES[n]))
            
            if len(self.RelaxationDataSet) % 2 == 0:
                self.AddButton.grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=8)
            else:
                self.AddButton.grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=4, columnspan=2)
                
        def Add(self):
            N = ([self.Nlabel])
            self.Nlabel += 1
            self.AddButton.grid_forget()
            
            self.RelaxationDataSet.append(tk.Button(self.param, text="High field rates", command = lambda: BrowseRates(self, N)))
    
            self.RelaxationDataType.append(tk.StringVar())
            self.RelaxationDataType[-1].set("Data Type")
            self.RelaxationType.append(tk.OptionMenu(self.param, self.RelaxationDataType[-1], *self.MenuRelaxType.keys()))
    
            Field = tk.StringVar()
            Field.set("High field (T)")
            self.Fields.append(tk.Entry(self.param, textvariable = Field, width=10))
            
            if len(self.RelaxationDataSet) % 2 == 0:
                self.RelaxationDataSet[-1].grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=4, columnspan=2)
                self.RelaxationType[-1].grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=6)
                self.Fields[-1].grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=7)
                
            else:
                self.RelaxationDataSet[-1].grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=0, columnspan=2)
                self.RelaxationType[-1].grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=2)
                self.Fields[-1].grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=3)
                
                
        #function to load previous parameters
        def LoadParam(self):
            root.withdraw()
            self.parameters = askopenfilename()
            
            paramFile = open(self.parameters, 'r')
            nparam = sum(1 for line in paramFile)-1
            paramFile.close()
            paramFile = open(self.parameters, 'r')
            
            paramFile.readline()
    
            previousSetUp = []
            for i in range(nparam):
                Line = paramFile.readline()
                Line = [x for x in Line.split('\n')]
                Line = Line[0]
                Line = [x for x in Line.split('\t')]
                previousSetUp.append(Line)
                
            self.TAUC.set(previousSetUp[0][1])
                    
            self.AccelerationTYPE.set(previousSetUp[1][1])
            self.NWALKER.set(previousSetUp[2][1])
            self.NMCMC.set(previousSetUp[3][1])
            if previousSetUp[4][1] != 'Not available':
                self.PDBID.set(previousSetUp[4][1])
                        
            if os.path.exists(previousSetUp[5][1]):
                self.ExperimentalSetUp = previousSetUp[5][1]
                self.ExpSetUp.config(text=os.path.basename(self.ExperimentalSetUp))
            else:
                self.ErrorText = tk.Label(self.param, text = "Status: Missing file(s)", fg="orange")
                self.ErrorText.grid(row=0, column=3)
                
            if os.path.exists(previousSetUp[6][1]):
                self.FieldCalibration = previousSetUp[6][1]
                self.FieldCal.config(text=os.path.basename(self.FieldCalibration))
            else:
                self.ErrorText = tk.Label(self.param, text = "Status: Missing file(s)", fg="orange")
                self.ErrorText.grid(row=0, column=3)
                
            if os.path.exists(previousSetUp[7][1]):
                self.Intrelax = previousSetUp[7][1]
                self.Irelaxom.config(text=os.path.basename(self.Intrelax))
            else:
                self.ErrorText = tk.Label(self.param, text = "Status: Missing file(s)", fg="orange")
                self.ErrorText.grid(row=0, column=3)
                    
            if os.path.exists(previousSetUp[8][1]):
                self.InputFile = previousSetUp[8][1]
                self.Input.config(text=os.path.basename(self.InputFile))
            else:
                self.ErrorText = tk.Label(self.param, text = "Status: Missing file(s)", fg="orange")
                self.ErrorText.grid(row=0, column=3)
                
                    
                    
            if len(previousSetUp[9]) == 3*len(self.mins2)+1:
                for i in range(len(self.mins2)):
                    self.mins2[i].set(previousSetUp[9][3*i+1])
                    self.maxs2[i].set(previousSetUp[9][3*i+2])
            if len(previousSetUp[10]) == 3*len(self.mintau)+1:
                for i in range(len(self.mintau)):
                    self.mintau[i].set(previousSetUp[10][3*i+1])
                    self.maxtau[i].set(previousSetUp[10][3*i+2])
            if len(previousSetUp[11]) == 3*len(self.minOthers)+1:
                for i in range(len(self.minOthers)):
                    self.minOthers[i].set(previousSetUp[11][3*i+1])
                    self.maxOthers[i].set(previousSetUp[11][3*i+2])
                    
    
            self.AddButton.grid_forget()
            for i in range(len(self.RelaxationDataSet)):
                self.RelaxationDataSet[i].grid_forget()
                self.RelaxationType[i].grid_forget()
                self.Fields[i].grid_forget()  
                
            self.RATES = []
            self.RelaxationDataSet = []
            self.RelaxationDataType = []
            self.RelaxationType = []
            self.Fields = []
            self.Nlabel = 0
            for i in range(12, nparam):
                if os.path.exists(previousSetUp[i][1]):
                    Add(self)
                    self.RelaxationDataSet[-1].config(text=os.path.basename(previousSetUp[i][1]))
                    self.RATES.append(previousSetUp[i][1])
                    self.RelaxationDataType[-1].set(previousSetUp[i][2])
                    self.Fields[-1].delete(0, "end")
                    self.Fields[-1].insert(0, previousSetUp[i][3])
                        
                    if len(self.RelaxationDataSet) % 2 == 0:
                        self.RelaxationDataSet[-1].grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=4, columnspan=2)
                        self.RelaxationType[-1].grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=6)
                        self.Fields[-1].grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=7)
                    else:
                        self.RelaxationDataSet[-1].grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=0, columnspan=2)
                        self.RelaxationType[-1].grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=2)
                        self.Fields[-1].grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=3)
                else:
                    self.ErrorText = tk.Label(self.param, text = "Status: Missing file(s)", fg="orange")
                    self.ErrorText.grid(row=0, column=3)
                
            if len(self.RelaxationDataSet) % 2 == 0:
                self.AddButton.grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=8, columnspan=2)
            else:
                self.AddButton.grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=4, columnspan=2)
            
############# GUI coding
        self.param = tk.Toplevel()
        self.param.title("MINOTAUR")
        self.param.geometry("1300x800")
        self.grid(row=0, column=0)
        
        
        self.img = ImageTk.PhotoImage(Image.open("lib/Logo.jpg"))
        panel = tk.Label(self.param, image=self.img)
        panel.grid(row=0, column=2, columnspan=3, rowspan=2)
        
#Button part
        self.quitButton = tk.Button(self.param, text="Quit", command = lambda: Quit(self))
        self.startButton = tk.Button(self.param, text="Start", bg="green", command=self.GetData)
        
        self.quitButton.grid(row=0, column=6)
        self.startButton.grid(row=0, column=7)

        self.Load = tk.Button(self.param, text="Load previous parameters", command = lambda: LoadParam(self))
        self.Load.grid(row=0, column=0, columnspan=3)
        
        self.Save = tk.Button(self.param, text="Save current parameters", command = lambda: SaveParam(self))
        self.Save.grid(row=1, column=6, columnspan=2)
        
#Fitting parameters part
        ParamText = tk.Label(self.param, text="Set fitting parameters:", fg="blue")
        
        self.TAUC = tk.IntVar()
        self.TAUC.set("Correlation time TauC (sec)")
        self.TauC = tk.Entry(self.param, textvariable = self.TAUC)

        self.AccelerationTYPE = tk.StringVar()
        self.AccelerationTYPE.set("Type of shuttling")
        AccelerationType = tk.OptionMenu(self.param, self.AccelerationTYPE, "Constant Speed", "Constant Acceleration")
        
        self.NMCMC= tk.IntVar()
        self.NMCMC.set("MCMC - Number of steps")
        self.Nmcmc = tk.Entry(self.param, textvariable = self.NMCMC, width=25)
        
        self.NWALKER = tk.IntVar()
        self.NWALKER.set("MCMC - Number of chains")
        self.Nwalker = tk.Entry(self.param, textvariable = self.NWALKER, width=25)
        
        
        ParamText.grid(row=1, column=0)
        self.TauC.grid(row=2, column=0, columnspan=3)
        
        AccelerationType.grid(row=2, column=3, columnspan=2)

        self.Nwalker.grid(row=3, column=5, columnspan=3)
        self.Nmcmc.grid(row=2, column=5, columnspan=3)
        
#Limits Part
        LimitText = tk.Label(self.param, text="Set Limits for MCMC initial guess:", fg="blue")

        MinS2Text = tk.Label(self.param, text="Min")
        MaxS2Text = tk.Label(self.param, text="Max")
        MinTauText = tk.Label(self.param, text="Min")
        MaxTauText = tk.Label(self.param, text="Max")
        MinOtherText = tk.Label(self.param, text="Min")
        MaxOtherText = tk.Label(self.param, text="Max")
        
        LimitText.grid(row=5, column=0)
        MinS2Text.grid(row=6, column=1)
        MaxS2Text.grid(row=6, column=2)
        MinTauText.grid(row=6, column=4)
        MaxTauText.grid(row=6, column=5)
        MinOtherText.grid(row=6, column=7)
        MaxOtherText.grid(row=6, column=8)
        
        
        self.OP = ParamFile.Names['OrderParam']
        self.MinS2 = []
        self.MaxS2 = []
        self.mins2 = []
        self.maxs2 = []
        for i in range(len(self.OP)):
            s2Text = tk.Label(self.param, text=self.OP[i])
            self.mins2.append(tk.IntVar())
            self.mins2[-1].set("0.0")
            self.MinS2.append(tk.Entry(self.param, textvariable = self.mins2[-1], width=10))
            self.maxs2.append(tk.IntVar())
            self.maxs2[-1].set("1.0")
            self.MaxS2.append(tk.Entry(self.param, textvariable = self.maxs2[-1], width=10))
            
            s2Text.grid(row=7+i, column=0)
            self.MinS2[i].grid(row=7+i, column=1)
            self.MaxS2[i].grid(row=7+i, column=2)

        
        self.CT = ParamFile.Names['CorrTimes']
        self.MinTau = []
        self.MaxTau = []
        self.mintau = []
        self.maxtau = []
        for i in range(len(self.CT)):
            tauText = tk.Label(self.param, text=self.CT[i])
            self.mintau.append(tk.IntVar())
            self.mintau[-1].set("")
            self.MinTau.append(tk.Entry(self.param, textvariable = self.mintau[-1], width=10))
            self.maxtau.append(tk.IntVar())
            self.maxtau[-1].set("")
            self.MaxTau.append(tk.Entry(self.param, textvariable = self.maxtau[-1], width=10))
            
            tauText.grid(row=7+i, column=3)
            self.MinTau[i].grid(row=7+i, column=4)
            self.MaxTau[i].grid(row=7+i, column=5)
            
            
        self.Others = ParamFile.Names['others']
        self.MinOthers = []
        self.MaxOthers = []
        self.minOthers = []
        self.maxOthers = []
        for i in range(len(self.Others)):
            OthersText = tk.Label(self.param, text=self.Others[i])
            self.minOthers.append(tk.IntVar())
            self.minOthers[-1].set("")
            self.MinOthers.append(tk.Entry(self.param, textvariable = self.minOthers[-1], width=10))
            self.maxOthers.append(tk.IntVar())
            self.maxOthers[-1].set("")
            self.MaxOthers.append(tk.Entry(self.param, textvariable = self.maxOthers[-1], width=10))
            
            OthersText.grid(row=7+i, column=6)
            self.MinOthers[i].grid(row=7+i, column=7)
            self.MaxOthers[i].grid(row=7+i, column=8)
            
        incr = max(len(self.OP), len(self.CT), len(self.Others))
        self.nline = 6 + incr



#Load files part
        TextFiles = tk.Label(self.param, text = "Load files:", fg="blue")
        
        self.FieldCal = tk.Button(self.param, text="Field calibration file", command = lambda: BrowseFieldCal(self))
        self.ExpSetUp = tk.Button(self.param, text="Experimental setup file", command = lambda: BrowseExpSetUp(self))
        self.Irelaxom = tk.Button(self.param, text="Relaxometry Intensity folder", command = lambda: BrowseRelaxomInt(self))
        self.Input = tk.Button(self.param, text="Other inputs file", command = lambda: BrowseInputs(self))

        self.RelaxationDataSet = []
        self.Nlabel = 1
        self.RelaxationDataSet.append(tk.Button(self.param, text="High field rates", command = lambda: BrowseRates(self, ([0]))))

        
        self.MenuRelaxType = {}
        for i in range(len(ParamFile.RelaxationRates)):
            self.MenuRelaxType[str(ParamFile.RelaxationRates[i])] = i+1
        self.RelaxationDataType = []
        self.RelaxationType = []
        self.RelaxationDataType.append(tk.StringVar())
        self.RelaxationDataType[-1].set("Data type")
        self.RelaxationType.append(tk.OptionMenu(self.param, self.RelaxationDataType[-1], *self.MenuRelaxType.keys()))
        
        
        self.Fields = []
        Field = tk.StringVar()
        Field.set("High field (T)")
        self.Fields.append(tk.Entry(self.param, textvariable = Field, width=10))
        
        
        self.AddButton = tk.Button(self.param, text="Add high field rates", command = lambda: Add(self))
        
        TextFiles.grid(row=self.nline+1, column=0)
        self.FieldCal.grid(row=self.nline+2, column=0, columnspan=2)
        self.ExpSetUp.grid(row=self.nline+2, column=2)
        self.Irelaxom.grid(row=self.nline+3, column=0, columnspan=2)
        self.Input.grid(row=self.nline+3, column=2)
        
        for i in range(len(self.RelaxationDataSet)):
            self.RelaxationDataSet[i].grid(row=self.nline+4+i, column=0, columnspan=2)
            self.RelaxationType[i].grid(row=self.nline+4+i, column=2)
            self.Fields[i].grid(row=self.nline+4+i, column=3)
            
        
#PDB ID
        TextPDB = tk.Label(self.param, text = "4-letter PDB ID (if available):", fg = "blue")
        
        self.PDBID = tk.IntVar()
        self.PDBID.set("")
        self.PDBid = tk.Entry(self.param, textvariable = self.PDBID, width=10)
        
        TextPDB.grid(row=self.nline+1, column=6, columnspan=2)
        self.PDBid.grid(row=self.nline+2, column=6, columnspan=2)
                
        

    def GetData(self):

        GUI.begin = time.strftime("%H") + "H" + time.strftime("%M")
        
#################################################### Create the results folder and copy the input files in it ####################################################
        workingDirectiory = os.path.dirname(os.path.abspath(__file__))
        ResultDirectory = workingDirectiory + "/Results"

        if not os.path.exists(ResultDirectory):
            os.makedirs(ResultDirectory)
            
        GUI.directoryName = ResultDirectory + "/" + time.strftime("%Y-%m-%d")
        t = 0
        n = 2
        while t == 0:
            if not os.path.exists(GUI.directoryName):
                os.makedirs(GUI.directoryName)
                t = 1
            else:
                GUI.directoryName = ResultDirectory + "/" + time.strftime("%Y-%m-%d") + "_" + str(n)
                n += 1
                

        dirInput = GUI.directoryName + "/InputFiles"
        os.makedirs(dirInput)
                
        dirOutput = GUI.directoryName + "/FitAllResidues"
        os.makedirs(dirOutput)
        
        dirFitOutput = GUI.directoryName + "/FittingResults"
        os.makedirs(dirFitOutput)
        
        
        #Copy the directory functions
        fromDirectory = sys.argv[1]
        toDirectory = dirInput + "/ExpressionsAndConstraints"
        os.makedirs(toDirectory)
        copy_tree(fromDirectory, toDirectory)
        
        #Copy the input files
        ExperimentalSetUpPath = self.ExperimentalSetUp
        FieldCalibrationPath = self.FieldCalibration
        GUI.InputPath = self.InputFile

        
        dirHFRelax = dirInput + "/HFRelaxationRates"
        os.makedirs(dirHFRelax)
        for i in range(len(ParamFile.RelaxationRates)):
            os.makedirs(dirHFRelax + "/" + ParamFile.RelaxationRates[i])
            os.makedirs(dirOutput + "/" + ParamFile.RelaxationRates[i])
            
        for i in range(len(self.RATES)):
            RelaxFilePath = self.RATES[i]
            destFile = dirHFRelax + "/" + str(self.RelaxationDataType[i].get()) + "/" + str(os.path.basename(self.RATES[i]))
            copyfile(RelaxFilePath, destFile)
            

        FileNames = ["/ExpSetUp.txt", "/FieldCalibration.txt", "/OtherInputs.txt"]
        n = 0
        for i in [ExperimentalSetUpPath, FieldCalibrationPath, GUI.InputPath]:
            destFile = dirInput + FileNames[n]
            copyfile(i, destFile)
            n+=1

        
#################################################### Read the data ####################################################
        GUI.Nmcmc = int(self.Nmcmc.get())
        Nwalker_b = int(self.Nwalker.get())
        GUI.Nwalker = MCMC.nWalkerCheck(Nwalker_b, len(ParamFile.Names['OrderParam']) + len(ParamFile.Names['CorrTimes'])+len(ParamFile.Names['others'])+1)

        GUI.TauC = float(self.TauC.get())
            
        GUI.OP = self.OP
        GUI.CT = self.CT
        GUI.Others = self.Others
        GUI.TotParam = []
        for i in self.OP:
            GUI.TotParam.append(i)
        for i in self.CT:
            GUI.TotParam.append(i)
        for i in self.Others:
            GUI.TotParam.append(i)
                
        
        GUI.AccelerationType = str(self.AccelerationTYPE.get())
        
        GUI.PDB = self.PDBid.get()
        if len(GUI.PDB) == 4:
            GUI.checkPDB = True
        else:
            GUI.checkPDB = False

        
        ExperimentalSetUp = open(ExperimentalSetUpPath, "r")
        #this is a file containing the experiment number, the height (LF value), d22 (stab LF), d25 (stab HF), WTHF (response time), shuttle LF (time of shuttling to LF), WTLF, shuttle HF and VC (time spent at LF)
        
        FieldCalibration = open(FieldCalibrationPath, "r")
        #This is a file containing the height and the corresponding field
        
        with FieldCalibration as input:
            FC = list(zip(*(line.strip().split("\t") for line in input)))
        FC = [[float(FC[col][line]) for col in range(2)] for line in range(len(FC[0]))]   #it is a one column vector containing the height and the field at each position
        heights = [FC[col][0] for col in range(len(FC))]
        fields = [FC[col][1] for col in range(len(FC))]
        
        higherHeights = []
        lowerHeights = []
        middleHeights = []
        middleFields = []
        higherFields = []
        lowerFields = []
        for i in range(len(heights)):
            if heights[i] <= 0.47:
                higherHeights.append(heights[i])
                higherFields.append(fields[i])
            else:
                if heights[i] <= 0.6:
                    middleHeights.append(heights[i])
                    middleFields.append(fields[i])
                else:
                    lowerHeights.append(heights[i])
                    lowerFields.append(fields[i])

        GUI.HigherCoefs = np.polyfit(higherHeights, higherFields, 10)
        GUI.MiddleCoefs = np.polyfit(middleHeights, middleFields, 5)
        GUI.LowerCoefs = np.polyfit(lowerHeights, lowerFields, 10)


        
        #BO Files
        HFfiles = []
        HFtypes = []
        AllFields = []
        for i in range(len(self.RelaxationDataSet)):
            HFfiles.append(self.RATES[i])
            HFtypes.append(str(self.RelaxationDataType[i].get()))
            AllFields.append(float(self.Fields[i].get()))
        
        UsedFields = []
        for i in AllFields:
            if i not in UsedFields:
                UsedFields.append(i)
        
        GUI.B0HFields_int = []
        targetLength = len(UsedFields)
        while len(GUI.B0HFields_int) != targetLength:
            GUI.B0HFields_int.append(max(UsedFields))
            UsedFields.remove(max(UsedFields))

        

        numberExp = sum(1 for line in ExperimentalSetUp)
        ExperimentalSetUp.close()
        ExperimentalSetUp = open(ExperimentalSetUpPath, "r")
        
        GUI.ExperimentNumber = [[] for i in range(numberExp)]
        GUI.Height = [[] for i in range(numberExp)]
        GUI.d22 = [[] for i in range(numberExp)]
        GUI.d25 = [[] for i in range(numberExp)]
        GUI.WTHF = [[] for i in range(numberExp)]
        GUI.SLF = [[] for i in range(numberExp)]
        GUI.WTLF = [[] for i in range(numberExp)]
        GUI.SHF = [[] for i in range(numberExp)]
        GUI.VC = [[] for i in range(numberExp)]
        
        for i in range(numberExp):
            Line = ExperimentalSetUp.readline()
            Line = [x for x in Line.split('\n')]
            Line = Line[0]
            Line = [x for x in Line.split('\t')]
            Line = [float(col) for col in Line]
            
            GUI.ExperimentNumber[i] = int(Line[0])
            GUI.Height[i] = Line[1]
            GUI.d22[i] = Line[2] * 1e-3
            GUI.d25[i] = Line[3] * 1e-3
            GUI.WTHF[i] = Line[4] * 1e-3
            GUI.SLF[i] = Line[5] * 1e-3
            GUI.WTLF[i] = Line[6] * 1e-3
            GUI.SHF[i] = Line[7] * 1e-3
            
            vc = []
            for j in range(len(Line)-8):
                vc.append(Line[j+8] * 1e-3)
            GUI.VC[i] = vc
                
        GUI.B0LFields_Int = [[] for i in range(numberExp)]
        for i in range(len(GUI.Height)):
            GUI.B0LFields_Int[i] = FitF.B0Fit(GUI.Height[i], GUI.LowerCoefs, GUI.MiddleCoefs, GUI.HigherCoefs)
            
        GUI.MagField = FitF.B0Fit(0.0, GUI.LowerCoefs, GUI.MiddleCoefs, GUI.HigherCoefs)
        
        #Plot the field profile with the considered fields during relaxometry
        FigOut.PlotFieldProfile(heights, GUI.LowerCoefs, GUI.MiddleCoefs, GUI.HigherCoefs, GUI.B0LFields_Int, fields, dirFitOutput)
        

        #HF Relaxation rates Files
        HFfiles_orga = [[[] for j in GUI.B0HFields_int] for i in ParamFile.RelaxationRates]
        HFields_orga = [[[] for j in GUI.B0HFields_int] for i in ParamFile.RelaxationRates]
        for i in range(len(HFfiles)):
            posType = ParamFile.RelaxationRates.index(HFtypes[i])
            posField = GUI.B0HFields_int.index(AllFields[i])
            
            HFfiles_orga[posType][posField].append(HFfiles[i])
            HFields_orga[posType][posField].append(AllFields[i])
        #Amino Acids number list
        nAA = sum(1 for line in open(HFfiles[0]))
        
        GUI.AAList = []
        for line in open(HFfiles[0]):
            line = line.split("\n")
            line = line[0]
            line = line.split("\t")
            GUI.AAList.append(int(line[0]))
            
#create the average HFdata list
        GUI.HFdata = [[[] for i in ParamFile.RelaxationRates] for k in GUI.AAList]
        GUI.B0HFields = [[[] for i in ParamFile.RelaxationRates] for k in GUI.AAList]
        for RelaxType in range(len(ParamFile.RelaxationRates)):
            for HField in range(len(HFfiles_orga[RelaxType])):
                ndata = [0.0 for i in range(nAA)]
                for i in range(len(HFfiles_orga[RelaxType][HField])):
                    file = open(HFfiles_orga[RelaxType][HField][i], 'r')
                    with file as input:
                        Filedata = list(zip(*(line.strip().split("\t") for line in input)))
                        
                    for AA in range(nAA):
                        if Filedata[1][AA] != 'NA':
                            ndata[AA] += 1.0
                            if i == 0:
                                GUI.HFdata[AA][RelaxType].append([int(Filedata[0][AA]), float(Filedata[1][AA]), float(Filedata[2][AA])])
                                GUI.B0HFields[AA][RelaxType].append(float(HFields_orga[RelaxType][HField][0]))
                            else:
                                GUI.HFdata[AA][RelaxType][-1][1] = float(GUI.HFdata[AA][RelaxType][-1][1]) + float(Filedata[1][AA])
                                GUI.HFdata[AA][RelaxType][-1][2] = float(GUI.HFdata[AA][RelaxType][-1][2]) + float(Filedata[2][AA])
                    file.close()
                for AA in range(nAA):
                    if ndata[AA] != 0.0:
                        GUI.HFdata[AA][RelaxType][-1][1] = float(GUI.HFdata[AA][RelaxType][-1][1])/ndata[AA]
                        GUI.HFdata[AA][RelaxType][-1][2] = float(GUI.HFdata[AA][RelaxType][-1][2])/ndata[AA]
                        
                        

        #LF intensities
        dirIntRelax = dirInput + "/RelaxometryIntensities"
        os.makedirs(dirIntRelax)
        
        GUI.Intensities = [[] for AA in GUI.AAList]
        GUI.B0LFields = [[] for AA in GUI.AAList]
        for n in range(len(GUI.ExperimentNumber)):
            InterensityFile = self.Intrelax + "/" + str(int(GUI.ExperimentNumber[n])) + ".txt"
            if os.path.isfile(InterensityFile):
                destFile = dirIntRelax + "/" + str(int(GUI.ExperimentNumber[n])) + ".txt"
                copyfile(InterensityFile, destFile)
                
                AA = 0
                for line in open(InterensityFile):
                    L = line.split("\n")
                    L = L[0]
                    L = L.split("\t")
                    if int(2*len(GUI.VC[n])+1) != len(L):
                        print("")
                        print("Eperiment number " + str(int(GUI.ExperimentNumber[n])) + " does not have the correct number of intensities")
                        print("")
                        sys.exit()
                    else:
                        ToAdd = []
                        for vc in range(int((len(L)-1)/2)):
                            if L[vc*2+1] != "NA":
                                ToAdd.append([int(L[0]), float(L[vc*2+1]), float(L[vc*2+2])])
                            else:
                                ToAdd.append([int(L[0]), "NA", "NA"])
                        if len(ToAdd) != 0:
                            GUI.B0LFields[AA].append(FitF.B0Fit(GUI.Height[n], GUI.LowerCoefs, GUI.MiddleCoefs, GUI.HigherCoefs))
                            GUI.Intensities[AA].append(ToAdd)
                        AA += 1
            else:
                print("")
                print("Missing the eperiment number " + str(int(GUI.ExperimentNumber[n])) + " file")
                print("")
                sys.exit()
        
                
                
        
#Other inputs
        InputF = open(GUI.InputPath, "r")
        with InputF as input:
            inputL = list(zip(*(line.strip().split('\t') for line in input)))
        
        GUI.OtherInputs = [[] for AA in GUI.AAList]
        for AA in range(len(GUI.AAList)):
            listOtherInputs = []
            for i in range(len(inputL)-1):
                listOtherInputs.append(float(inputL[1+i][AA]))
            GUI.OtherInputs[AA] = [float(inputL[0][AA]), listOtherInputs]
            
            
                        
#Write the Parameters file
        loadFile = dirInput + "/Parameters.txt"
        GUI.bnds = Out.writeParam(self, loadFile, GUI.TauC, GUI.AccelerationType, GUI.Nmcmc, GUI.Nwalker, GUI.checkPDB, GUI.PDB, ExperimentalSetUpPath, FieldCalibrationPath, self.Intrelax, GUI.InputPath)
        
        
        
        self.PositionR1 = 0
        for i in range(len(self.RelaxFunc)):
            if ParamFile.RelaxationRates[i] != "R1":
                self.PositionR1 += 1
                break
        
        
        GUI.PositionR1 = self.PositionR1
        GUI.RelaxFunc = self.RelaxFunc

        Calculations()


class Calculations(GUI):
    def __init__(self):
        

            self.FieldLists()
            
            
    def FieldLists(self):
#################################################### Preparation before MCMC ####################################################
        print("")
        print("Optimizing shuttling time increment...")
       
        self.LFtimes = [[self.d22[wtlf] + self.WTLF[wtlf] + self.VC[wtlf][vc] for vc in range(len(self.VC[wtlf]))] for wtlf in range(len(self.WTLF))]
        posLowField = self.B0LFields[0].index(min(self.B0LFields[0]))
        
        RandomParam = [[] for i in range(10)]
        for i in range(10):
            for P in range(len(self.TotParam)):
                RandomParam[i].append(uniform(self.bnds[P][0], self.bnds[P][1]))
                
                
        self.Increment = ShSim.optShuttling(self, posLowField, np.array(RandomParam), ParamFile.PositionAuto)
        self.FieldListUp, self.FieldListDown = ShSim.FieldList(self, self.Increment)
        print("Final used increment: ", self.Increment, " s")
        
        
        
        
        print("")
        print("Choosing propagator calculation method...")
        refRM = np.asarray(_RelaxMat.RelaxMat(5.0, RandomParam[0], self.TauC, self.OtherInputs[0][1])[0])
        

        print(" Method 1")
        start = time.time()
        for i in range(10000):
            linalg.expm(-refRM)
        end = time.time()
        Duration1 = end-start
        print("  time for 10,000 iterations: ", round(Duration1, 1), " s")
        
        print(" Method 2")
        start = time.time()
        for i in range(10000):
            eig, eigvec = np.linalg.eig(refRM)
            eigvec @ np.diag(np.exp(-eig)) @ np.linalg.inv(eigvec)
        end = time.time()
        Duration2 = end-start
        print("  time for 10,000 iterations: ", round(Duration2, 1), " s")
        
        if min(Duration1, Duration2) == Duration1:
            self.PropFunction = ShSim.PropCalExp
            print("Choosing calculation method done. Method 1 chosen.")
        else:
            self.PropFunction = ShSim.PropCalDiag
            print("Choosing calculation method done. Method 2 chosen.")
        
        self.MCMCcalculations()
        
        
        
    def MCMCcalculations(self):
        FigMCMCcorrFolder = self.directoryName + "/FittingResults/Correlations"
        FigMCMCtrajFolder = self.directoryName + "/FittingResults/Trajectories"
        FigIntensities  = self.directoryName + "/FittingResults/Intensities"
        os.makedirs(FigMCMCcorrFolder)
        os.makedirs(FigMCMCtrajFolder)
        os.makedirs(FigIntensities)

        self.MCMCparam = [[] for AA in self.AAList]
        self.Acceptance = [[] for AA in self.AAList]
        self.FinalSimulatedIntensities = [[] for AA in self.AAList]
        
        self.R1LFDataForCurve_BackCalc = [[] for i in self.AAList]
        self.R1LFDataForCurve_Fitted = [[] for i in self.AAList]
        self.ScalingIntensities = [[] for i in self.AAList]
        
        nParam = len(ParamFile.Names['OrderParam']) + len(ParamFile.Names['CorrTimes'])+len(ParamFile.Names['others'])+1
        
        
        
        
        print("")
        print("Monte Carlo")
        
        
        for AA in range(len(self.AAList)):
            print("")
            print(" Residue ", self.AAList[AA])
            print("")
            
            
            self.MCMCparam[AA], self.Acceptance[AA], FullTraj = MCMC.MarkovChainMonteCarlo(FigMCMCcorrFolder, FigMCMCtrajFolder, self.Intensities[AA], self.HFdata[AA], self.B0HFields[AA], self.TauC, self.OtherInputs[AA][1], self.MagField, self.Increment, self.FieldListUp, self.FieldListDown, self.ExperimentNumber, self.WTHF, self.d25, self.LFtimes, self.B0LFields[AA], self.Nwalker, self.Nmcmc, self.bnds, self.TotParam, self.AAList[AA], self.PropFunction, nParam)
                
            self.FinalSimulatedIntensities[AA] = self.PropFunction(np.array(self.MCMCparam[AA][0][:-1]), self.TauC, np.array(self.OtherInputs[AA][1]), self.MagField, self.Increment, self.FieldListUp, self.FieldListDown, self.ExperimentNumber, self.WTHF, self.d25, self.LFtimes, self.B0LFields[AA], ParamFile.PositionAuto)
    
            Out.WriteMCMCTraj(self, FullTraj, self.AAList[AA])         #file containing the MCMC trajectories
                

                
            print("")
            print(" Making figures")
            print("")
            
            #Draw the intensities
            self.ScalingIntensities[AA] = MCMC.ScalingFactor(self.FinalSimulatedIntensities[AA], self.Intensities[AA])
            
            timeForSim = [np.linspace(min(self.LFtimes[LField])-min(self.VC[LField]), max(self.LFtimes[LField]), 100) for LField in range(len(self.VC))]
            timeForSim_Plot = [np.linspace(0.0, max(self.VC[LField]), 100) for LField in range(len(self.VC))]
            FinalSimulatedIntensitiesFull = self.PropFunction(np.array(self.MCMCparam[AA][0][:-1]), self.TauC, np.array(self.OtherInputs[AA][1]), self.MagField, self.Increment, self.FieldListUp, self.FieldListDown, self.ExperimentNumber, self.WTHF, self.d25, timeForSim, self.B0LFields[AA], ParamFile.PositionAuto)
            
            BackIntensities = [[self.ScalingIntensities[AA][LField]*FinalSimulatedIntensitiesFull[LField][t] for t in range(len(FinalSimulatedIntensitiesFull[LField]))] for LField in range(len(self.B0LFields[AA]))]
                
            IntensitiesFigFolder = FigIntensities + "/Residue" + str(self.AAList[AA])
            os.makedirs(IntensitiesFigFolder)
                
                
            IntensitiesForPlot = [[] for LF in range(len(self.Intensities[AA]))]
            IntensitiesErrForPlot = [[] for LF in range(len(self.Intensities[AA]))]
            DelaysForPlot = [[] for LF in range(len(self.Intensities[AA]))]
            for LF in range(len(self.Intensities[AA])):
                for VC in range(len(self.Intensities[AA][LF])):
                    if self.Intensities[AA][LF][VC][1] != "NA":
                        IntensitiesForPlot[LF].append(self.Intensities[AA][LF][VC][1])
                        IntensitiesErrForPlot[LF].append(self.Intensities[AA][LF][VC][2])
                        DelaysForPlot[LF].append(self.VC[LF][VC])
                
            FigOut.PlotIntensities(self, BackIntensities, IntensitiesFigFolder, IntensitiesForPlot, IntensitiesErrForPlot, DelaysForPlot, timeForSim_Plot, AA)
                
            #Draw the relaxation rates 
            for LField in range(len(self.B0LFields[AA])):
                self.R1LFDataForCurve_BackCalc[AA].append(self.RelaxFunc[self.PositionR1](self.B0LFields[AA][LField], self.MCMCparam[AA][0][:-1], self.TauC, self.OtherInputs[AA][1])[0])
                ParamOpt, ParamCov = curve_fit(FitF.exp, np.asarray(self.VC[LField]), np.asarray(self.Intensities[AA][LField])[:,1])
                
                self.R1LFDataForCurve_Fitted[AA].append([self.B0LFields[AA][LField], ParamOpt[0]])
                        
                        
            yRateMesHF = [[] for i in ParamFile.RelaxationRates]
            yRateErrMesHF = [[] for i in ParamFile.RelaxationRates]
            ResiHF = [[] for i in ParamFile.RelaxationRates]
            xB0HF = [[] for i in ParamFile.RelaxationRates]
            for RelaxRate in range(len(ParamFile.RelaxationRates)):
                for HField in range(len(self.B0HFields[AA][RelaxRate])):
                    yRateMesHF[RelaxRate].append(self.HFdata[AA][RelaxRate][HField][1])
                    yRateErrMesHF[RelaxRate].append(self.HFdata[AA][RelaxRate][HField][2])
                    xB0HF[RelaxRate].append(self.B0HFields[AA][RelaxRate][HField])
                    ResiHF[RelaxRate].append(self.HFdata[AA][RelaxRate][HField][1] - self.RelaxFunc[RelaxRate](self.B0HFields[AA][RelaxRate][HField], self.MCMCparam[AA][0][:-1], self.TauC, self.OtherInputs[AA][1])[0])
                        
                            
            minLF = min(self.B0LFields[AA])
            maxHF = 25.0
            xFields = np.logspace(np.log(minLF)/np.log(10.), np.log(maxHF)/np.log(10.), 100)
            xFieldsHF = np.linspace(8., maxHF, 100)
            for RelaxRate in range(len(ParamFile.RelaxationRates)):
                if RelaxRate == self.PositionR1:
                    yback = [self.RelaxFunc[self.PositionR1](B0, self.MCMCparam[AA][0][:-1], self.TauC, self.OtherInputs[AA][1])[0] for B0 in xFields]
                    FigOut.PlotR1(self, xB0HF[self.PositionR1], ResiHF[self.PositionR1], yRateMesHF[self.PositionR1], yRateErrMesHF[self.PositionR1], yback, xFields, AA)
                else:
                    yback = [self.RelaxFunc[RelaxRate](B0, self.MCMCparam[AA][0][:-1], self.TauC, self.OtherInputs[AA][1])[0] for B0 in xFieldsHF]
                    FigOut.PlotRate(self, xB0HF[RelaxRate], ResiHF[RelaxRate], yRateMesHF[RelaxRate], yRateErrMesHF[RelaxRate], ParamFile.RelaxationRates[RelaxRate], yback, xFieldsHF, AA)
                
            print("")
                

        self.WriteResult()

    def WriteResult(self):
        LAAList = len(self.AAList)
        
        print("")
        print("Writing final results...")
        
        
        
        #Put together figures already done
        F = self.directoryName + "/FittingResults"
        if len(self.AAList) > 1:
            FigOut.Convert(F + "/Correlations", "AllCorrelations.pdf", "png", False)
            FigOut.Convert(F + "/Trajectories", "AllTrajectories.pdf", "png", False)
        
        for AA in range(len(self.AAList)):
            F1 = F + "/Intensities/Residue" + str(self.AAList[AA])
            FigOut.Convert(F1, "AllDecays_Residue" + str(self.AAList[AA]) + ".pdf", "png", True)
            
        
        dirFigs = self.directoryName + "/PlotParameters"
        os.makedirs(dirFigs)
        
        Out.WriteMCMCParam(self)        #file containing the parameters of the spectral density function extracted from the MCCM

        #Draw the Chi2
        AllChi2 = [[] for i in self.AAList]
        for AA in range(LAAList):
            AllChi2[AA] = FitF.Chi2TOT(self.MCMCparam[AA][0][:-1], self.FinalSimulatedIntensities[AA], self.Intensities[AA], self.HFdata[AA], self.B0HFields[AA], self.TauC, self.OtherInputs[AA][1])

        FigOut.PlotChi2(dirFigs, AllChi2, self.AAList)

        #Draw the spetral density function parameters
        for param in range(len(self.TotParam)):
            paramForPlot = [self.MCMCparam[AA][0][param] for AA in range(len(self.AAList))]
            ErrForPlot = [(self.MCMCparam[AA][1][param]+self.MCMCparam[AA][2][param])/2.0 for AA in range(len(self.AAList))]
            FigOut.PlotDynParam(dirFigs, paramForPlot, ErrForPlot, self.TotParam[param], self.AAList)
            
            

        #Write the LF R1
        Out.WriteLFR1(self)        #file containing the scaling factors for intensities, back-calculated and fitted low field R1
             


        if len(self.AAList) > 1:
            F = self.directoryName + "/FitAllResidues/"
            for Rate in range(len(ParamFile.RelaxationRates)):
                F1 = F + str(ParamFile.RelaxationRates[Rate])
                FigOut.Convert(F1, "All" + str(ParamFile.RelaxationRates[Rate]) + ".pdf", "png", True)


                        
#Write the PDB files
        if self.checkPDB:
            Out.WritePDB(self, AllChi2)
                

        print("    Writing results: Done")
        print("")
        
        end = time.strftime("%H") + "H" + time.strftime("%M")
        print("Started at: " + self.begin)
        print("Ended at: " + end)
        
        sys.exit(0)
        
    
if __name__ == "__main__":
    root = tk.Tk()
    GUI(root)
    root.mainloop()
    
    
    
    
    