#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tkinter as tk
from tkinter.filedialog import askopenfilename, askdirectory
from PIL import ImageTk, Image
from shutil import copyfile
from distutils.dir_util import copy_tree
from random import uniform
from scipy import linalg
from matplotlib.backends.backend_pdf import PdfPages

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
import ReadInput


os.environ["OMP_NUM_THREADS"] = "1"


#################################################### Graphical User Interface ####################################################
class GUI(tk.Frame):
    def __init__(self, parent):
        
        tk.Frame.__init__(self, parent)
        
        self.parent = parent
        
        self.RelaxFunc = ParamFile.ImportFunc()
        self.RATES = []
        self.browsing_dir = os.path.dirname(os.path.abspath(__file__))
        
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
            self.FieldCalibration = askopenfilename(initialdir = self.browsing_dir)
            self.FieldCal.config(text=os.path.basename(self.FieldCalibration))
            self.browsing_dir = self.FieldCalibration
            
        def BrowseExpSetUp(self):
            self.ExperimentalSetUp = askopenfilename(initialdir = self.browsing_dir)
            self.ExpSetUp.config(text=os.path.basename(self.ExperimentalSetUp))
            self.browsing_dir = self.ExperimentalSetUp
            
        def BrowseRelaxomInt(self):
            self.Intrelax = askdirectory(initialdir = self.browsing_dir)
            self.Irelaxom.config(text=os.path.basename(self.Intrelax))
            self.browsing_dir = self.Intrelax
            
        def BrowseInputs(self):
            self.InputFile = askopenfilename(initialdir = self.browsing_dir)
            self.Input.config(text=os.path.basename(self.InputFile))
            self.browsing_dir = self.InputFile
            
        def BrowseRates(self, n):
            n = n[0]
            if n+1 > len(self.RATES):
                self.RATES.append(askopenfilename(initialdir = self.browsing_dir))
            else:
                self.RATES[n] = askopenfilename(initialdir = self.browsing_dir)
                
            self.RelaxationDataSet[n].config(text=os.path.basename(self.RATES[n]))
            self.browsing_dir = self.RATES[n]
            
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
            self.parameters = askopenfilename(initialdir = self.browsing_dir)
            
            UpdateGUI(self, self.parameters)
        
        def UpdateGUI(self, path):
            nparam = sum(1 for line in open(path, 'r'))-1
            paramFile = open(path, 'r')
            
            paramFile.readline()
    
            previousSetUp = []
            for i in range(nparam):
                Line = paramFile.readline()
                Line = [x for x in Line.split('\n')]
                Line = Line[0]
                Line = [x for x in Line.split('\t')]
                previousSetUp.append(Line)
                
            self.TAUC.set(previousSetUp[0][1])
                    
            self.Shuttling_TYPE.set(previousSetUp[1][1])
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
                
            #if the user indicated it in the command line, immediately start
            if len(sys.argv) > 3 and str(sys.argv[3]) == '-start':
                    self.GetData()
            
            
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

        self.Shuttling_TYPE = tk.StringVar()
        self.Shuttling_TYPE.set("Type of shuttling")
        Shuttling_Type = tk.OptionMenu(self.param, self.Shuttling_TYPE, "Constant Speed", "Constant Acceleration")
        
        self.NMCMC= tk.IntVar()
        self.NMCMC.set("MCMC - Number of steps")
        self.Nmcmc = tk.Entry(self.param, textvariable = self.NMCMC, width=25)
        
        self.NWALKER = tk.IntVar()
        self.NWALKER.set("MCMC - Number of chains")
        self.Nwalker = tk.Entry(self.param, textvariable = self.NWALKER, width=25)
        
        
        ParamText.grid(row=1, column=0)
        self.TauC.grid(row=2, column=0, columnspan=3)
        
        Shuttling_Type.grid(row=2, column=3, columnspan=2)

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
        
        
        #option to indicate a load-file in the command line
        if len(sys.argv) > 2:
            UpdateGUI(self, str(sys.argv[2]))
        

    def GetData(self):

        GUI.begin = time.strftime("%H") + "H" + time.strftime("%M")
        
#################################################### Create the results folder and copy the input files in it ####################################################
        result_directory = os.path.dirname(os.path.abspath(__file__)) + "/Results"

        if not os.path.exists(result_directory):
            os.makedirs(result_directory)
            
        GUI.directory_name = result_directory + "/" + time.strftime("%Y-%m-%d")
        n = 2
        while True:
            if not os.path.exists(GUI.directory_name):
                os.makedirs(GUI.directory_name)
                break
            else:
                GUI.directory_name = result_directory + "/" + time.strftime("%Y-%m-%d") + "_" + str(n)
                n += 1

        dirInput = GUI.directory_name + "/InputFiles"
        os.makedirs(dirInput)
                
        dirOutput = GUI.directory_name + "/FitAllResidues"
        os.makedirs(dirOutput)
        
        dirFitOutput = GUI.directory_name + "/FittingResults"
        os.makedirs(dirFitOutput)
        
        #Copy the directory functions
        fromDirectory = sys.argv[1]
        toDirectory = dirInput + "/ExpressionsAndConstraints"
        os.makedirs(toDirectory)
        copy_tree(fromDirectory, toDirectory)
        
        #Copy the input files
        GUI.InputPath = self.InputFile
        
        dirHFRelax = dirInput + "/HFRelaxationRates"
        os.makedirs(dirHFRelax)
        for i in range(len(ParamFile.RelaxationRates)):
            os.makedirs(dirHFRelax + "/" + ParamFile.RelaxationRates[i])
            
        for i in range(len(self.RATES)):
            RelaxFilePath = self.RATES[i]
            destFile = dirHFRelax + "/" + str(self.RelaxationDataType[i].get()) + "/" + str(os.path.basename(self.RATES[i]))
            copyfile(RelaxFilePath, destFile)

        file_names = {}
        file_names["ExpSetUp.txt"] = self.ExperimentalSetUp
        file_names["FieldCalibration.txt"] = self.FieldCalibration
        file_names["OtherInputs.txt"] = self.InputFile
        for f in file_names.keys():
            destFile = dirInput + "/" + f
            copyfile(file_names[f], destFile)

        
#################################################### Read the data ####################################################
        GUI.Nmcmc = int(self.Nmcmc.get())
        Nwalker_input = int(self.Nwalker.get())
        GUI.Nwalker = MCMC.nWalkerCheck(Nwalker_input, len(ParamFile.Names['OrderParam']) + len(ParamFile.Names['CorrTimes'])+len(ParamFile.Names['others'])+1)

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
        
        GUI.Shuttling_type = str(self.Shuttling_TYPE.get())
        
        GUI.PDB = self.PDBid.get()
        if len(GUI.PDB) == 4:
            GUI.checkPDB = True
        else:
            GUI.checkPDB = False
        

        #Field calibration
        field_cal = ReadInput.Read_Field_Calibration(self.FieldCalibration)
        
        field_cal_high = {}
        field_cal_middle = {}
        field_cal_low = {}
        for h in field_cal.keys():
            if h <= 0.47:
                field_cal_high[h] = field_cal[h]
            else:
                if h <= 0.6:
                    field_cal_middle[h] = field_cal[h]
                else:
                    field_cal_low[h] = field_cal[h]

        GUI.HigherCoefs = np.polyfit(list(field_cal_high.keys()), list(field_cal_high.values()), 10)
        GUI.MiddleCoefs = np.polyfit(list(field_cal_middle.keys()), list(field_cal_middle.values()), 5)
        GUI.LowerCoefs = np.polyfit(list(field_cal_low.keys()), list(field_cal_low.values()), 10)

        GUI.MagField = FitF.B0Fit(0.0, GUI.LowerCoefs, GUI.MiddleCoefs, GUI.HigherCoefs)

        #experimental set up
        GUI.set_up = ReadInput.Read_Exp_Setup(self.ExperimentalSetUp)
        #this is a file containing the experiment number, the height (LF value), d22 (stab LF), d25 (stab HF), WTHF (response time), shuttle LF (time of shuttling to LF), WTLF, shuttle HF and VC (time spent at LF)

        #Plot the field profile with the considered fields during relaxometry
        GUI.B0LFields = {}
        for exp in GUI.set_up.keys():
            GUI.B0LFields[exp] = FitF.B0Fit(GUI.set_up[exp]['height'], GUI.LowerCoefs, GUI.MiddleCoefs, GUI.HigherCoefs)
        FigOut.PlotFieldProfile(field_cal, GUI.LowerCoefs, GUI.MiddleCoefs, GUI.HigherCoefs, GUI.B0LFields, dirFitOutput)
        
        
        #high field Files
        HF_data_all = {}
        HF_err_all = {}
        for RelaxType in ParamFile.RelaxationRates:
            HF_data_all[RelaxType] = {}
            HF_err_all[RelaxType] = {}
            
        for dataset in range(len(self.RelaxationDataSet)):
            relax_type = str(self.RelaxationDataType[dataset].get())
            b0 = float(self.Fields[dataset].get())
                
            if b0 not in HF_data_all[relax_type].keys():
                HF_data_all[relax_type][b0] = []
                HF_err_all[relax_type][b0] = []
                    
            data, err = ReadInput.Read_High_Field_Data(self.RATES[dataset])
            HF_data_all[relax_type][b0].append(data)
            HF_err_all[relax_type][b0].append(err)        
            
        #create the average HF data
        GUI.HF_data = {}
        GUI.HF_err = {}
        for RelaxType in HF_data_all.keys():
            GUI.HF_data[RelaxType] = {}
            GUI.HF_err[RelaxType] = {}
            
            for field in HF_data_all[RelaxType].keys():
                GUI.HF_data[RelaxType][field] = {}
                GUI.HF_err[RelaxType][field] = {}
                
                for aa in HF_data_all[RelaxType][field][0].keys():
                    GUI.HF_data[RelaxType][field][aa] = 0.
                    GUI.HF_err[RelaxType][field][aa] = 0.
                    
                    n_repeat = len(HF_data_all[RelaxType][field])
                    for repeat in range(n_repeat):
                        GUI.HF_data[RelaxType][field][aa] += HF_data_all[RelaxType][field][repeat][aa] / n_repeat
                        GUI.HF_err[RelaxType][field][aa] += HF_err_all[RelaxType][field][repeat][aa]**2
                    GUI.HF_err[RelaxType][field][aa] = np.sqrt(GUI.HF_err[RelaxType][field][aa]) / n_repeat                        
                        

        #LF intensities
        dirIntRelax = dirInput + "/RelaxometryIntensities"
        os.makedirs(dirIntRelax)
        
        valid_file_ext = ['txt', 'dat', 'csv']
        
        GUI.Intensities = {}
        GUI.Err_Int = {}
        for exp in GUI.set_up.keys():
            file_found = False
            for extension in valid_file_ext:
                InterensityFile = f'{self.Intrelax}/{exp}.{extension}'
                if os.path.isfile(InterensityFile):
                    destFile = f'{dirIntRelax}/{exp}.{extension}'
                    copyfile(InterensityFile, destFile)
                    file_found = True
                    break
                
            if not file_found:
                print()
                print(f"Missing the eperiment number {exp} file")
                print()
                sys.exit()
                
            GUI.Intensities[exp], GUI.Err_Int[exp] = ReadInput.Read_Relaxometry_Decay(InterensityFile, GUI.set_up[exp]['vc'], exp)
                
        
#Other inputs
        GUI.OtherInputs = ReadInput.Read_Other_Input(GUI.InputPath)
            
                        
#Write the Parameters file
        loadFile = dirInput + "/Parameters.txt"
        GUI.bnds = Out.writeParam(self, loadFile, GUI.TauC, GUI.Shuttling_type, GUI.Nmcmc, GUI.Nwalker, GUI.checkPDB, GUI.PDB, self.ExperimentalSetUp, self.FieldCalibration,
                                  self.Intrelax, GUI.InputPath)
        
        self.PositionR1 = np.where(np.asarray(ParamFile.RelaxationRates) == "R1")[0][0]
        
        GUI.PositionR1 = self.PositionR1
        GUI.RelaxFunc = self.RelaxFunc


        Calculations()


class Calculations(GUI):
    def __init__(self):

            self.FieldLists()
            
    def FieldLists(self):
#################################################### Preparation before MCMC ####################################################
        print()
        print("Optimizing shuttling increment...")
        lowest_field_idx = list(self.B0LFields.values()).index(min(list(self.B0LFields.values())))
        exp_lowest_field = list(self.B0LFields.keys())[lowest_field_idx]
        
        RandomParam = [[] for i in range(10)]
        for i in range(10):
            for P in range(len(self.TotParam)):
                RandomParam[i].append(uniform(self.bnds[P][0], self.bnds[P][1]))
                
        self.Increment = ShSim.optShuttling(self, exp_lowest_field, np.array(RandomParam), ParamFile.PositionAuto)
        self.shuttling_fields, self.shuttling_delays = ShSim.FieldList(self, self.Increment)
        print("Final used increment: ", self.Increment, " m")
        
        print()
        print("Choosing propagator calculation method...")
        aa = list(self.OtherInputs.keys())[0]
        relax_mat = np.asarray(_RelaxMat.RelaxMat(5.0, RandomParam[0], self.TauC, self.OtherInputs[aa]))
        dt = 1e-3

        print(" Method 1")
        start = time.time()
        for i in range(10000):
            linalg.expm(-dt * relax_mat)
        end = time.time()
        Duration1 = end-start
        print(f"  time for 10,000 iterations: {round(Duration1, 1)} s")
        
        print(" Method 2")
        start = time.time()
        for i in range(10000):
            eig, eigvec = np.linalg.eig(relax_mat)
            eigvec @ np.diag(np.exp(-dt * eig)) @ np.linalg.inv(eigvec)
        end = time.time()
        Duration2 = end-start
        print(f"  time for 10,000 iterations: {round(Duration2, 1)} s")
        
        if min(Duration1, Duration2) == Duration1:
            self.PropFunction = ShSim.PropCalExp
            print("Choosing calculation method done. Method 1 chosen.")
        else:
            self.PropFunction = ShSim.PropCalDiag
            print("Choosing calculation method done. Method 2 chosen.")
            
            
            
        # param = np.asarray([0.9, 62475078762.51121, 8.811421637299617, 0.8890747323612564])
        # calc_diag = ShSim.PropCalDiag_print(param, 5e-9, [], 14.1, 6e-3, self.shuttling_fields, self.shuttling_delays, self.set_up, self.B0LFields, ParamFile.PositionAuto)
        # calc_exp = ShSim.PropCalDiag(param, 5e-9, [], 14.1, 6e-3, self.shuttling_fields, self.shuttling_delays, self.set_up, self.B0LFields, ParamFile.PositionAuto)
        # calc_diag = ShSim.PropCalDiag(param, 5e-9, [], 14.1, 6e-3, self.shuttling_fields, self.shuttling_delays, self.set_up, self.B0LFields, ParamFile.PositionAuto)
        
        # exp = list(self.set_up.keys())[-1]
        
        # print(exp)
        # print(calc_exp[exp])
        # print(calc_diag[exp]) 

        # print()
        # print()
        # print(self.shuttling_delays['up'][exp])
        # print(self.shuttling_delays['down'][exp])
        # sys.exit()
        
        self.MCMCcalculations()
        
        
    def MCMCcalculations(self):
        FigMCMCcorrFolder = self.directory_name + "/FittingResults/Correlations"
        FigMCMCtrajFolder = self.directory_name + "/FittingResults/Trajectories"
        FigIntensities  = self.directory_name + "/FitAllResidues/Intensities"
        os.makedirs(FigMCMCcorrFolder)
        os.makedirs(FigMCMCtrajFolder)
        os.makedirs(FigIntensities)
        
        pdf_trajectories = PdfPages(f'{FigMCMCtrajFolder}/All_trajectories.pdf' )
        pdf_correlations = PdfPages(f'{FigMCMCcorrFolder}/All_correlations.pdf' )
        
        pdf_rates = {}
        for RelaxRate in ParamFile.RelaxationRates:
            pdf_rates[RelaxRate] = PdfPages(f'{self.directory_name}/FitAllResidues/{RelaxRate}.pdf')

        self.MCMCparam = {}
        self.Acceptance = {}
        self.FinalSimulatedIntensities = {}
        
        self.R1LF_BackCalc = {}
        self.R1LF_Fitted = {}
        self.ScalingIntensities = {}
        
        nParam = len(ParamFile.Names['OrderParam']) + len(ParamFile.Names['CorrTimes'])+len(ParamFile.Names['others'])+1

        e = list(self.Intensities.keys())[0]
        AAList = list(self.Intensities[e].keys())
        
        print()
        print("Monte Carlo")
        
        for AA in AAList:
            print()
            print(f" Residue {AA}")
            print()
            
            self.MCMCparam[AA], self.Acceptance[AA], FullTraj = MCMC.MarkovChainMonteCarlo(pdf_trajectories, pdf_correlations, AA, nParam, FigMCMCcorrFolder, FigMCMCtrajFolder,
                                                                self.TauC, self.MagField, self.Nwalker, self.Nmcmc, self.bnds, self.TotParam,
                                                                self.Increment, self.PropFunction, self.shuttling_fields, self.shuttling_delays,
                                                                self.Intensities, self.Err_Int, self.HF_data, self.HF_err, self.OtherInputs, self.set_up, self.B0LFields)
                
            self.FinalSimulatedIntensities[AA] = self.PropFunction(np.asarray(self.MCMCparam[AA][0][:-1]), self.TauC, self.OtherInputs[AA], self.MagField,
                                                                   self.Increment, self.shuttling_fields, self.shuttling_delays, self.set_up,
                                                                   self.B0LFields, ParamFile.PositionAuto)
    
            Out.WriteMCMCTraj(self, FullTraj, AA)         #file containing the MCMC trajectories
                
            print()
            print(" Making figures")
            print()
            
            self.ScalingIntensities[AA] = MCMC.ScalingFactor(self.FinalSimulatedIntensities[AA], self.Intensities, self.Err_Int, AA)
            self.R1LF_BackCalc[AA] = FitF.CalcR1_LF(self.MCMCparam[AA][0][:-1], self.TauC, self.OtherInputs[AA], self.B0LFields)
            self.R1LF_Fitted[AA] = FitF.FitR1_LF(self.Intensities, AA)
            
            #Draw the intensities
            FigOut.PlotIntensities(self, FigIntensities, AA)
                
            #Draw the relaxation rates 
            for RelaxRate in ParamFile.RelaxationRates:
                if RelaxRate == 'R1':
                    FigOut.PlotR1(self, pdf_rates['R1'], AA)
                else:
                    FigOut.PlotRate(self, pdf_rates[RelaxRate], RelaxRate, AA)
                    
            print("")
                
        #close all the pdf figures
        for RelaxRate in ParamFile.RelaxationRates:
            pdf_rates[RelaxRate].close()
        pdf_trajectories.close()
        pdf_correlations.close()
            
            
        self.WriteResult()

    def WriteResult(self):
        r = list(self.HF_data.keys())[0]
        b = list(self.HF_data[r].keys())[0]
        AAList = list(self.HF_data[r][b].keys())
        
        print("")
        print("Writing final results...")
        
        dirFigs = self.directory_name + "/PlotParameters"
        os.makedirs(dirFigs)
        
        Out.WriteMCMCParam(self, AAList)        #file containing the parameters of the spectral density function extracted from the MCMC

        #Draw the Chi2
        AllChi2 = {}
        for AA in AAList:
            AllChi2[AA] = FitF.Chi2TOT(self.MCMCparam[AA][0][:-1], self.FinalSimulatedIntensities[AA], self.ScalingIntensities[AA], self.Intensities,
                                       self.Err_Int, self.HF_data, self.HF_err, self.TauC, self.OtherInputs[AA], AA)
        FigOut.PlotChi2(dirFigs, AllChi2)

        #Draw the spectral density function parameters
        for count, param in enumerate(self.TotParam):
            FigOut.PlotDynParam(dirFigs, self.MCMCparam, param, count)

        #Write the LF R1
        Out.WriteLFR1(self)        #file containing the scaling factors for intensities, back-calculated and fitted low field R1

        #Write the PDB files
        if self.checkPDB:
            Out.WritePDB(self, AllChi2)

        print("    Writing results: Done")
        print("")
        
        try:
            os.system('rm *.so')
        except:
            pass
        
        end = time.strftime("%H") + "H" + time.strftime("%M")
        print("Started at: " + self.begin)
        print("Ended at: " + end)
        
        sys.exit(0)
        
    
if __name__ == "__main__":
    root = tk.Tk()
    GUI(root)
    root.mainloop()
