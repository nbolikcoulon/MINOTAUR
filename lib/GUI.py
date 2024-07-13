#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################################
#                                                 #
#                       GUI                       #
#                                                 #
###################################################
import numpy as np
import tkinter as tk
from PIL import ImageTk, Image
from tkinter.filedialog import askopenfilename, askdirectory

import Parameters as ParamFile

import Outputs as Out

import os
import sys


def Quit(self):
    """
    quit MINOTAUR
    
    """
    sys.exit()
    
    
def set_float(entry, value):
    """
    sets the value of an entry in the GUI to a float

    Parameters
    ----------
    entry : TYPE: tkinter entry
        DESCRIPTION: entry in the GUI accepting a float.
    value : TYPE: float
        DESCRIPTION: taken from an input file.
        
    """
    try:
        entry.delete(0, len(entry.get()))
    except:
        pass
    entry.insert(0, float(value[1]))
        
    
def set_pdb(entry, value):
    """
    sets the PDB ID entry value

    Parameters
    ----------
    entry : TYPE: tkinter entry
        DESCRIPTION: entry in the GUI for the PDB ID.
    value : TYPE: float
        DESCRIPTION: taken from an input file.
        
    """
    if len(value[1]) == 4:
        try:
            entry.delete(0, len(entry.get()))
        except:
            pass
        entry.insert(0, value[1])
        

def set_lim(entry, value):
    """
    sets the limits for the initial parameters guesses

    Parameters
    ----------
    entry : TYPE: array of tkinter entry
        DESCRIPTION: entries in the GUI for the limits.
    value : TYPE: array
        DESCRIPTION: contains the parameter name, lower and upper of the given parameter 
        
    """
    for i in range(len(entry[0])):
        set_float(entry[0][i], [0, value[3*i+1]])
        set_float(entry[1][i], [0, value[3*i+2]])
    
    
def set_selection(entry, value):
    """
    sets a scroling menu.

    Parameters
    ----------
    entry : TYPE: tkinter entry
        DESCRIPTION: entry in the GUI for a menu.
    value : TYPE: str
        DESCRIPTION: taken from an input file.
        
    """
    entry.set(value[1])
    
    
def set_button(entry, value):
    """
    sets a checkbox button.

    Parameters
    ----------
    entry : TYPE: tkinter entry
        DESCRIPTION: entry in the GUI for a menu.
    value : TYPE: str
        DESCRIPTION: taken from an input file.
        
    """
    if value[1] == 'True':
        entry.set(1)
    else:
        entry.set(0)
    
    
def new_entry(self, label, width = None):
    """
    creates a new entry accepting text in the GUI

    Parameters
    ----------
    label : TYPE: str
        DESCRIPTION: instructions to write in the entry.
    width : TYPE, optional: float
        DESCRIPTION: width of the entry box. The default is None.

    Returns
    -------
    entry : TYPE: tkinter entry
        DESCRIPTION: new entry in the GUI.

    """
    text_var = tk.IntVar()
    text_var.set(label)
    entry = tk.Entry(self.param, textvariable = text_var, width = width)
    
    return entry
    
    
def load_old_param(self, filename = None):
    """
    loads parameters from an input file to setup the GUI.
    creates a dictionnary containing labels for each entry in the GUI with the 
    corresponding line in the input GUI setup file, the entry in the GUI and
    configuration functions.

    Parameters
    ----------
    filename : TYPE, optional: str
        DESCRIPTION: path to file containing parameters. The default is None.
        
    """
    if filename == None:
        filename = askopenfilename(initialdir = self.browsing_dir)
        
    self.experimental_setup = None
    self.field_calibration = None
    self.int_relax_folder = None
    self.other_input_file = None
    previous_set_up = {'tc': {'entry': self.TauC_entry, 'set function': set_float},
                      'data_scaling': {'entry' : self.Data_scaling_var, 'set function' : set_button},
                     'shuttling_type': {'entry' : self.Shuttling_TYPE, 'set function' : set_selection},
                     'n_walker': {'entry': self.number_walker_entry, 'set function': set_float},
                     'n_mcmc': {'entry': self.number_mcmc_steps_entry, 'set function': set_float},
                     'pdb': {'entry': self.PDB_id, 'set function': set_pdb},
                     'exp_setup': {'text config': self.ExpSetUp_config},
                     'field_cal': {'text config': self.FieldCal_config},
                     'int_relaxometry': {'text config': self.Irelaxom_config},
                     'input_file': {'text config': self.OtherInput_config},
                     'OrderParam': {'entry': [self.MinS2_entry, self.MaxS2_entry], 'set function': set_lim},
                     'CorrTimes': {'entry': [self.MinTau_entry, self.MaxTau_entry], 'set function': set_lim},
                     'others': {'entry': [self.MinOthers_entry, self.MaxOthers_entry], 'set function': set_lim}}
    
    paramFile = open(filename, 'r')
    paramFile.readline()
    
    for entry in previous_set_up.keys():
        Line = paramFile.readline().split('\n')[0].split('\t')
        Line = list(filter(lambda x: x!='', Line))
        previous_set_up[entry]['data'] = Line
    previous_set_up['high_field'] = []
    while True:
        Line = paramFile.readline().split('\n')[0].split('\t')
        Line = list(filter(lambda x: x!='', Line))
        if len(Line) < 2:
            break
        previous_set_up['high_field'].append(Line)
        
    paramFile.close()
        
    set_gui(self, previous_set_up)
    
    
def declare_gui(self):
    """
    creates all the entries for the GUI.
    
    """
    #main buttons
    self.Quit = tk.Button(self.param, text="Quit", command = lambda: Quit(self))
    self.Start = tk.Button(self.param, text="Start", bg="green", command = self.Get_Data)
    self.Load = tk.Button(self.param, text="Load previous parameters", command = lambda: load_old_param(self))
    self.Save = tk.Button(self.param, text="Save current parameters", command = lambda: Out.write_gui_pre_saved_param(self))
    
    #Fitting parameters
    self.TauC_entry = new_entry(self, "Correlation time TauC (sec)")
    
    self.Data_scaling_var = tk.IntVar()
    self.Data_scaling = tk.Checkbutton(self.param, text='Scale data (av=0, std=1)', variable = self.Data_scaling_var, onvalue=1, offvalue=0)

    self.Shuttling_TYPE = tk.StringVar()
    self.Shuttling_TYPE.set("Type of shuttling")
    self.shuttling_type = tk.OptionMenu(self.param, self.Shuttling_TYPE, "Constant Speed", "Constant Acceleration", "Bruker 2024 design")
    
    self.number_mcmc_steps_entry = new_entry(self, "MCMC - Number of steps", width=25)
    self.number_walker_entry = new_entry(self, "MCMC - Number of chains", width=25)
    
    #Limits Part
    self.MinS2_entry = np.empty(shape = len(ParamFile.Names['OrderParam']), dtype = tk.Entry)
    self.MaxS2_entry = np.empty(shape = len(ParamFile.Names['OrderParam']), dtype = tk.Entry)
    for c, order_param in enumerate(ParamFile.Names['OrderParam']):
        self.MinS2_entry[c] = new_entry(self, '0.0', width=10)
        self.MaxS2_entry[c] = new_entry(self, '1.0', width=10)
    
    self.MinTau_entry = np.empty(shape = len(ParamFile.Names['CorrTimes']), dtype = tk.Entry)
    self.MaxTau_entry = np.empty(shape = len(ParamFile.Names['CorrTimes']), dtype = tk.Entry)
    for c, corr_time in enumerate(ParamFile.Names['CorrTimes']):
        self.MinTau_entry[c] = new_entry(self, '', width=10)
        self.MaxTau_entry[c] = new_entry(self, '', width=10)
        
    self.MinOthers_entry = np.empty(shape = len(ParamFile.Names['others']), dtype = tk.Entry)
    self.MaxOthers_entry = np.empty(shape = len(ParamFile.Names['others']), dtype = tk.Entry)
    for c, others in enumerate(ParamFile.Names['others']):
        self.MinOthers_entry[c] = new_entry(self, '', width=10)
        self.MaxOthers_entry[c] = new_entry(self, '', width=10)
        
    #Load files part
    self.FieldCal_config = tk.Button(self.param, text="Field calibration file", command = lambda: browse_field_calibration(self))
    self.ExpSetUp_config = tk.Button(self.param, text="Experimental setup file", command = lambda: browse_exp_setup(self))
    self.Irelaxom_config = tk.Button(self.param, text="Relaxometry Intensity folder", command = lambda: browse_relaxometry_int(self))
    self.OtherInput_config = tk.Button(self.param, text="Other inputs file", command = lambda: browse_other_inputs(self))

    self.RelaxationDataSet = []
    self.Nlabel = 1
    self.RelaxationDataSet.append(tk.Button(self.param, text="High field rates", command = lambda: browse_rates(self, ([0]))))
    
    self.MenuRelaxType = {}
    for i in range(len(ParamFile.RelaxationRates)):
        self.MenuRelaxType[str(ParamFile.RelaxationRates[i])] = i+1
    self.relaxation_data_type = []
    self.RelaxationType = []
    self.relaxation_data_type.append(tk.StringVar())
    self.relaxation_data_type[-1].set("Data type")
    self.RelaxationType.append(tk.OptionMenu(self.param, self.relaxation_data_type[-1], *self.MenuRelaxType.keys()))
    
    self.fields = ['']
    self.fields[0] = new_entry(self, "High field (T)", width=10)
    
    self.AddButton = tk.Button(self.param, text="Add high field rates", command = lambda: add_high_field(self))
    
    #PDB ID
    self.PDB_id = new_entry(self, "", width=10)
    

def set_gui(self, previous_set_up):
    """
    sets the GUI according to previous setup.

    Parameters
    ----------
    previous_set_up : TYPE: dictionnary
        DESCRIPTION: created by the function load_old_param.
        
    """
    no_error = True
    for entry in ['tc', 'data_scaling', 'shuttling_type', 'n_walker', 'n_mcmc', 'OrderParam', 'CorrTimes', 'others']:
        try:
            previous_set_up[entry]['set function'](previous_set_up[entry]['entry'],
                                                   previous_set_up[entry]['data'])
        except:
            no_error = False
    try:
        previous_set_up[entry]['set function'](previous_set_up['pdb']['entry'],
                                               previous_set_up['pdb']['data'])
    except:
        pass
    for entry in ['exp_setup', 'field_cal', 'int_relaxometry', 'input_file']:
        if os.path.exists(previous_set_up[entry]['data'][1]):
            previous_set_up[entry]['text config'].config(text=os.path.basename(previous_set_up[entry]['data'][1]))
        else:
            no_error = False
            
    self.experimental_setup = previous_set_up['exp_setup']['data'][1]
    self.field_calibration = previous_set_up['field_cal']['data'][1]
    self.int_relax_folder = previous_set_up['int_relaxometry']['data'][1]
    self.other_input_file = previous_set_up['input_file']['data'][1]
    
    self.AddButton.grid_forget()
    for i in range(len(self.RelaxationDataSet)):
        self.RelaxationDataSet[i].grid_forget()
        self.RelaxationType[i].grid_forget()
        self.fields[i].grid_forget()  
        
    self.RATES = []
    self.RelaxationDataSet = []
    self.relaxation_data_type = []
    self.RelaxationType = []
    self.fields = []
    for n_high_field in range(len(previous_set_up['high_field'])):
        try:
            if os.path.exists(previous_set_up['high_field'][n_high_field][1]):
                add_high_field(self)
                self.RelaxationDataSet[-1].config(text=os.path.basename(previous_set_up['high_field'][n_high_field][1]))
                self.RATES.append(previous_set_up['high_field'][n_high_field][1])
                self.relaxation_data_type[-1].set(previous_set_up['high_field'][n_high_field][2])
                self.fields[-1].delete(0, "end")
                self.fields[-1].insert(0, float(previous_set_up['high_field'][n_high_field][3]))
                    
                if len(self.RelaxationDataSet) % 2 == 0:
                    self.RelaxationDataSet[-1].grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=4, columnspan=2)
                    self.RelaxationType[-1].grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=6)
                    self.fields[-1].grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=7)
                else:
                    self.RelaxationDataSet[-1].grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=0, columnspan=2)
                    self.RelaxationType[-1].grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=2)
                    self.fields[-1].grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=3)
            else:
                no_error = False
        except:
            no_error = False
        
    if len(self.RelaxationDataSet) % 2 == 0:
        self.AddButton.grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=8, columnspan=2)
    else:
        self.AddButton.grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=4, columnspan=2)
        
    if not no_error:
        self.ErrorText = tk.Label(self.param, text = "Incomplete loading", fg="orange")
        self.ErrorText.grid(row=0, column=3)
        
    #if the user indicated it in the command line, immediately start
    if no_error and len(sys.argv) > 3 and str(sys.argv[3]) == '-start':
        self.Get_Data()
    

def browse_field_calibration(self):
    """
    browsing function for the field calibration file
    
    """
    self.field_calibration = askopenfilename(initialdir = self.browsing_dir)
    self.FieldCal_config.config(text=os.path.basename(self.field_calibration))
    self.browsing_dir = self.field_calibration
    
    
def browse_exp_setup(self):
    """
    browsing function for the experimental setup file
    
    """
    self.experimental_setup = askopenfilename(initialdir = self.browsing_dir)
    self.ExpSetUp_config.config(text=os.path.basename(self.experimental_setup))
    self.browsing_dir = self.experimental_setup
    
    
def browse_relaxometry_int(self):
    """
    browsing function for the relaxometry intensity decay folder
    
    """
    self.int_relax_folder = askdirectory(initialdir = self.browsing_dir)
    self.Irelaxom_config.config(text=os.path.basename(self.int_relax_folder))
    self.browsing_dir = self.int_relax_folder
    
    
def browse_other_inputs(self):
    """
    browsing function for the other input file
    
    """
    self.other_input_file = askopenfilename(initialdir = self.browsing_dir)
    self.OtherInput_config.config(text=os.path.basename(self.other_input_file))
    self.browsing_dir = self.other_input_file
    
    
def browse_rates(self, n):
    """
    browsing function for the high-field rate files

    Parameters
    ----------
    n : TYPE: list of one int
        DESCRIPTION: index of the rate being considered in the array of high-field entries.
        
    """
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
        
        
def add_high_field(self):
    """
    adds entries for another high-field data set.
    
    """
    N = ([self.Nlabel])
    self.Nlabel += 1
    self.AddButton.grid_forget()
    
    self.RelaxationDataSet.append(tk.Button(self.param, text="High field rates", command = lambda: browse_rates(self, N)))

    self.relaxation_data_type.append(tk.StringVar())
    self.relaxation_data_type[-1].set("Data Type")
    self.RelaxationType.append(tk.OptionMenu(self.param, self.relaxation_data_type[-1], *self.MenuRelaxType.keys()))

    self.fields.append('')
    self.fields[-1] = new_entry(self, "High field (T)", width=10)
    
    if len(self.RelaxationDataSet) % 2 == 0:
        self.RelaxationDataSet[-1].grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=4, columnspan=2)
        self.RelaxationType[-1].grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=6)
        self.fields[-1].grid(row=int(self.nline+3+len(self.RelaxationDataSet)/2), column=7)
        
    else:
        self.RelaxationDataSet[-1].grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=0, columnspan=2)
        self.RelaxationType[-1].grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=2)
        self.fields[-1].grid(row=int(self.nline+3+(1+len(self.RelaxationDataSet))/2), column=3)
        
        
def pack(self):
    """
    places the entries in the GUI grid.
    
    """
    self.param.title("MINOTAUR")
    self.param.geometry("1300x800")
    self.grid(row=0, column=0)
    
    self.img = ImageTk.PhotoImage(Image.open("lib/Logo.jpg"))
    panel = tk.Label(self.param, image=self.img)
    panel.grid(row=0, column=2, columnspan=3, rowspan=2)
        
    ParamText = tk.Label(self.param, text="Set fitting parameters:", fg="blue")

    LimitText = tk.Label(self.param, text="Set Limits for MCMC initial guess:", fg="blue")
    MinS2Text = tk.Label(self.param, text="Min")
    MaxS2Text = tk.Label(self.param, text="Max")

    MinTauText = tk.Label(self.param, text="Min")
    MaxTauText = tk.Label(self.param, text="Max")
    
    MinOtherText = tk.Label(self.param, text="Min")
    MaxOtherText = tk.Label(self.param, text="Max")
    
    TextFiles = tk.Label(self.param, text = "Load files:", fg="blue")
    
    TextPDB = tk.Label(self.param, text = "4-letter PDB ID (if available):", fg = "blue")
    
    self.Quit.grid(row=0, column=6)
    self.Start.grid(row=0, column=7)
    self.Load.grid(row=0, column=0, columnspan=3)
    self.Save.grid(row=1, column=6, columnspan=2)
    
    ParamText.grid(row=1, column=0)
    self.TauC_entry.grid(row=2, column=0, columnspan=3)
    self.Data_scaling.grid(row=3, column=0, columnspan=3)
    self.shuttling_type.grid(row=2, column=3, columnspan=2)
    self.number_walker_entry.grid(row=3, column=5, columnspan=3)
    self.number_mcmc_steps_entry.grid(row=2, column=5, columnspan=3)
    
    LimitText.grid(row=5, column=0)
    if len(ParamFile.Names['OrderParam']) > 0:
        MinS2Text.grid(row=6, column=1)
        MaxS2Text.grid(row=6, column=2)
    if len(ParamFile.Names['CorrTimes']) > 0:
        MinTauText.grid(row=6, column=4)
        MaxTauText.grid(row=6, column=5)
    if len(ParamFile.Names['others']) > 0:
        MinOtherText.grid(row=6, column=7)
        MaxOtherText.grid(row=6, column=8)
    
    s2Text = []
    for c, order_param in enumerate(ParamFile.Names['OrderParam']):
        s2Text.append(tk.Label(self.param, text=order_param))
        s2Text[c].grid(row=7+c, column=0)
        self.MinS2_entry[c].grid(row=7+c, column=1)
        self.MaxS2_entry[c].grid(row=7+c, column=2)
        
    tauText = []
    for c, corr_time in enumerate(ParamFile.Names['CorrTimes']):
        tauText.append(tk.Label(self.param, text=corr_time))
        tauText[c].grid(row=7+c, column=3)
        self.MinTau_entry[c].grid(row=7+c, column=4)
        self.MaxTau_entry[c].grid(row=7+c, column=5)
        
    OthersText = []
    for c, others in enumerate(ParamFile.Names['others']):
        OthersText.append(tk.Label(self.param, text=others))
        OthersText[c].grid(row=7+c, column=6)
        self.MinOthers_entry[c].grid(row=7+c, column=7)
        self.MaxOthers_entry[c].grid(row=7+c, column=8)
        
    TextPDB.grid(row=self.nline+1, column=6, columnspan=2)
    self.PDB_id.grid(row=self.nline+2, column=6, columnspan=2)
        
    TextFiles.grid(row=self.nline+1, column=0)
    self.FieldCal_config.grid(row=self.nline+2, column=0, columnspan=2)
    self.ExpSetUp_config.grid(row=self.nline+2, column=2)
    self.Irelaxom_config.grid(row=self.nline+3, column=0, columnspan=2)
    self.OtherInput_config.grid(row=self.nline+3, column=2)
    for i in range(len(self.RelaxationDataSet)):
        self.RelaxationDataSet[i].grid(row=self.nline+4+i, column=0, columnspan=2)
        self.RelaxationType[i].grid(row=self.nline+4+i, column=2)
        self.fields[i].grid(row=self.nline+4+i, column=3)