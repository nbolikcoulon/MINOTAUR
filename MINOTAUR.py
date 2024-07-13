#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tkinter as tk
import numpy as np
import time
import sys
import os

from matplotlib.backends.backend_pdf import PdfPages

CSFolder = str(sys.argv[1])
RatesFolder = CSFolder + "/Rates"

sys.path.append(os.path.abspath(CSFolder))
sys.path.append(os.path.abspath(RatesFolder))
import Parameters as ParamFile

sys.path.append(os.path.abspath('lib'))
import FitFunctions as FitF
import ShuttlingSimulation as ShSim
import MonteCarlo as MCMC
import Outputs as Out
import FigOutputs as FigOut
import ReadInput
import GUI

os.environ["OMP_NUM_THREADS"] = "1"


class MINOTAUR(tk.Frame):
    def __init__(self, parent):
        
        tk.Frame.__init__(self, parent)
        
        self.parent = parent
        
        self.relaxation_function = ParamFile.ImportFunc()
        
        self.browsing_dir = os.path.dirname(os.path.abspath(__file__))
        
        self.RATES = []
        self.nline = 6 + max(len(ParamFile.Names['OrderParam']), len(ParamFile.Names['CorrTimes']), len(ParamFile.Names['others']))
        
        self.GUI()
        

############ GUI coding - main window
    def GUI(self):
        """
        Deals with the GUI coding:
            - get info from paramters.txt file
            - setup the GUI
            
        """
        self.param = tk.Toplevel()
   
        #declare GUI
        GUI.declare_gui(self)
        
        #setting up in the GUI
        GUI.pack(self)
        
        #option to indicate a load-file in the command line
        if len(sys.argv) > 2:
            load_file = str(sys.argv[2])
            GUI.load_old_param(self, load_file)
        

############ GUI coding - read data
    def Get_Data(self):
        """
        Read the data as provided in the GUI
        Data (relaxation rates and intensities) are scaled here as well. See function for info
        Field profile is created. See function for info
        
        """
        self.begin = time.strftime("%H") + "H" + time.strftime("%M")
        
        ReadInput.copy_input(self)
        
        #Read the data
        ReadInput.read_gui_parameters(self)

        #Field calibration
        field_cal, self.B0_cal_coeff, self.Static_MagField, self.tunnel_position, self.tunnel_field = ReadInput.read_field_calibration(self.field_calibration)
        
        #experimental set up
        self.set_up = ReadInput.read_exp_setup(self.experimental_setup, self.shuttling_type, field_cal)
        ReadInput.copy_relaxometry_decays(self)
        
        #high field Files
        self.HF_data = ReadInput.read_high_field_data(self.RATES, self.relaxation_data_type, self.fields)
        
        #LF intensities
        self.intensities = ReadInput.read_relaxometry_decays(self.int_relax_folder, self.set_up)
        
        #Other inputs
        self.other_inputs = ReadInput.read_other_input(self.other_input_file)
        
        #residue list
        ReadInput.get_residue_list(self)
        
        #scale data for MCMC
        ReadInput.scale_data(self)
        
        #Plot the field profile with the considered fields during relaxometry
        self.B0_low_field = FitF.Get_B0_low_field(self.set_up, self.B0_cal_coeff, self.tunnel_position, self.tunnel_field)
        FigOut.plot_field_profile(field_cal, self.set_up, self.B0_cal_coeff, self.B0_low_field, self.tunnel_position, self.tunnel_field, self.dir_fit_output)
        
        #Write the Parameters file
        Out.write_gui_param(self)
        
        self.Initialize()


############ data analysis
    def Initialize(self):
        """
        initialization prior to the run:
            - optimize the shuttling distance increment
            - create the field lists used during the shuttling simulations
            - optimize the propagator calculation methods
            
        """
        ## Preparation before MCMC
        # trajectories
        self.Increment = ShSim.optimize_shuttling_increment(self, ParamFile.PositionAuto)
        self.shuttling_fields, self.shuttling_delays = ShSim.make_field_list(self, self.Increment)
        
        #propagator calculations
        ShSim.optimize_calc_propagator(self)
        
        #variable declation
        self.MCMC_param = {}
        self.Acceptance = {}
        self.final_simulated_intensities = {}
        self.back_calculated_R1_LF = {}
        self.fitted_R1_LF = {}
        self.scaling_factors_intensities = {}

        #create pdf for figures  
        os.makedirs(f'{self.dir_fit_output}/Correlations')
        os.makedirs(f'{self.dir_fit_output}/Trajectories')
        os.makedirs(f'{self.dir_output_res}/Intensities')
        
        self.pdf_trajectories = PdfPages(f'{self.dir_fit_output}/Trajectories/All_trajectories.pdf' )
        self.pdf_correlations = PdfPages(f'{self.dir_fit_output}/Correlations/All_correlations.pdf' )
        self.pdf_rates = {RelaxRate: PdfPages(f'{self.dir_output_res}/{RelaxRate}.pdf') for RelaxRate in ParamFile.RelaxationRates}
        
        self.Calculations()
        
        
    def Calculations(self):
        """
        Perform the MCMC calculations
        makes residue-specific figures based on MCMC results
        
        """
        print("\nMonte Carlo")
        
        for AA in self.residue_list:
            print(f"\n Residue {AA}\n")
            
            self.MCMC_param[AA], self.Acceptance[AA], Full_Trajectory = MCMC.Markov_Chain_Monte_Carlo(self, AA)
            Out.write_trajectory(self, Full_Trajectory, AA)
            self.final_simulated_intensities[AA] = ShSim.Expected_Values(np.asarray(self.MCMC_param[AA]['Mean'][:-1]), self.TauC, self.other_inputs[AA], self.Static_MagField,
                                                                    self.shuttling_fields, self.shuttling_delays, self.set_up,
                                                                    self.B0_low_field, ParamFile.PositionAuto, self.PropFunction)
            print("\n Making figures")
            
            self.scaling_factors_intensities[AA] = MCMC.scaling_factor(self.final_simulated_intensities[AA], self.intensities, AA)
            self.back_calculated_R1_LF[AA] = FitF.Calc_R1_LF(self.MCMC_param[AA]['Mean'][:-1], self.TauC, self.other_inputs[AA], self.B0_low_field)
            self.fitted_R1_LF[AA] = FitF.Fit_R1_LF(self.intensities, AA)
        
            #Plot the intensities
            FigOut.plot_intensities(self, AA)
                        
            #Plot the relaxation rates 
            for RelaxRate in ParamFile.RelaxationRates:
                if RelaxRate == 'R1':
                    FigOut.plot_R1(self, AA)
                else:
                    FigOut.plot_rate(self, RelaxRate, AA)
                    
        #close all the pdf figures
        for RelaxRate in ParamFile.RelaxationRates:
            self.pdf_rates[RelaxRate].close()
        self.pdf_trajectories.close()
        self.pdf_correlations.close()
            
        self.WriteResult()
        

    def WriteResult(self):
        """
        Creates outputs (rest of the figures and text files)
        
        """
        print("\nWriting final results...")
        
        os.makedirs(f'{self.directory_name}/PlotParameters')
        
        Out.write_MCMC_parameters(self)        #file containing the parameters of the spectral density function extracted from the MCMC

        #Draw the Chi2
        All_Chi2 = {}
        for AA in self.residue_list:
            All_Chi2[AA] = FitF.Chi2_TOT(self.MCMC_param[AA]['Mean'][:-1], self.final_simulated_intensities, self.scaling_factors_intensities,
                                         self.intensities, self.HF_data, self.TauC, self.other_inputs, AA)
        FigOut.plot_chi2(f'{self.directory_name}/PlotParameters', All_Chi2)

        #Draw the spectral density function parameters
        for count, param in enumerate(self.list_parameters):
            FigOut.plot_dynamic_parameters(f'{self.directory_name}/PlotParameters', self.MCMC_param, param, count)

        #Write the LF R1
        Out.write_LF_R1(self)        #file containing the scaling factors for intensities, back-calculated and fitted low field R1

        #Write the PDB files
        if self.check_PDB:
            Out.write_PDB(self, All_Chi2)

        print("    Writing results: Done\n")
        
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
    MINOTAUR(root)
    root.mainloop()