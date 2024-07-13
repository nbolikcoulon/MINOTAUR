#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#########################################################################
#                                                                       #
#                       Creates output text files                       #
#                                                                       #
#########################################################################
import numpy as np

import Parameters as ParamFile

import os

            
def write_gui_pre_saved_param(self):
    """
    writes the status of the GUI once the 'Save current parameters' button is pressed. This file can be used to populate the GUI for a later run.
    
    """
    n = 1
    while True:
        pre_saved_filename = f'PreSaved_Parameters_{n}.txt'
        if not os.path.exists(pre_saved_filename):
            break
        n += 1
    pre_saved_file = open(pre_saved_filename, 'w')
        
    experimental_setup_path = "Not provided"
    field_calibration_path = "Not provided"
    int_relax_path = "Not provided"
    other_input_path = "Not provided"
    if hasattr(self, "experimental_setup"):
        experimental_setup_path = self.experimental_setup
    if hasattr(self, "field_calibration"):
        field_calibration_path = self.field_calibration
    if hasattr(self, "int_relax_folder"):
        int_relax_path = self.int_relax_folder
    if hasattr(self, "other_input_file"):
        other_input_path = self.other_input_file
        
    Nchains = int(float(self.number_walker_entry.get()))
    Nstep = int(float(self.number_mcmc_steps_entry.get()))
    TauC = float(self.TauC_entry.get())
    if self.Data_scaling_var.get() == 1:
        scale_data = 'True'
    else:
        scale_data = 'False'
                
    pre_saved_file.write("PRE SAVED\n")
        
    pre_saved_file.write(f"Global tumbling correlation time\t{TauC}\n")
    pre_saved_file.write(f"Scale Data\t{scale_data}\n")
    pre_saved_file.write(f"Shuttling type\t{self.Shuttling_TYPE.get()}\n")
    pre_saved_file.write(f"Number of chains\t{Nchains}\n")
    pre_saved_file.write(f"Number of MCMC steps\t{Nstep}\n")
    pre_saved_file.write(f"PDB ID\t{self.PDB_id.get()}\n")
    pre_saved_file.write(f"Setup path\t{experimental_setup_path}\n")
    pre_saved_file.write(f"Field calibration path\t{field_calibration_path}\n")
    pre_saved_file.write(f"Relaxometry data path\t{int_relax_path}\n")
    pre_saved_file.write(f"Ohter Input path\t{other_input_path}\n")
        
    for i in range(len(ParamFile.Names['OrderParam'])):
        mins2 = self.MinS2_entry[i].get()
        maxs2 = self.MaxS2_entry[i].get()
        pre_saved_file.write(f"{ParamFile.Names['OrderParam'][i]}\t{mins2}\t{maxs2}\t")
    
    pre_saved_file.write("\n")
    for i in range(len(ParamFile.Names['CorrTimes'])):
        mintau = self.MinTau_entry[i].get()
        maxtau = self.MaxTau_entry[i].get()
        pre_saved_file.write(f"{ParamFile.Names['CorrTimes'][i]}\t{mintau}\t{maxtau}\t")
        
    pre_saved_file.write("\n")
    for i in range(len(ParamFile.Names['others'])):
        minOthers = self.MinOthers_entry[i].get()
        maxOthers = self.MaxOthers_entry[i].get()
        pre_saved_file.write(f"{ParamFile.Names['others'][i]}\t{minOthers}\t{maxOthers}\t")
            
    for DataSet in range(len(self.RATES)):
        pre_saved_file.write("\n")
        pre_saved_file.write(f"HF data set {DataSet+1}\t{self.RATES[DataSet]}\t{self.relaxation_data_type[DataSet].get()}\t{self.fields[DataSet].get()}")

    pre_saved_file.close()
        
        
def write_gui_param(self):
    """
    writes the status of the GUI once calculations are about to start. This file can be used to populate the GUI for a later run.
    
    """
    saved_file = open(f'{self.directory_name}/InputFiles/Parameters.txt', 'w')
        
    saved_file.write("UNCORRECT\n")
    saved_file.write(f"Global tumbling correlation time\t{self.TauC_entry.get()}\n")
    if self.Data_scaling_var.get() == 1:
        scale_data = 'True'
    else:
        scale_data = 'False'
    saved_file.write(f"Scale Data\t{scale_data}\n")
    saved_file.write(f"Shuttling type\t{self.shuttling_type}\n")
    saved_file.write(f"Number of Chains\t{self.number_walker_entry.get()}\n")
    saved_file.write(f"Number of MCMC steps\t{self.number_mcmc_steps_entry.get()}\n")
    saved_file.write("PDB ID\t")
    if self.check_PDB:
        saved_file.write(self.PDB)
    else:
        saved_file.write('Not available')
    saved_file.write("\n")
    saved_file.write(f"Setup Path\t{self.experimental_setup}\n")
    saved_file.write(f"Field calibration path\t{self.field_calibration}\n")
    saved_file.write(f"Relaxometry data path\t{self.int_relax_folder}\n")
    saved_file.write(f"Other inputs path\t{self.other_input_file}")
        
    for name in ParamFile.Names.keys():
        saved_file.write("\n")
        for n, limits in enumerate(self.bonds_starting_point[name]):
            saved_file.write(f"{ParamFile.Names[name][n]}\t{limits[0]}\t{limits[1]}\t")
            
    for DataSet in range(len(self.RelaxationDataSet)):
        saved_file.write("\n")
        saved_file.write(f"HF data set {DataSet+1}\t{self.RATES[DataSet]}\t{self.relaxation_data_type[DataSet].get()}\t{self.fields[DataSet].get()}")
        
    saved_file = open(f'{self.directory_name}/InputFiles/Parameters.txt', 'w')
    saved_file.write("CONFIRMED\n")
    saved_file.close()
    
        
def write_trajectory(self, trajectory, AA):
    """
    Writes the MCMC trajectory.

    Parameters
    ----------
    trajectory : TYPE: array
        DESCRIPTION: contains the evolution of parameters for each chains at each steps.
    AA : TYPE: str
        DESCRIPTION: residue.
        
    """
    traj_file = open(f'{self.dir_fit_output}/Trajectories/Trajectory_Residue_{AA}.txt', 'w')
    
    traj_file.write('Step')
    for param in self.list_parameters:
        traj_file.write(f'\t{param}')
    traj_file.write('\tf')
        
    for count, step in enumerate(trajectory):
        traj_file.write(f'\n{count}')
        for param in step:
            traj_file.write(f'\t{param}')
    traj_file.close()
    
    
def write_MCMC_parameters(self):
    """
    writes the output file containing the MCMC parameters (mean, positive and negative errors) along with the mean acceptance fraction.
    
    """
    out_file = open(f'{self.dir_fit_output}/MCMCparameters.txt', 'w')
    
    out_file.write("Residue\tMAF")
    for param in self.list_parameters:
        out_file.write(f"\t{param}\t+ error\t- error")
    out_file.write("\tlnf\t+ error\t- error")
    
    for AA in self.residue_list:
        out_file.write(f"\n{AA}")
        for param in range(len(self.list_parameters)+1):
            out_file.write(f'\t{self.MCMC_param[AA]["Mean"][param]}\t{self.MCMC_param[AA]["+"][param]}\t{self.MCMC_param[AA]["-"][param]}')
        
    out_file.close()
    
    
def check_PDB_relevant(residue_list):
    """
    checks whether a PDB file is relevant, i.e. more than one residue is analyzed.

    Parameters
    ----------
    residue_list : TYPE: array
        DESCRIPTION: contains the list of residues.

    Returns
    -------
    bool
        DESCRIPTION: True if a file can be written.

    """
    if len(residue_list) == 1:
        print('Only one residue provided: PDB coloring files will not be created')
        return False
    return True


def normalize_mcmc_parameters(mcmc_output, chi2_list, labels):
    """
    normalizes the MCMC parameters and chi2 from 0 to 1 to define colors in the PDB file.

    Parameters
    ----------
    mcmc_output : TYPE: dictionnary
        DESCRIPTION: contains all the mcmc parameters.
    chi2_list : TYPE: dictionnary
        DESCRIPTION: contains the chi2 values.
    labels : TYPE: array
        DESCRIPTION: name of the parameters.

    Returns
    -------
    normalized_param : TYPE: dictionnary
        DESCRIPTION: normalized MCMC parameters and chi2.

    """
    maxParam, minParam = {}, {}
    all_param = {}
    
    for count, param in enumerate(labels[:-1]):
        all_param[param] = {}
        for AA in mcmc_output.keys():
            all_param[param][AA] = mcmc_output[AA]['Mean'][count]
        minParam[param], maxParam[param] = min(all_param[param].values()), max(all_param[param].values())
        
    all_param['Chi2'] = chi2_list
    minParam['Chi2'], maxParam['Chi2'] = min(chi2_list.values()), max(chi2_list.values())
    
    normalized_param = {}
    for param in labels:
        normalized_param[param] = {}
        for AA in mcmc_output.keys():
            normalized_param[param][AA] = (all_param[param][AA] - minParam[param]) / (maxParam[param] - minParam[param])
            
    return normalized_param
        
        
def write_PDB(self, chi2_list):
    """
    writes the PDB file, if relevant.

    Parameters
    ----------
    chi2_list : TYPE: dictionnary
        DESCRIPTION: list of chi2.
        
    """
    if not check_PDB_relevant(self.residue_list):
        return
        
    os.makedirs(f'{self.dir_fit_output}/PDBFiles')
    
    labels = np.copy(self.list_parameters)
    labels = np.append(labels, 'Chi2')
    normalized_param = normalize_mcmc_parameters(self.MCMC_param, chi2_list, labels)
                    
    for param in labels:
        PDB_file = open(f'{self.dir_fit_output}/PDBFiles/{self.PDB}_{param}.pml', 'w')
                
        PDB_file.write(f"load {self.PDB}.pdb\n")
        PDB_file.write("hide\n")
        PDB_file.write("show cartoon, all\n\n")
                
        PDB_file.write("#Color the whole structure in grey as default color\n")
        PDB_file.write(f"color grey, {self.PDB}\n\n")
                
        PDB_file.write("#Define a color in the blue to red range for each residue we have data for\n")
        for AA in self.residue_list:
            PDB_file.write(f"set_color resi_{AA} = [{normalized_param[param][AA]}, 0.0, {1.0 - normalized_param[param][AA]}]\n")

        PDB_file.write("\n#Coloring each residue according to its color\n")
        for AA in self.residue_list:
            PDB_file.write("color resi_{AA}, resi {AA}\n")
                    
        PDB_file.close()
        
        
def write_LF_R1(self):
    """
    writes the back-calculated and fitted from an exponential decay low field R1. Scaling factors to reproduce the decays are also written here.
    
    """
    frelax = open(f'{self.dir_output_res}/LowField_R1.txt', 'w')
    
    frelax.write("\t\tField\tIntensity scaling\tBack-calculated R1\tFitted R1\n")
            
    for AA in self.residue_list:
        frelax.write(f"Residue {AA}\n")
        for exp in self.set_up.keys():
            frelax.write(f"{self.B0_low_field[exp]}\t{self.scaling_factors_intensities[AA][exp]}\t{self.back_calculated_R1_LF[AA][exp]}\t{self.fitted_R1_LF[AA][exp]}\n")
            
    frelax.close()