#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##########################################################
#                                                        #
#                       Read files                       #
#                                                        #
##########################################################
from shutil import copyfile
from distutils.dir_util import copy_tree

import Parameters as ParamFile
import FitFunctions as FitF
import MonteCarlo as MCMC

import numpy as np
import time
import sys
import os


def copy_input(self):
    """
    copies all the input files into the InputFiles directory
    relaxometry decay files are not copied at this point as they need to be checked
    
    """
    #create main output directory
    result_directory = "Results"

    if not os.path.exists(result_directory):
        os.makedirs(result_directory)
        
    self.directory_name = f'{result_directory}/{time.strftime("%Y-%m-%d")}'
    n = 2
    while True:
        if not os.path.exists(self.directory_name):
            os.makedirs(self.directory_name)
            break
        else:
            self.directory_name = f'{result_directory}/{time.strftime("%Y-%m-%d")}_{n}'
            n += 1

    #create subdirectories
    dir_input = f'{self.directory_name}/InputFiles'
    os.makedirs(dir_input)
            
    self.dir_output_res = f'{self.directory_name}/FitAllResidues'
    os.makedirs(self.dir_output_res)
    
    self.dir_fit_output = f'{self.directory_name}/FittingResults'
    os.makedirs(self.dir_fit_output)
    
    #Copy the directory containing the C-script functions
    expression_folder = f'{dir_input}/Expressions_and_constraints'
    os.makedirs(expression_folder)
    copy_tree(f'{sys.argv[1]}', expression_folder)
    
    #Copy the input files
    os.makedirs(f'{dir_input}/HFRelaxationRates')
    for rate_name in ParamFile.RelaxationRates:
        os.makedirs(f'{dir_input}/HFRelaxationRates/{rate_name}')
        
    for c, rate_file in enumerate(self.RATES):
        copyfile(rate_file, f'{dir_input}/HFRelaxationRates/{self.relaxation_data_type[c].get()}/{os.path.basename(rate_file)}')

    file_names = {'ExpSetUp.txt': self.experimental_setup,
                  'FieldCalibration.txt': self.field_calibration,
                  'OtherInputs.txt': self.other_input_file}
    for f in file_names.keys():
        copyfile(file_names[f], f'{dir_input}/{f}')
        
        
def copy_relaxometry_decays(self):
    """
    copies intensity relaxometry decay files
    makes sure that each experiment defined in the setup file has an associated decay file
    
    """
    os.makedirs(f'{self.directory_name}/InputFiles/RelaxometryIntensities')
    
    valid_file_ext = ['txt', 'dat', 'csv']
    
    for exp in self.set_up.keys():
        file_found = False
        for extension in valid_file_ext:
            interensity_file = f'{self.int_relax_folder}/{exp}.{extension}'
            if os.path.isfile(interensity_file):
                copyfile(interensity_file, f'{self.directory_name}/InputFiles/RelaxometryIntensities/{exp}.{extension}')
                file_found = True
                break
            
        if not file_found:
            print()
            print(f"Missing the eperiment number {exp} file")
            print()
            sys.exit()
            
            
def read_gui_parameters(self):
    """
    reads all the gui entry parameters
    """
    self.number_mcmc_steps = int(float(self.number_mcmc_steps_entry.get()))
    number_walker_input = int(float(self.number_walker_entry.get()))
    self.number_walker = MCMC.number_walker_check(number_walker_input,
                                                  len(ParamFile.Names['OrderParam']) + len(ParamFile.Names['CorrTimes'])+len(ParamFile.Names['others'])+1)

    self.TauC = float(self.TauC_entry.get())
    if self.Data_scaling_var.get() == 1:
        self.perform_scaling  = 'True'
    else:
        self.perform_scaling  = 'False'
        
    self.list_parameters = []
    for S2 in ParamFile.Names['OrderParam']:
        self.list_parameters.append(S2)
    for tau in ParamFile.Names['CorrTimes']:
        self.list_parameters.append(tau)
    for other in ParamFile.Names['others']:
        self.list_parameters.append(other)
        
    self.bonds_starting_point = {'OrderParam': [], 'CorrTimes': [], 'others': []}
    for i in range(len(ParamFile.Names['OrderParam'])):
        mins2 = float(self.MinS2_entry[i].get())
        maxs2 = float(self.MaxS2_entry[i].get())
        self.bonds_starting_point['OrderParam'].append([mins2, maxs2])
    for i in range(len(ParamFile.Names['CorrTimes'])):
        mintau = float(self.MinTau_entry[i].get())
        maxtau = float(self.MaxTau_entry[i].get())
        self.bonds_starting_point['CorrTimes'].append([mintau, maxtau])
    for i in range(len(ParamFile.Names['others'])):
        minOthers = float(self.MinOthers_entry[i].get())
        maxOthers = float(self.MaxOthers_entry[i].get())
        self.bonds_starting_point['others'].append([minOthers, maxOthers])
    self.bonds_starting_point['lnf'] = [[-10.0, 1.0]]
    
    self.shuttling_type = str(self.Shuttling_TYPE.get())
    
    self.PDB = self.PDB_id.get()
    if len(self.PDB) == 4:
        self.check_PDB = True
    else:
        self.check_PDB = False
    

def read_line(line):
    """
    custom readline function to account for different possible column separations
    possible separations are tabs (\t), columns (,) and spaces (' ')
    we prefer to avoid the use of the pandas library which assumes homogeneous formating
    
    """
    line = line.split('\n')[0]
    if '\t' in line:
        line = line.split('\t')
    elif ',' in line:
        line = line.split(',')
    elif ' ' in line:
        line = line.split(' ')
        
    line = list(filter(lambda x: x!='', line))
    
    return line

            
def read_field_calibration(filename):
    """
    read the field calibration file

    Parameters
    ----------
    filename : TYPE: str
        DESCRIPTION: path to the field calibration file. Two-column file containing the height and field
        
    Returns
    -------
    field_cal : TYPE: dictionnary
        DESCRIPTION: associates the height with the field.
    B0_cal_coeff : TYPE: dictionnary
        DESCRIPTION: contains the polynomial decomposition to compute the field along the sample trajectory.
    B0_HF : TYPE: float
        DESCRIPTION: static magnetic field at which relaxometry experiments are performed.
    tunnel_position: TYPE: float
        DESCRIPTION: position at which the tunnel starts.
    tunnel_field : TYPE: float
        DESCRIPTION: field of the tunnel.
        
    """
    field_cal = {}
    
    f = open(filename, 'r')
    while True:
        line = read_line(f.readline())
        
        try:
            field_cal[float(line[0])] = float(line[1])
        except:
            f.close()
            break
    
    B0_cal_coeff, tunnel_position, tunnel_field = FitF.calibrate_B0(field_cal)
    B0_HF = FitF.Calc_B0(0.0, B0_cal_coeff, tunnel_position, tunnel_field)
        
    return field_cal, B0_cal_coeff, B0_HF, tunnel_position, tunnel_field
        
        
def read_exp_setup(filename, shuttling_type, field_cal):
    """
    reads the setup file

    Parameters
    ----------
    filename : TYPE: str
        DESCRIPTION: path to the setup file.
        multiple-column file containing:
            1) the experiment numer
            2) the shuttling height (or shuttling field in the case of FCC)
            3) the waiting time at low field after shuttling up
            4) the waiting time at high field after shuttling down
            5) the waiting time at high field before shuttling up
            6) the shuttling time to the low field position
            7) the waiting time at low field before shuttling down
            8) the shuttling time to the high field position
            9) all the relaxation delays
    All delays have to be given in milliseconds.
    shuttling_type : TYPE: str
        DESCRIPTION: type of shuttling device
    field_cal : TYPE: dictionnary
        DESCRIPTION: matches distance with field. Required if shuttling_type is 'Bruker 2024 design'.

    Returns
    -------
    set_up : TYPE: dictionnary
        DESCRIPTION: contains all the parameters for each experiments.

    """
    set_up = {}
    
    f = open(filename, 'r')
    while True:
        line = read_line(f.readline())
        
        try:
            exp = int(float(line[0]))
            set_up[exp] = {}
            
            if shuttling_type == 'Bruker 2024 design':
                if line[1] == 'FCC':
                    set_up[exp]['field'] = float(line[2])
                    set_up[exp]['height'] = max(list(field_cal.keys()))
                else:
                    set_up[exp]['height'] = float(line[2])
                set_up[exp]['d22'] = float(line[3]) * 1e-3
                set_up[exp]['d25'] = float(line[4]) * 1e-3
                set_up[exp]['WTHF'] = float(line[5]) * 1e-3
                set_up[exp]['SLF'] = float(line[6]) * 1e-3
                set_up[exp]['WTLF'] = float(line[7]) * 1e-3
                set_up[exp]['SHF'] = float(line[8]) * 1e-3
                
                set_up[exp]['LF_times'] = {}
                set_up[exp]['vc'] = []
                for j in range(len(line)-9):
                    delay = float(line[j+9]) * 1e-3
                    set_up[exp]['vc'].append(delay)
                    set_up[exp]['LF_times'][delay] = set_up[exp]['d22'] + set_up[exp]['WTLF'] + delay
            else:
                set_up[exp]['height'] = float(line[1])
                set_up[exp]['d22'] = float(line[2]) * 1e-3
                set_up[exp]['d25'] = float(line[3]) * 1e-3
                set_up[exp]['WTHF'] = float(line[4]) * 1e-3
                set_up[exp]['SLF'] = float(line[5]) * 1e-3
                set_up[exp]['WTLF'] = float(line[6]) * 1e-3
                set_up[exp]['SHF'] = float(line[7]) * 1e-3
                
                set_up[exp]['LF_times'] = {}
                set_up[exp]['vc'] = []
                for j in range(len(line)-8):
                    delay = float(line[j+8]) * 1e-3
                    set_up[exp]['vc'].append(delay)
                    set_up[exp]['LF_times'][delay] = set_up[exp]['d22'] + set_up[exp]['WTLF'] + delay
            
        except:
            f.close()
            return set_up
        
        
def read_high_field_file(filename):
    """
    read high-field relaxation rate file

    Parameters
    ----------
    filename : TYPE: str
        DESCRIPTION: path to a file containing relaxation rates. A 3-column file with the residue, rate and error.

    Returns
    -------
    hf_data : TYPE: str
        DESCRIPTION: contains the values and rates of one particular rate.

    """
    hf_data = {}
    
    f = open(filename, 'r')
    while True:
        line = read_line(f.readline())
        
        if len(line) > 0:
            aa = line[0]
            
            try:
                hf_data[aa] = {'data': float(line[1]), 'error': float(line[2])}
            except:
                pass
        else:
            f.close()
            return hf_data
        
        
def read_high_field_data(rates_file, data_type, fields):
    """
    reads all the high field rate files
    uses the read_high_field_file function to read each file
    calculate the rate averages in the case of repeats

    Parameters
    ----------
    rates_file : TYPE: array of files
        DESCRIPTION: contains the list of files.
    data_type : TYPE: array of entries
        DESCRIPTION: contains the list of types of rates.
    fields : TYPE: array of entries
        DESCRIPTION: contains the list of fields.

    Returns
    -------
    HF_data : TYPE: dictionnary
        DESCRIPTION: contains all the high field relaxation rates.

    """
    HF_data_all = {relax_type: {} for relax_type in ParamFile.RelaxationRates}
    for dataset in range(len(data_type)):
        relax_type = str(data_type[dataset].get())
        b0 = float(fields[dataset].get())
            
        if b0 not in HF_data_all[relax_type].keys():
            HF_data_all[relax_type][b0] = []
                
        data = read_high_field_file(rates_file[dataset])
        HF_data_all[relax_type][b0].append(data)
        
    #create the average HF data
    HF_data = {}
    for relax_type in HF_data_all.keys():
        HF_data[relax_type] = {}
        
        for field in HF_data_all[relax_type].keys():
            HF_data[relax_type][field] = {aa: {'data': 0.,
                                               'error': 0.} for aa in HF_data_all[relax_type][field][0].keys()}
            
            for aa in HF_data_all[relax_type][field][0].keys():
                n_repeat = len(HF_data_all[relax_type][field])
                for repeat in range(n_repeat):
                    HF_data[relax_type][field][aa]['data'] += HF_data_all[relax_type][field][repeat][aa]['data'] / n_repeat
                    HF_data[relax_type][field][aa]['error'] += HF_data_all[relax_type][field][repeat][aa]['error']**2
                HF_data[relax_type][field][aa]['error'] = np.sqrt(HF_data[relax_type][field][aa]['error']) / n_repeat
                
    return HF_data
        
        
def read_relaxometry_decay_file(filename, delays, exp):
    """
    reads a relaxometry decay file.
    exits if the number of intensities does not match the number of delay times

    Parameters
    ----------
    filename : TYPE: str
        DESCRIPTION: path to the decay file. A multiple-column file containing the residue and
        then pairs of {intensity error} for each delays
    delays : TYPE: array
        DESCRIPTION: contains the list of delays associated to the file.
    exp : TYPE: str
        DESCRIPTION: experiment number.

    Returns
    -------
    data_av : TYPE: dictionnary
        DESCRIPTION: intensities and error, averaged in the case of repeat for one or more relaxation delays.

    """
    #read
    data = {}
    f = open(filename, 'r')
    while True:
        line = read_line(f.readline())
        
        if len(line) > 0:
            aa = line[0]
            data[aa] = {}
            
            if len(line) != int(2*len(delays) + 1):
                    print("")
                    print(f"Experiment number {exp} does not have the correct number of intensities")
                    print("")
                    sys.exit()
                    
            for c, vc in enumerate(delays):
                try:
                    if vc in data[aa].keys() and data[aa][vc] != 'NA':
                        data[aa][vc]['data'].append(float(line[1 + 2*c]))
                        data[aa][vc]['error'].append(float(line[2 + 2*c])**2)
                    else:
                        data[aa][vc] = {'data': [float(line[1 + 2*c])], 'error': [float(line[2 + 2*c])**2]}
                except:
                    if vc not in data[aa].keys():
                        data[aa][vc] = {'data': 'NA', 'error': 'NA'}
        else:
            f.close()
            break
        
    #calculate average if some delays are repeated
    data_av = {aa: {} for aa in data.keys()}
    for aa in data.keys():
        for vc in data[aa].keys():
            if data[aa][vc]['data'] != 'NA':
                data_av[aa][vc] = {'data': np.average(data[aa][vc]['data']),
                                   'error': np.sqrt(sum(data[aa][vc]['error'])) / len(data[aa][vc]['error'])}
            else:
                data_av[aa][vc] = {'data': 'NA', 'error': 'NA'}
        
    return data_av


def read_relaxometry_decays(folder, set_up):
    """
    master function to read relaxometry decay files

    Parameters
    ----------
    folder : TYPE: str
        DESCRIPTION: path to the folder containing the decay files.
    set_up : TYPE: dictionnary
        DESCRIPTION: contains the relaxometry experimental parameters .

    Returns
    -------
    Intensities : TYPE: dictionnary
        DESCRIPTION: intensities and error for all the relaxometry experiments.

    """
    valid_file_ext = ['txt', 'dat', 'csv']
    
    Intensities = {}
    for exp in set_up.keys():
        for extension in valid_file_ext:
            interensity_file = f'{folder}/{exp}.{extension}'
            if os.path.isfile(interensity_file):
                break
            
        Intensities[exp] = read_relaxometry_decay_file(interensity_file, set_up[exp]['vc'], exp)
            
    return Intensities


def scale_intensities(data, perform_scaling):
    """
    scaling of the intensities in preparation for the MCMC. Intensities are scaled such that:
        - the average of each decay equals 0
        - the standard deviation of each decay equals 1.
    Intensities are then transformed as: I' = (I - average(I)) / std(I).

    Parameters
    ----------
    data : TYPE: dictionnary
        DESCRIPTION: contains the intensities and error for one experiment
    perform_scaling : TYPE: Bool
        DESCRIPTION: true if scaling has to be performed.

    Returns
    -------
    scaled_data : TYPE: dictionnary
        DESCRIPTION: scaled intensities and errors.
    average : TYPE: float
        DESCRIPTION: average of the recorded intensities. Equals 0 when perform_scaling=False.
    std : TYPE: float
        DESCRIPTION: standard deviation of the recorded intensities. Equals 1 when perform_scaling=False.

    """
    if perform_scaling == 'True':
        average, av_2 = 0., 0.
        count = 0
        for vc in data.keys():
            if data[vc]['data'] != 'NA':
                count += 1
                average += data[vc]['data']
                av_2 += data[vc]['data']**2
        average /= count
        std = np.sqrt(av_2/count - average**2)
    else:
        average = 0.
        std = 1.
    
    scaled_data = {}
    for vc in data.keys():
        if data[vc]['data'] != 'NA':
            scaled_data[vc] = {'data': (data[vc]['data'] - average) / std,
                                'error': data[vc]['error'] / std}
        else:
            scaled_data[vc] = {'data': 'NA', 'error': 'NA'}
        
    return scaled_data, average, std


def scale_rates(data, aa, perform_scaling):
    """
    scaling of the rates in preparation for the MCMC. Rates are scaled such that:
        - the average of one particular rate over multiple fields equals 0
        - the standard deviation of one particular rate over multiple fields equals 1.
    Rates are then transformed as: R' = (R - average(R)) / std(R).
    If a rate is recorded at a single field, std(R) is replaced by 1

    Parameters
    ----------
    data : TYPE: dictionnary
        DESCRIPTION: contains the relaxation rates and errors.
    aa : TYPE: str
        DESCRIPTION: residue.
    perform_scaling : TYPE: Bool
        DESCRIPTION: true if scaling has to be performed.

    Returns
    -------
    scaled_rates : TYPE: dictionary
        DESCRIPTION: contains the scaled rates and errors.
    average : TYPE: float
        DESCRIPTION: average of the recorded rates over multiple fields. Equals 0 when perform_scaling=False.
    std : TYPE: float
        DESCRIPTION: standard deviation of the recorded rates over multiple fields. Equals 1 when perform_scaling=False.

    """
    rate_list = []
    for field in data.keys():
        if aa in data[field].keys() and data[field][aa]['data'] != 'NA':
            rate_list.append(data[field][aa]['data'])
    
    scaled_rates = {field: {} for field in data.keys()}
    
    if len(rate_list) > 1.:
        if perform_scaling == 'True':
            average = np.average(rate_list)
            std = np.std(rate_list)
        else:
            average = 0.
            std = 1.
        
        for field in data.keys():
            if aa in data[field].keys():
                scaled_rates[field] = {'data': (data[field][aa]['data'] - average) / std,
                                        'error': data[field][aa]['error'] / std}
        
        return scaled_rates, average, std
    
    return scaled_rates, 0., 1.


def scale_data(self):
    """
    master function for scaling the data in order to properly perform the MCMC and avoid having some type of data over-weighting the cost function.

    """
    
    self.intensities_scaled, self.HF_data_scaled = {}, {}
    self.average = {'intensities': {}, 'high field': {}}
    self.std = {'intensities': {}, 'high field': {}}
    
    #scale intensities
    for exp in self.intensities.keys():
        self.intensities_scaled[exp] = {}
        self.average['intensities'][exp], self.std['intensities'][exp] = {}, {}
        
        for aa in self.residue_list:
            int_scaled, av, std = scale_intensities(self.intensities[exp][aa], self.perform_scaling)
            
            self.intensities_scaled[exp][aa] = int_scaled
            self.average['intensities'][exp][aa] = av
            self.std['intensities'][exp][aa] = std
            
    # scale high-field data
    for rate_type in self.HF_data.keys():
        self.HF_data_scaled[rate_type] = {}
        self.average['high field'][rate_type], self.std['high field'][rate_type] = {}, {}
        
        for aa in self.residue_list:
            rate_scaled, av, std = scale_rates(self.HF_data[rate_type], aa, self.perform_scaling)
            
            self.HF_data_scaled[rate_type][aa] = rate_scaled
            self.average['high field'][rate_type][aa] = av
            self.std['high field'][rate_type][aa] = std
                

def read_other_input(filename):
    """
    reads the other input file. 

    Parameters
    ----------
    filename : TYPE: str
        DESCRIPTION: path to the file. A multiple-column file containing first the residue, and then as many columns as necessary.
        User need to ensure that the order of the columns is consistent with their use in the C-script files.

    Returns
    -------
    other : TYPE: dictionnary
        DESCRIPTION: contains the other inputs.

    """
    other = {}

    f = open(filename, 'r')
    while True:
        line = read_line(f.readline())

        if len(line) > 0:
            aa = line[0]
            other[aa] = np.empty(shape=(len(line)-1), dtype=float)
            
            for n in range(1, len(line)):
                other[aa][n-1] = float(line[n])
        else:
            f.close()
            return other
        
        
def get_residue_list(self):
    """
    get the residue list. It is obtained from the list of intensity decays.
    
    """
    self.residue_list = []
    
    for exp in self.intensities.keys():
        for aa in self.intensities[exp].keys():
            if aa not in self.residue_list:
                self.residue_list.append(aa)