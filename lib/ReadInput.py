#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################################
#                                                         #
#                       Reads files                       #
#                                                         #
###########################################################
import numpy as np
import sys


def read_line(line):
        line = line.split('\n')[0]
        if '\t' in line:
            line = line.split('\t')
        elif ',' in line:
            line = line.split(',')
        elif ' ' in line:
            line = line.split(' ')
            
        line = list(filter(lambda x: x!='', line))
        
        return line
    
            
def Read_Field_Calibration(path):
    field_cal = {}
    
    f = open(path, 'r')
    while True:
        line = read_line(f.readline())
        
        try:
            field_cal[float(line[0])] = float(line[1])
            
        except:
            f.close()
            return field_cal
        
        
def Read_Exp_Setup(path):
    set_up = {}
    
    f = open(path, 'r')
    while True:
        line = read_line(f.readline())
        
        try:
            exp = int(float(line[0]))
            set_up[exp] = {}
            
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
        
        
def Read_High_Field_Data(path):
    hf_data = {}
    hf_err = {}
    
    f = open(path, 'r')
    while True:
        line = read_line(f.readline())
        
        if len(line) > 0:
            aa = line[0]
            
            try:
                hf_data[aa] = float(line[1])
                hf_err[aa] = float(line[2])
            except:
                pass
            
        else:
            f.close()
            return hf_data, hf_err
        
        
def Read_Relaxometry_Decay(path, delays, exp):
    data = {}
    err = {}
    
    f = open(path, 'r')
    while True:
        line = read_line(f.readline())
        
        if len(line) > 0:
            aa = line[0]
            data[aa] = {}
            err[aa] = {}
            
            if len(line) != int(2*len(delays) + 1):
                    print(len(line))
                    print(2*len(delays) + 1)
                    print(delays)
                    print("")
                    print(f"Experiment number {exp} does not have the correct number of intensities")
                    print("")
                    sys.exit()
                    
            for vc in range(int((len(line) - 1) / 2)):
                try:
                    if delays[vc] in data[aa].keys():
                        data[aa][delays[vc]].append(float(line[1 + 2*vc]))
                        err[aa][delays[vc]].append(float(line[2 + 2*vc]))
                    else:
                        data[aa][delays[vc]] = [float(line[1 + 2*vc])]
                        err[aa][delays[vc]] = [float(line[2 + 2*vc])]
                except:
                    if delays[vc] not in data[aa].keys():
                        data[aa][delays[vc]] = 'NA'
                        err[aa][delays[vc]] = 'NA'
                
        else:
            f.close()
            break
        
    data_av = {}
    err_av = {}
    for aa in data.keys():
        data_av[aa] = {}
        err_av[aa] = {}
        for vc in data[aa].keys():
            if data[aa][vc] != 'NA':
                data_av[aa][vc] = np.average(data[aa][vc])
                err_av[aa][vc] = np.average(err[aa][vc])
            else:
                data_av[aa][vc] = 'NA'
                err_av[aa][vc] = 'NA'
        
    return data_av, err_av
                

def Read_Other_Input(path):
    other = {}

    f = open(path, 'r')
    while True:
        line = read_line(f.readline())

        if len(line) > 0:
            aa = line[0]
            other[aa] = []
            
            for n in range(1, len(line)):
                other[aa].append(float(line[n]))
            other[aa] = np.asarray(other[aa])
            
        else:
            f.close()
            break
    
    return other