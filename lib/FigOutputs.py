#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################################
#                                                    #
#                   Figure outputs                   #
#                                                    #
######################################################
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
import corner

import ShuttlingSimulation as ShSim
import FitFunctions as FitF
import Parameters as ParamFile

        
Navy = (0.0, 0.0, 0.501961)
Deeppink = (1.0, 0.0784314, 0.576471)
DarkCyan = (0.0, 0.545098, 0.545098)
Cornflowerblue = (0.392157, 0.584314, 0.929412)
BlueViolet = (0.541176, 0.168627, 0.886275)
Mediumvioletred = (0.780392, 0.0823529, 0.521569)
DarkOrange = (1.0, 0.54902, 0.0)
Crimson = (0.862745, 0.0784314, 0.235294)
Mediumseegreen = (0.235294, 0.701961, 0.443137)
Indigo = (0.294118, 0.0, 0.509804)
Darkgreen = (0.0, 0.392157, 0.0)
blue = (0.0, 0.0, 1.0)
        
colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf',
          "saddlebrown", Cornflowerblue, DarkOrange, Darkgreen, Deeppink, Navy, Mediumvioletred, BlueViolet, Mediumseegreen, Indigo, Mediumvioletred, blue, Crimson]
for i in np.linspace(0, 1, 100):
    colors.append((i, i, i))

    
def plot_field_profile(field_cal, set_up, B0_cal_coeff, B0_low_fields, tunnel_position, tunnel_field, directory_name):
    """
    plots the calibration field profile, and indicates the low fields at which experiments are recorded

    Parameters
    ----------
    field_cal : TYPE: dictionnary
        DESCRIPTION: field vs height from the calibration file.
    set_up : TYPE: dictionnary
        DESCRIPTION: relaxometry setup info
    B0_cal_coeff : TYPE: dictionnary
        DESCRIPTION: result of the fit of the field calibration.
    B0_low_fields : TYPE: dictionnary
        DESCRIPTION: low fields at which experiments are recorded.
    tunnel_position: TYPE: float
        DESCRIPTION: position at which the tunnel starts.
    tunnel_field : TYPE: float
        DESCRIPTION: field of the tunnel.
    directory_name : TYPE: str
        DESCRIPTION: output directory for the figure.
        
    """
    h_plot = list(filter(lambda x: x < tunnel_position, field_cal.keys()))
    dist_plot = [FitF.Calc_B0(h, B0_cal_coeff, tunnel_position, tunnel_field) for h in h_plot]

    fig = plt.figure(figsize=(16.18/2,5))
    ax = fig.add_subplot(111)
    
    ax.plot(dist_plot, h_plot, '-b', label='Polynomial fit')
    ax.scatter(field_cal.values(), field_cal.keys(), label='Measured field')
    for exp in B0_low_fields.keys():
        if 'field' not in set_up[exp].keys():
            ax.axvline(x=B0_low_fields[exp], color=Darkgreen)
        else:
            ax.axvline(x=B0_low_fields[exp], color=Mediumseegreen)
    
    ax.set_xlabel("field (T)", fontsize=15)
    ax.set_ylabel("heigth (m)", fontsize=15)
    ax.set_xscale('log')
    plt.legend(loc='best')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')    
    plt.savefig(f'{directory_name}/Fit_of_B0.pdf', format='pdf')
    plt.close()
    
    
def figure_trajectory(pdf, traj_directory, sampler, mean, labels, cutoff, AA):
    """
    plots the MCMC trajectory.

    Parameters
    ----------
    pdf : TYPE: str
        DESCRIPTION: name of the pdf file containing the trajectory for all the residues..
    traj_directory : TYPE: str
        DESCRIPTION: directory where individual trajectory figures are saved.
    sampler : TYPE: array
        DESCRIPTION: MCMC trajectory.
    mean : TYPE: dictionnary
        DESCRIPTION: mean value of each parameters.
    labels : TYPE: list
        DESCRIPTION: name of parameters.
    cutoff : TYPE: int
        DESCRIPTION: limit before which MCMC points are not considered.
    AA : TYPE: str
        DESCRIPTION: residue.
        
    """
    nParam = len(mean)
    
    fig, axes = plt.subplots(nParam, 1, sharex=True)
    for param in range(nParam):
        axes[param].plot(sampler.chain[:, :, param].T, color='b', alpha=0.4)
        axes[param].yaxis.set_major_locator(MaxNLocator(5))
        axes[param].axvline(cutoff, color='k')
        axes[param].axhline(mean[param], color='r')
        if param < nParam-1:
            axes[param].set_ylabel(labels[param])
        else:
            axes[param].set_ylabel("lnf")
        
    fig.suptitle(f"Residue {AA}")
    fig.tight_layout(h_pad=0.0)
    
    pdf.savefig( fig )
    plt.savefig(f'{traj_directory}/line-time_Residue_{AA}.pdf', format='pdf')
    plt.close()
            
    
def figure_correlation_plot(pdf, correlation_directory, samples, mean, labels, AA):
    """
    plots the MCMC correlation plots.

    Parameters
    ----------
    pdf : TYPE: str
        DESCRIPTION: name of the pdf file containing the correlation plot for all the residues...
    correlation_directory : TYPE: str
        DESCRIPTION.: directory where individual correlation plots are saved.
    samples : TYPE: array
        DESCRIPTION: MCMC trajectory.
    mean : TYPE: dictionnary
        DESCRIPTION: mean value of each parameters..
    labels : TYPE: list
        DESCRIPTION: name of parameters..
    AA : TYPE: str
        DESCRIPTION: residue.
        
    """
    units = [" (ms)", " (us)", " (ns)", " (ps)"]
    labelsCorr = list(labels)
    
    for param in range(len(ParamFile.Names['OrderParam']), len(mean)-1-len(ParamFile.Names['others'])):
        if mean[param] * 100.0 > 1.0 or mean[param] <= 0.0:
            pass
        else:
            div = 1
            u = -1
            while True:
                div = div*1e3
                u += 1
                if mean[param] * div > 1 or div == 1e12:
                    break
            labelsCorr[param] = labels[param] + units[u]
            for i in range(len(samples)):
                samples[i][param] = samples[i][param] * div
            
    fig = corner.corner(samples, color = 'blue', labels=labelsCorr, quantiles = [0.16, 0.5, 0.84], show_titles=True, title_fmt=u'.2f', plot_contours = True)
    fig.gca().annotate(f"Residue {AA}" , (1.0, 1.0), xycoords="figure fraction", xytext=(-20, -10), textcoords="offset points", ha="right", va="top")
    
    pdf.savefig( fig )
    plt.savefig(f'{correlation_directory}/Correlations_Residue_{AA}.pdf', format='pdf')
    plt.close()
    
    
def bar_plot(figname, xaxis, data, error = False):
    """
    master function to create bar plots.

    Parameters
    ----------
    figname : TYPE: str
        DESCRIPTION: figure name.
    xaxis : TYPE: array
        DESCRIPTION: x-values.
    data : TYPE: array
        DESCRIPTION: y-values.
    error : TYPE, optional: array
        DESCRIPTION. y-value errors. The default is False.
        
    """
    ind = np.arange(len(xaxis))
    width = 1.0/(len(xaxis)+1.0)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.bar(ind+width, data, width, color = "b", yerr=error)
    
    ax.set_ylabel(r"$\chi^2$", fontsize=15)
    ax.set_xlabel('residue number', fontsize=15)
    xTickMarks = [str(i) for i in xaxis]
    ax.set_xticks(ind+width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=8)
                    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.savefig(figname, format='pdf')
    plt.close()


def plot_chi2(figure_directory, All_Chi2):
    """
    bar plots of the chi2 along the protein sequence.

    Parameters
    ----------
    figure_directory : TYPE: str
        DESCRIPTION: figure directory.
    All_Chi2 : TYPE: dictionnary
        DESCRIPTION: contains the chi2 value for each residue.
        
    """
    residue_list = list(All_Chi2.keys())
    chi2_plot = list(All_Chi2.values())
    
    bar_plot(f'{figure_directory}/Chi2.pdf', residue_list, chi2_plot)
    
    
def plot_dynamic_parameters(figure_directory, all_param, param, param_count):
    """
    bar plots of the MCMC parameters along the protein sequence.

    Parameters
    ----------
    figure_directory : TYPE: str
        DESCRIPTION: figure directory.
    all_param : TYPE: dictionnary
        DESCRIPTION: contains the MCMC results.
    param : TYPE: str
        DESCRIPTION: name of the parameter for which the bar plot will be created.
    param_count : TYPE: int
        DESCRIPTION: counter for param in arrays.
        
    """
    residue_list = list(all_param.keys())
    mean_param = [all_param[aa]['Mean'][param_count] for aa in residue_list]
    err_param = [(all_param[aa]['+'][param_count] + all_param[aa]['-'][param_count]) / 2. for aa in residue_list]
    
    bar_plot(f'{figure_directory}/{param}.pdf', residue_list, mean_param, err_param)
    
    
def format_data_plot(data_dict, AA):
    """
    data formating before plotting.
    Valid both for high-field rates and intensity decays, hense the check on the dictionnary keys.

    Parameters
    ----------
    data_dict : TYPE: dictionnary
        DESCRIPTION: data to plot.
    AA : TYPE: str
        DESCRIPTION: residue of interest.

    Returns
    -------
    f_x : TYPE: array
        DESCRIPTION: x-data.
    f_data : TYPE:array
        DESCRIPTION: y-data.
    f_err : TYPE: array
        DESCRIPTION: y-error.

    """
    f_x = []
    f_data = []
    f_err = []
    
    if AA in data_dict.keys():
        for x in data_dict[AA].keys():
            try:
                f_data.append(float(data_dict[AA][x]['data']))
                f_err.append(float(data_dict[AA][x]['error']))
                f_x.append(x)
            except:
                pass
    else:
        for x in data_dict.keys():
            try:
                f_data.append(float(data_dict[x][AA]['data']))
                f_err.append(float(data_dict[x][AA]['error']))
                f_x.append(x)
            except:
                pass
        
    return f_x, f_data, f_err


def plot_R1(self, AA):
    """
    plots the R1, both high-fields and low-fields.

    Parameters
    ----------
    AA : TYPE: str
        DESCRIPTION: residue.
        
    """
    x_data_HF, data_HF, err_HF = format_data_plot(self.HF_data['R1'], AA)
    residuals = FitF.residuals(x_data_HF, data_HF, self.relaxation_function['R1'], self.MCMC_param[AA]['Mean'][:-1], self.TauC, self.other_inputs[AA])
    
    x_data_LF = list(self.B0_low_field.values())
    y_data_LF, y_err_LF = [], []
    for exp in self.set_up.keys():
        time, data, err = format_data_plot(self.intensities[exp], AA)
        rate_fit, err_fit, _, _ = FitF.fit_intensity_decay(time, data)
        y_data_LF.append(rate_fit)
        y_err_LF.append(err_fit)
        
    x_back_min = np.log(0.9 * min(x_data_LF)) / np.log(10.)
    x_back_max = np.log(1.1 * max(x_data_HF)) / np.log(10.) if len(x_data_HF) != 0 else np.log(30.) / np.log(10.)
    x_back = np.logspace(x_back_min, x_back_max, 100)
    y_back = [self.relaxation_function['R1'](b0, self.MCMC_param[AA]['Mean'][:-1], self.TauC, self.other_inputs[AA]) for b0 in x_back]
        
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 8])
    ax_residuals = plt.subplot(gs[0])
    ax_rate = plt.subplot(gs[1], sharex = ax_residuals)
    
    ax_residuals.axhline(0.0, color = 'k', dashes = [3, 3])

    ax_rate.plot(x_back, y_back, color = colors[0], label='Back calculated R1')
    ax_rate.errorbar(x_data_LF, y_data_LF, yerr = y_err_LF, fmt = '+', color = colors[3], label="R1 from exponential fit of intensity decays")
    if len(x_data_HF) != 0:
        ax_residuals.errorbar(x_data_HF, residuals, yerr = err_HF, fmt='o', color = colors[4])
        ax_rate.errorbar(x_data_HF, data_HF, yerr = err_HF, fmt='o', color = colors[4], label="Measured R1 at HF")
        
    ax_rate.legend(loc='best')
    start, end = ax_residuals.get_ylim()
    ax_residuals.yaxis.set_ticks(np.array([-max(abs(start), abs(end)), 0.0, max(abs(start), abs(end))]))
    plt.xscale('log')
    plt.xlim(10**x_back_min, 10**x_back_max)
    ax_residuals.set_ylabel("residuals")
    ax_rate.set_ylabel(r"$R_1$ ($s^{-1}$)", fontsize=15)
    ax_residuals.set_title(r"$R_1$ Residue " + str(AA), fontsize=18)
    plt.xlabel("field (T)", fontsize=15)
    ax_residuals.xaxis.set_ticks_position('bottom')
    ax_residuals.yaxis.set_ticks_position('left')
    ax_rate.xaxis.set_ticks_position('bottom')
    ax_rate.yaxis.set_ticks_position('left')
    fig.tight_layout(h_pad=0.0)
    
    self.pdf_rates['R1'].savefig( fig )
    plt.close()


def plot_rate(self, RelaxRate, AA):
    """
    plots relaxation rates, others than R1.

    Parameters
    ----------
    RelaxRate : TYPE: str
        DESCRIPTION: label of relaxation rate to be plotted.
    AA : TYPE: str
        DESCRIPTION: residue.
        
    """
    x_data_HF, data_HF, err_HF = format_data_plot(self.HF_data[RelaxRate], AA)
    residuals = FitF.residuals(x_data_HF, data_HF, self.relaxation_function[RelaxRate], self.MCMC_param[AA]['Mean'][:-1], self.TauC, self.other_inputs[AA])
    
    x_back_min = 0.9 * min(x_data_HF) if len(x_data_HF) != 0 else 10.
    x_back_max = 1.1 * max(x_data_HF) if len(x_data_HF) != 0 else 30.
    x_back = np.linspace(x_back_min, x_back_max, 100)
    y_back = [self.relaxation_function[RelaxRate](b0, self.MCMC_param[AA]['Mean'][:-1], self.TauC, self.other_inputs[AA]) for b0 in x_back]
        
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 8])
    ax_residuals = plt.subplot(gs[0])
    ax_rate = plt.subplot(gs[1], sharex = ax_residuals)
                    
    ax_residuals.axhline(0.0, color = 'k', dashes = [3, 3])

    ax_rate.plot(x_back, y_back , color = colors[0], label=f'Back calculated {RelaxRate}')
    if len(x_data_HF) != 0:
        ax_residuals.errorbar(x_data_HF, residuals, yerr=err_HF, fmt='+', color = (0.0, 0.392157, 0.0))
        ax_rate.errorbar(x_data_HF, data_HF, yerr = err_HF, fmt='+', color = colors[1], label=f'Measured {RelaxRate} at HF')
        
    ax_rate.legend(loc='best')
    start, end = ax_residuals.get_ylim()
    ax_residuals.yaxis.set_ticks(np.array([-max(abs(start), abs(end)), 0.0, max(abs(start), abs(end))]))
    plt.xlabel("field (T)", fontsize=15)
    ax_residuals.set_ylabel("residuals")
    ax_rate.set_ylabel(RelaxRate + r" ($s^{-1}$)", fontsize=15)
    fig.tight_layout(h_pad=0.0)
    ax_residuals.set_title(f"{RelaxRate} of residue {AA}", fontsize=18)
    ax_residuals.xaxis.set_ticks_position('bottom')
    ax_residuals.yaxis.set_ticks_position('left')
    ax_rate.xaxis.set_ticks_position('bottom')
    ax_rate.yaxis.set_ticks_position('left')
    
    self.pdf_rates[RelaxRate].savefig( fig )
    plt.close()
    
    
def scaling_intensities(simulated_intensity, intensity_relaxometry, intensity_error):
    """
    scaling factor for the simulated intensity decay

    Parameters
    ----------
    simulated_intensity : TYPE: array
        DESCRIPTION: simulated intensity decay.
    intensity_relaxometry : TYPE: array
        DESCRIPTION: experimental intensity decay.
    intensity_error : TYPE: array
        DESCRIPTION: experimental error for the intensities.

    Returns
    -------
    TYPE: float
        DESCRIPTION: scaling factor.

    """
    ScalingNum = 0.0
    ScalingDen = 0.0
    for vc in range(len(simulated_intensity)):
        ScalingNum += intensity_relaxometry[vc]*simulated_intensity[vc] / (intensity_error[vc])**2
        ScalingDen += (simulated_intensity[vc] / intensity_error[vc])**2
            
    return ScalingNum / ScalingDen
    
    
def plot_intensities(self, AA):
    """
    plots intensity decays, along with exponential fit and simulated decay.

    Parameters
    ----------
    AA : TYPE: str
        DESCRIPTION: residue.
        
    """
    pdf = PdfPages(f'{self.dir_output_res}/Intensities/Intensity_decays_Residue_{AA}.pdf')
    
    for exp in self.set_up.keys():
        time, data, err = format_data_plot(self.intensities[exp], AA)
        rate_fit, err_fit, pre_exp, pre_exp_err = FitF.fit_intensity_decay(time, data)
        
        time_for_scaling = np.asarray(time) + self.set_up[exp]['d22'] + self.set_up[exp]['WTLF']
        x_back = np.linspace(0., 1.2 * max(self.set_up[exp]['vc']), 100)
        x_back_calc = x_back + self.set_up[exp]['d22'] + self.set_up[exp]['WTLF']
        

        y_back_exp_fit = FitF.exp(x_back, rate_fit, pre_exp)
        y_back = ShSim.Expected_Values_For_Plot(np.asarray(self.MCMC_param[AA]['Mean'][:-1]), self.TauC, self.other_inputs[AA], self.Static_MagField,
                                           self.shuttling_fields, self.shuttling_delays, self.set_up, self.B0_low_field, ParamFile.PositionAuto, x_back_calc, exp, self.PropFunction)
        back_for_scaling = ShSim.Expected_Values_For_Plot(np.asarray(self.MCMC_param[AA]['Mean'][:-1]), self.TauC, self.other_inputs[AA], self.Static_MagField,
                                                     self.shuttling_fields, self.shuttling_delays, self.set_up, self.B0_low_field, ParamFile.PositionAuto, time_for_scaling, exp, self.PropFunction)
        scaling = scaling_intensities(back_for_scaling, data, err)
        
        fig, ax = plt.subplots()
        ax.plot(x_back, scaling * y_back, color = colors[1], label="Intensity decay from MCMC simulation")
        ax.plot(x_back, y_back_exp_fit, color = colors[2], dashes = [3, 3], label="Exponential fit")
        ax.errorbar(time, data, yerr = err, fmt = 'o', color = colors[0], label="Measured intensities")
        
        ax.legend(loc='best')
        plt.title(f"Intensity decay for residue {AA} at {round(self.B0_low_field[exp], 2)} T", fontsize=18)
        plt.xlabel("VC time (s)", fontsize=15)
        plt.ylabel("Intensity", fontsize=15)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
        pdf.savefig( fig )
        plt.close()
        
    pdf.close()