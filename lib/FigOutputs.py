#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################################
#                                                    #
#                   Figure outputs                   #
#       - Relaxation rates profiles                  #
#       - Relaxation dispersion plots                #
#       - Dynamics parameters                        #
#                                                    #
######################################################
import FitFunctions as FitF

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
import corner

import ShuttlingSimulation as ShSim
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

    
def PlotFieldProfile(field_cal, LowerCoefs, MiddleCoefs, HigherCoefs, B0LFields, directoryName):
    yax = np.linspace(min(field_cal.keys()), max(field_cal.keys()), num=len(field_cal.keys())*100)
    xax = [FitF.B0Fit(yax[i], LowerCoefs, MiddleCoefs, HigherCoefs) for i in range(len(yax))]
    fig = plt.figure(figsize=(16.18/2,5))
    ax = fig.add_subplot(111)
    plt.xlabel("field (T)", fontsize=15)
    plt.ylabel("heigth (m)", fontsize=15)
    plt.xscale('log')
    for exp in B0LFields.keys():
        plt.axvline(x=B0LFields[exp], color=Darkgreen)
    ax.plot(xax, yax, '-b', label='Polynomial fit')
    ax.scatter(field_cal.values(), field_cal.keys(), label='Measured field')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.legend(loc='best')
    figname = directoryName + "/Fit_of_B0.pdf"
    plt.tight_layout()
    plt.savefig(figname, format='pdf')
    plt.close()
    
    
def FigTraj(pdf, TrajDirName, sampler, Mean, Labels, Lim, AA):
    nParam = len(Mean)
    
    fig, axes = plt.subplots(nParam, 1, sharex=True)
    for param in range(nParam):
        axes[param].plot(sampler.chain[:, :, param].T, color="b", alpha=0.4)
        axes[param].yaxis.set_major_locator(MaxNLocator(5))
        axes[param].axvline(Lim, color='k')
        axes[param].axhline(Mean[param], color=(0.780392, 0.0823529, 0.521569))
        if param < nParam-1:
            axes[param].set_ylabel(Labels[param])
        else:
            axes[param].set_ylabel("lnf")
        
    fig.suptitle("Residue " + str(AA))
    fig.tight_layout(h_pad=0.0)
    
    pdf.savefig( fig )
    plt.savefig(f'{TrajDirName}/line-time_Residue_{AA}.pdf', format='pdf')
    plt.close()
            
    
def FigCorr(pdf, CorrDirName, samples, Mean, Labels, AA):
    units = [" (ms)", " (us)", " (ns)", " (ps)"]
    labelsCorr = list(Labels)
    
    for param in range(len(ParamFile.Names['OrderParam']), len(Mean)-1-len(ParamFile.Names['others'])):
        if Mean[param] * 100.0 > 1.0 or Mean[param] <= 0.0:
            pass
        else:
            check = 0
            div = 1
            u = -1
            while check == 0:
                div = div*1e3
                u += 1
                if Mean[param] * div > 1 or div == 1e12:
                    labelsCorr[param] = Labels[param] + units[u]
                    for i in range(len(samples)):
                        samples[i][param] = samples[i][param] * div
                    check = 1
                else:
                    pass
            
    fig = corner.corner(samples, color = 'blue', labels=labelsCorr, quantiles = [0.16, 0.5, 0.84], show_titles=True, title_fmt=u'.2f', plot_contours = True)
    fig.gca().annotate("Residue " + str(AA), (1.0, 1.0), xycoords="figure fraction", xytext=(-20, -10), textcoords="offset points", ha="right", va="top")
    
    pdf.savefig( fig )
    plt.savefig(f'{CorrDirName}/Correlations_Residue_{AA}.pdf', format='pdf')
    plt.close()


def PlotChi2(dirFigs, AllChi2):
    AAList = []
    chi2_plot = []
    for aa in AllChi2.keys():
        AAList.append(aa)
        chi2_plot.append(AllChi2[aa])
    
    ind = np.arange(len(AAList))
    width = 1.0/(len(AAList)+1.0)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(ind+width, chi2_plot, width, color = "b")
    ax.set_ylabel(r"$\chi^2$", fontsize=15)
    ax.set_xlabel('residue number', fontsize=15)
    xTickMarks = [str(i) for i in AAList]
    ax.set_xticks(ind+width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=8)
                    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.savefig(f'{dirFigs}/Chi2.pdf', format='pdf')
    plt.close()
               
    
def PlotDynParam(dirFigs, AllParam, param, param_count):
    AAList = []
    mean_param = []
    err_param = []
    for aa in AllParam.keys():
        AAList.append(aa)
        mean_param.append(AllParam[aa][0][param_count])
        err_param.append((AllParam[aa][1][param_count] + AllParam[aa][2][param_count]) / 2.)
    
    ind = np.arange(len(AAList))
    width = 1.0/(len(AAList)+1.0)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(ind+width, mean_param, width, color = "b", yerr=err_param)
    ax.set_ylabel(param, fontsize=15)
    ax.set_xlabel('residue number', fontsize=15)
    xTickMarks = [str(i) for i in AAList]
    ax.set_xticks(ind+width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=8)
                        
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.savefig(f'{dirFigs}/{param}.pdf', format='pdf')
    plt.close()
    
    
def Format_Rates_Plot(data_dict, err_dict, x_dict, AA):
    f_x = []
    f_data = []
    f_err = []
    
    if x_dict == 'keys':
        for x in data_dict.keys():
            try:
                f_data.append(data_dict[x][AA])
                f_err.append(err_dict[x][AA])       
                f_x.append(x)
            except:
                pass
    else:
        for x in data_dict.keys():
            try:
                f_data.append(data_dict[x][AA])
                f_err.append(err_dict[x][AA])
                f_x.append(x_dict[x][AA])
            except:
                pass
        
    return f_x, f_data, f_err


def Format_IntensityDecay_Plot(data_dict, err_dict, x_dict):
    f_x = []
    f_data = []
    f_err = []
    
    if x_dict == 'keys':
        for x in data_dict.keys():
            if data_dict[x] != 'NA':
                f_x.append(x)
                f_data.append(data_dict[x])
                f_err.append(err_dict[x])  
    else:
        for x in data_dict.keys():
            if data_dict[x] != 'NA':
                f_x.append(x_dict[x])
                f_data.append(data_dict[x])
                f_err.append(err_dict[x])
        
    return f_x, f_data, f_err


def PlotR1(self, pdf, AA):
    
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 8])

    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1], sharex = ax0)
    
    ax0.axhline(0.0, color = 'k', dashes = [3, 3])
    
    x_data_HF, data_HF, err_HF = Format_Rates_Plot(self.HF_data['R1'], self.HF_err['R1'], 'keys', AA)
    residuals = []
    for c, b0 in enumerate(x_data_HF):
        residuals.append(data_HF[c] - self.RelaxFunc['R1'](b0, self.MCMCparam[AA][0][:-1], self.TauC, self.OtherInputs[AA]))
    
    x_data_LF = list(self.B0LFields.values())
    y_data_LF, y_err_LF = [], []
    for exp in self.set_up.keys():
        time, data, err = Format_IntensityDecay_Plot(self.Intensities[exp][AA], self.Err_Int[exp][AA], 'keys')
        rate_fit, err_fit, pre_exp, pre_exp_err = FitF.fit_intensity_decay(time, data)
        y_data_LF.append(rate_fit)
        y_err_LF.append(err_fit)
        
    if len(x_data_HF) != 0:
        x_back_min, x_back_max = 0.9 * min(x_data_LF), 1.1 * max(x_data_HF)
    else:
        x_back_min, x_back_max = 0.9 * min(x_data_LF), 30.
    x_back = np.linspace(x_back_min, x_back_max, 100)
    y_back = []
    for b0 in x_back:
        y_back.append(self.RelaxFunc['R1'](b0, self.MCMCparam[AA][0][:-1], self.TauC, self.OtherInputs[AA]))
    

#Plot back-calculated R1
    ax1.plot(x_back, y_back, color = colors[0], label='Back calculated R1')
    
#Plot fitted R1 from intensities
    ax1.errorbar(x_data_LF, y_data_LF, yerr = y_err_LF, fmt = '+', color = colors[3], label="R1 from exponential fit of intensity decay")
    
#Plot HF measured R1
    if len(x_data_HF) != 0:
        ax0.errorbar(x_data_HF, residuals, yerr = err_HF, fmt='o', color = colors[4])
        ax1.errorbar(x_data_HF, data_HF, yerr = err_HF, fmt='o', color = colors[4], label="Measured R1 at HF")
        
#parameters for plot
    ax1.legend(loc='best')
    start, end = ax0.get_ylim()
    ax0.yaxis.set_ticks(np.array([-max(abs(start), abs(end)), 0.0, max(abs(start), abs(end))]))
    plt.xscale('log')
    plt.xlim(max(0.01, x_back_min-0.1), x_back_max * 1.1)
    axes = fig.get_axes()
    axes[0].set_ylabel("residuals")
    axes[1].set_ylabel(r"$R_1$ ($s^{-1}$)", fontsize=15)
    ax0.set_title(r"$R_1$ Residue " + str(AA), fontsize=18)
    plt.xlabel("field (T)", fontsize=15)
    ax0.xaxis.set_ticks_position('bottom')
    ax0.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    fig.tight_layout(h_pad=0.0)
    
    pdf.savefig( fig )
    plt.close()


def PlotRate(self, pdf, RelaxRate, AA):
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 8])
                    
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1], sharex = ax0)
                    
    ax0.axhline(0.0, color = 'k', dashes = [3, 3])
    
    x_data_HF, data_HF, err_HF = Format_Rates_Plot(self.HF_data[RelaxRate], self.HF_err[RelaxRate], 'keys', AA)
    residuals = []
    for c, b0 in enumerate(x_data_HF):
        residuals.append(data_HF[c] - self.RelaxFunc[RelaxRate](b0, self.MCMCparam[AA][0][:-1], self.TauC, self.OtherInputs[AA]))
    
    if len(x_data_HF) != 0:
        x_back_min, x_back_max = 0.9 * min(x_data_HF), 1.1 * max(x_data_HF)
    else:
        x_back_min, x_back_max = 9., 30.
    x_back = np.linspace(x_back_min, x_back_max, 100)
    y_back = []
    for b0 in x_back:
        y_back.append(self.RelaxFunc[RelaxRate](b0, self.MCMCparam[AA][0][:-1], self.TauC, self.OtherInputs[AA]))

#Plot back-calculated rate
    ax1.plot(x_back, y_back, color = colors[0], label=f'Back calculated {RelaxRate}')
    
#Plot HF measured rate
    if len(x_data_HF) != 0:
        ax0.errorbar(x_data_HF, residuals, yerr=err_HF, fmt='+', color = (0.0, 0.392157, 0.0))
        ax1.errorbar(x_data_HF, data_HF, yerr = err_HF, fmt='+', color = colors[1], label=f'Measured {RelaxRate} at HF')
        
#parameters for plot
    ax1.legend(loc='best')
    start, end = ax0.get_ylim()
    ax0.yaxis.set_ticks(np.array([-max(abs(start), abs(end)), 0.0, max(abs(start), abs(end))]))
    plt.xlabel("field (T)", fontsize=15)
        
    axes = fig.get_axes()
    axes[0].set_ylabel("residuals")
    axes[1].set_ylabel(RelaxRate + r" ($s^{-1}$)", fontsize=15)
    fig.tight_layout(h_pad=0.0)
    ax0.set_title(f"{RelaxRate} of residue {AA}", fontsize=18)
    ax0.xaxis.set_ticks_position('bottom')
    ax0.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    
    pdf.savefig( fig )
    plt.close()
    
    
def ScalingFactor(SimulatedIntensity, IntensityRelaxometry, ErrIntensity):
    
    ScalingNum = 0.0
    ScalingDen = 0.0
    for vc in range(len(SimulatedIntensity)):
        ScalingNum += IntensityRelaxometry[vc]*SimulatedIntensity[vc]/(ErrIntensity[vc])**2
        ScalingDen += (SimulatedIntensity[vc]/ErrIntensity[vc])**2
            
    return ScalingNum/ScalingDen
    
    
def PlotIntensities(self, FigIntensitiesFolder, AA):
    
    pdf = PdfPages(f'{FigIntensitiesFolder}/Intensity_decays_Residue_{AA}.pdf')
    
    for exp in self.set_up.keys():
        
        time, data, err = Format_IntensityDecay_Plot(self.Intensities[exp][AA], self.Err_Int[exp][AA], 'keys')
        time_scaling = np.asarray(time) + self.set_up[exp]['d22'] + self.set_up[exp]['WTLF']
        rate_fit, err_fit, pre_exp, pre_exp_err = FitF.fit_intensity_decay(time, data)
        
        max_vc = 1.2 * max(self.set_up[exp]['vc'])
        x_back = np.linspace(0., max_vc, 100)
        x_back_calc = x_back + self.set_up[exp]['d22'] + self.set_up[exp]['WTLF']
        y_back = ShSim.PropCalDiag_ForPlot(np.asarray(self.MCMCparam[AA][0][:-1]), self.TauC, self.OtherInputs[AA], self.MagField, self.Increment,
                                           self.shuttling_fields, self.shuttling_delays, self.set_up, self.B0LFields, ParamFile.PositionAuto, x_back_calc, exp)
        y_back_exp_fit = FitF.exp(x_back, rate_fit, pre_exp)
        
        back_for_scaling = ShSim.PropCalDiag_ForPlot(np.asarray(self.MCMCparam[AA][0][:-1]), self.TauC, self.OtherInputs[AA], self.MagField, self.Increment,
                                                     self.shuttling_fields, self.shuttling_delays, self.set_up, self.B0LFields, ParamFile.PositionAuto, time_scaling, exp)
        scaling = ScalingFactor(back_for_scaling, data, err)

        
        fig, ax = plt.subplots()
        ax.plot(x_back, scaling * y_back, color = colors[1], label="Intensity decay from MCMC simulation")
        ax.plot(x_back, y_back_exp_fit, color = colors[2], dashes = [3, 3], label="Exponential fit")
        ax.errorbar(time, data, yerr = err, fmt = 'o', color = colors[0], label="Measured intensities")
        
        ax.legend(loc='best')
        plt.title(f"Intensity decay for residue {AA} at {round(self.B0LFields[exp], 2)} T", fontsize=18)
        plt.xlabel("VC time (s)", fontsize=15)
        plt.ylabel("Intensity", fontsize=15)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
        pdf.savefig( fig )
        plt.close()
        
    pdf.close()
