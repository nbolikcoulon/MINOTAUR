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
from imageio import imread
import os, glob
import corner

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
        
colors = [DarkCyan, "saddlebrown", Cornflowerblue, DarkOrange, Darkgreen, Deeppink, Navy, Mediumvioletred, BlueViolet, Mediumseegreen, Indigo, Mediumvioletred, blue, Crimson]
for i in np.linspace(0, 1, 100):
    colors.append((i, i, i))

    
def PlotFieldProfile(heights, LowerCoefs, MiddleCoefs, HigherCoefs, B0LFields, fields, directoryName):
    yax = np.linspace(heights[0], heights[-1], num=len(heights)*100)
    xax = [FitF.B0Fit(yax[i], LowerCoefs, MiddleCoefs, HigherCoefs) for i in range(len(yax))]
    fig = plt.figure(figsize=(16.18/2,5))
    ax = fig.add_subplot(111)
    plt.xlabel("field (T)", fontsize=15)
    plt.ylabel("heigth (m)", fontsize=15)
    plt.xscale('log')
    for F in B0LFields:
        plt.axvline(x=F, color=Darkgreen)
    ax.plot(xax, yax, '-b', label='Polynomial fit')
    ax.plot(fields, heights, 'k+', label='Measured field')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.legend(loc='best')
    figname = directoryName + "/Fit_of_B0.pdf"
    plt.tight_layout()
    plt.savefig(figname, format='pdf')
    plt.close()
    
    
    
def FigTraj(TrajDirName, sampler, Mean, Labels, Lim, AA):
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
    figname = TrajDirName + "/line-time_Residue" + str(AA) + '.png'
    fig.savefig(figname, format='png')
    plt.close()
            
    
def FigCorr(CorrDirName, samples, Mean, Labels, AA):
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
    figname = CorrDirName + "/Correlations_Residue" + str(AA) + '.png'
    fig.savefig(figname, format='png')
    plt.close()
    



def PlotChi2(dirFigs, AllChi2, AAList):
    ind = np.arange(len(AAList))
    width = 1.0/(len(AAList)+1.0)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(ind+width, AllChi2, width, color = "b")
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
    figname = dirFigs + "/Chi2.pdf"
    plt.savefig(figname, format='pdf')
    plt.close()
               
    
def PlotDynParam(dirFigs, AllParam, ErrParam, parameters, AAList):
    ind = np.arange(len(AAList))
    width = 1.0/(len(AAList)+1.0)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(ind+width, AllParam, width, color = "b", yerr=ErrParam)
    ax.set_ylabel(parameters, fontsize=15)
    ax.set_xlabel('residue number', fontsize=15)
    xTickMarks = [str(i) for i in AAList]
    ax.set_xticks(ind+width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=8)
                        
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    figname = dirFigs + "/" + parameters + ".pdf"
    plt.savefig(figname, format='pdf')
    plt.close()


def PlotR1(self, xB0HF, ResiHF, yRateMeasHF, yRateErrMeasHF, yback, xFields, AA):
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 8])

    ax0 = plt.subplot(gs[0])
    ax0.plot(xFields, [0.0 for i in xFields], color = (0.0, 0.0, 0.501961))

#Plot back-calculated R1
    ax1 = plt.subplot(gs[1], sharex = ax0)
    ax1.plot(xFields, yback, color = colors[0], label='Back calculated R1')
    
#Plot fitted R1 from intensities
    xB0LF = []
    R1LFfittedForCurve = []
    for i in range(len(self.R1LFDataForCurve_Fitted[AA])):
        xB0LF.append(self.R1LFDataForCurve_Fitted[AA][i][0])
        if self.R1LFDataForCurve_Fitted[AA][i][1] != "NA":
            R1LFfittedForCurve.append(self.R1LFDataForCurve_Fitted[AA][i][1])
    ax1.scatter(xB0LF, R1LFfittedForCurve, color = colors[3], label="Fitted R1 at LF")
    
#Plot HF measured R1
    if len(yRateMeasHF) != 0:
        ax0.errorbar(xB0HF, ResiHF, yerr=yRateErrMeasHF, fmt='+', color = colors[4])
        ax1.errorbar(xB0HF, yRateMeasHF, yerr = yRateErrMeasHF, fmt='+', color = colors[4], label="Measured R1 at HF")
        
#parameters for plot
    ax1.legend(loc='best')
    start, end = ax0.get_ylim()
    ax0.yaxis.set_ticks(np.array([-max(abs(start), abs(end)), 0.0, max(abs(start), abs(end))]))
    plt.xscale('log')
    plt.xlim([min(xB0LF)-0.1, 25.0])
    axes = fig.get_axes()
    axes[0].set_ylabel("residuals")
    axes[1].set_ylabel(r"$R_1$ ($s^{-1}$)", fontsize=15)
    ax0.set_title(r"$R_1$ Residue " + str(self.AAList[AA]), fontsize=18)
    plt.xlabel("field (T)", fontsize=15)
    ax0.xaxis.set_ticks_position('bottom')
    ax0.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    fig.tight_layout(h_pad=0.0)
    Figname = self.directoryName + "/FitAllResidues/R1/R1_Residue" + str(self.AAList[AA])
    plt.savefig(Figname + ".pdf", format='pdf')
    plt.savefig(Figname + ".png", format='png')
    plt.close()



def PlotRate(self, xB0HF, ResiHF, yRateMesHF, yRateErrMesHF, RelaxType, yback, xFields, AA):
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 8])
                    
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1], sharex = ax0)
                    
    ax0.plot(xFields, [0.0 for i in xFields], color = (0.0, 0.0, 0.501961))

#Plot back-calculated rate
    ax1.plot(xFields, yback, color = colors[0], label='Back calculated ' + str(RelaxType))
    
#Plot HF measured rate
    if len(yRateMesHF) != 0:
        ax0.errorbar(xB0HF, ResiHF, yerr=yRateErrMesHF, fmt='+', color = (0.0, 0.392157, 0.0))
        ax1.errorbar(xB0HF, yRateMesHF, yerr = yRateErrMesHF, fmt='+', color = colors[1], label="Measured " + str(RelaxType) + " at HF")
        
#parameters for plot
    ax1.legend(loc='best')
    start, end = ax0.get_ylim()
    ax0.yaxis.set_ticks(np.array([-max(abs(start), abs(end)), 0.0, max(abs(start), abs(end))]))
    plt.xlabel("field (T)", fontsize=15)
    # if len(xB0HF) != 0:
    #     plt.xlim([min(min(xB0HF)*0.9, 8.0), max(max(xB0HF)*1.1, 25.0)])
    # else:
    #     plt.xlim([8.0, 25.0])
        # 
    # minY = min(yback)
    # maxY = max(yback)
    # plt.ylim(minY*0.8, maxY*1.2)
        
    axes = fig.get_axes()
    axes[0].set_ylabel("residuals")
    axes[1].set_ylabel(str(RelaxType) + r" ($s^{-1}$)", fontsize=15)
    fig.tight_layout(h_pad=0.0)
    ax0.set_title(str(RelaxType) + " of residue " + str(self.AAList[AA]), fontsize=18)
    ax0.xaxis.set_ticks_position('bottom')
    ax0.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    Figname = self.directoryName + "/FitAllResidues/" + str(RelaxType) + "/Measured" + str(RelaxType) + "_VS_BackCalculated" + str(RelaxType) + "_Residue" + str(self.AAList[AA])
    plt.savefig(Figname + ".pdf", format='pdf')
    plt.savefig(Figname + ".png", format='png')
    plt.close()
                    

    
    
    
    
    
def PlotIntensities(self, BackIntensities, FigIntensitiesFolder, IntensitiesForPlot, IntensitiesErrForPlot, DelaysForPlot, DelaysForPlot2, AA):
    for LF in range(len(self.B0LFields[AA])):
        figname = FigIntensitiesFolder + "/Residue" + str(self.AAList[AA]) + "_" + str(round(self.B0LFields[AA][LF], 2))
        
        t = 0
        n = 2
        while t == 0:
            if os.path.exists(figname + ".ps"):
                figname = figname + "_" + str(n)
                n += 1
            else:
                t = 1

        
        fig, ax = plt.subplots()
        ax.plot(DelaysForPlot2[LF], BackIntensities[LF], color = colors[1], label="Simulated intensity decay")
        ax.errorbar(DelaysForPlot[LF], IntensitiesForPlot[LF], yerr = IntensitiesErrForPlot[LF], fmt = 'o', color = colors[0], label="Measured intensities")
        
        ax.legend(loc='best')
        plt.title("Intensity decay for residue " + str(self.AAList[AA]) + " at " + str(round(self.B0LFields[AA][LF], 2)) + " T", fontsize=18)
        plt.xlabel("VC time (s)", fontsize=15)
        plt.ylabel("Intensity", fontsize=15)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        plt.savefig(figname + ".pdf", format='pdf')
        plt.savefig(figname + ".png", format='png')
        plt.close()
    
    
    
    
    
def plotImage(folder, file):
    im = imread(os.path.join(folder, file)).astype(np.float32) / 255
    plt.imshow(im)
    a = plt.gca()
    a.get_xaxis().set_visible(False)
    a.get_yaxis().set_visible(False)
    
    
    
    
def Convert(Folder, OutputName, Format, Remove):
    pp = PdfPages(Folder + "/" + OutputName)
    files = glob.glob(Folder + '/*.' + Format)
    plt.clf()
    for i in range(len(files)):
        plotImage(Folder, files[i])
        pp.savefig(plt.gcf())
        plt.clf()

    pp.close()
        
    if Remove:
        for f in files:
            os.remove(f)

