from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
import os
from pathlib import Path

def PlotProton(NMRData,Isomers,settings):
    
    xdata = NMRData.protondata["xdata"]
    
    ydata = NMRData.protondata["ydata"]
    
    centres = NMRData.protondata["centres"]
    
    exp_peaks = NMRData.protondata["exppeaks"]
    
    peak_regions = NMRData.protondata["peakregions"]
    
    cummulative_vectors = NMRData.protondata["cummulativevectors"]
    
    integral_sum = NMRData.protondata["integralsum"]
    
    integrals = NMRData.protondata["integrals"]
    
    sim_regions = NMRData.protondata["sim_regions"]

    if settings.OutputFolder == '':

        gdir = Path.cwd() / "Graphs" / settings.InputFiles[0]

    else:

        gdir = settings.OutputFolder / "Graphs" / settings.InputFiles[0]

    for isomerindex,isomer in enumerate(Isomers):

        assigned_shifts = isomer.Hshifts

        assigned_peaks = []

        for peak in isomer.Hexp:

            if peak != '':
                assigned_peaks.append(peak)

        assigned_labels = isomer.Hlabels

        #################################### will probs need to fix sorting here

        plt.close()
        
        fig1 = plt.figure(1)
        
        fig1.set_size_inches(30, 17)
        
        plt.xlim([10, 0])

        plt.xlabel("ppm")
        
        plt.plot(xdata, ydata, label='data', color='grey')
        
        set_exp = sorted(list(set(exp_peaks)))[::-1]
        
        simulate_spectrum(xdata, assigned_shifts, assigned_peaks, set_exp)
        
        plt.axhline(1.05, color='grey')
        
        # plt integral information
        
        prev = 15
        
        count = 0
        
        for index in range(0, len(peak_regions)):
        
            if abs(prev - xdata[centres[index]]) < 0.45:
                count += 1
            else:
                count = 0
                prev = xdata[centres[index]]
        
            plt.annotate(str(integrals[index]) + ' Hs',
                         xy=(xdata[centres[index]], -(0.1) - 0.1 * count), color='C' + str(index % 10), size=18)
        
            plt.annotate(str(round(xdata[centres[index]], 3)) + ' ppm',
                         xy=(xdata[centres[index]], -(0.15) - 0.1 * count), color='C' + str(index % 10), size=18)
        
            plt.plot(xdata[peak_regions[index]], cummulative_vectors[index] + integral_sum[index],
                     color='C' + str(index % 10),
                     linewidth=2)


        for index in range(0, len(peak_regions) - 1):
            plt.plot([xdata[peak_regions[index][-1]], xdata[peak_regions[index + 1][0]]],
                     [integral_sum[index + 1], integral_sum[index + 1]], color='grey')
        
        for index, region in enumerate(peak_regions):
            plt.plot(xdata[region], sim_regions[index], color='C' + str(index % 10))

        ### plotting assignment
        
        plt.yticks([], [])
        plt.title(str(settings.InputFiles[0]) +
                  "\nProton NMR of Isomer " + str(isomerindex + 1) + "\n Number of Peaks Found = " + str(len(exp_peaks)))
        
        # plot assignments
        
        for ind1, peak in enumerate(assigned_peaks):
            plt.plot([peak, assigned_shifts[ind1]],
                     [1, 1.05], linewidth=0.5, color='cyan')

        # annotate peak locations
        
        for x, txt in enumerate(exp_peaks):
        
            if exp_peaks[x] in assigned_peaks:

                color = 'C1'

            else:

                color = 'grey'
        
            plt.plot(txt, -0.02, 'o', color=color)
        
        # annotate shift positions
        
        prev = 0
        
        count = 0

        s = np.argsort(np.array(assigned_shifts))

        s_assigned_shifts = np.array(assigned_shifts)[s]

        s_assigned_labels = np.array(assigned_labels)[s]

        s_assigned_peaks = np.array(assigned_peaks)[s]

        
        for x, txt in enumerate(s_assigned_labels[::-1]):

            w = np.where( set_exp == s_assigned_peaks[::-1][x])[0][0]

            color = w % 10
        
            if abs(prev - s_assigned_shifts[::-1][x]) < 0.2:
                count += 1

            else:
                count = 0
                prev = s_assigned_shifts[::-1][x]
                
            plt.annotate(txt, (s_assigned_shifts[::-1][x], + 1.25 + 0.05 * count), size=18,color = 'C' + str(color))

        plt.ylim([-0.5, 2.0])

        f_name = "Proton_" + str(isomerindex + 1) + ".svg"

        plt.savefig(gdir /f_name,format = "svg", bbox_inches='tight')

        plt.close()

def simulate_spectrum(spectral_xdata_ppm,calc_shifts,assigned_peaks,set_exp):

    for ind,shift in enumerate(calc_shifts):

        exp_p = assigned_peaks[ind]

        ind2 = set_exp.index(exp_p)

        y = lorentzian(spectral_xdata_ppm,0.001,shift,0.2)

        plt.plot(spectral_xdata_ppm,y+1.05,color = 'C' + str(ind2 % 10))


def lorentzian(p, w, p0, A):
    x = (p0 - p) / (w / 2)
    L = A / (1 + x ** 2)

    return L
