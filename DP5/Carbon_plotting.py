from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
import os
from pathlib import Path

def PlotCarbon(NMRData,Isomers,settings):
    
    xdata = NMRData.carbondata["xdata"]
    
    ydata = NMRData.carbondata["ydata"]
    
    exppeaks = NMRData.carbondata["exppeaks"]

    simulated_ydata = NMRData.carbondata["simulated_ydata"]
    
    removed = NMRData.carbondata["removed"]

    if settings.OutputFolder == '':

        gdir = Path.cwd() / "Graphs" / settings.InputFiles[0] 

    else:

        gdir = settings.OutputFolder / "Graphs" / settings.InputFiles[0]

    for isomerindex,isomer in enumerate(Isomers):

        assigned_peaks = []

        assigned_shifts = []

        assigned_labels = []

        for i,peak in enumerate(isomer.Cexp):

            if peak != '':

                assigned_peaks.append(peak)

                assigned_shifts.append(isomer.Cshifts[i])

                assigned_labels.append(isomer.Clabels[i])

        plt.close()
    
        fig1 = plt.figure(1)
        fig1.set_size_inches(30, 17)
        ax1 = fig1.add_subplot(111)
    
        exppeaks_ppm = xdata[exppeaks].tolist()

        shiftlist = assigned_shifts
    
        totallist = exppeaks_ppm + shiftlist
    
        plt.xlim([max(totallist) + 10, min(totallist) - 10])
    
        plt.plot(xdata, ydata, color='grey', linewidth=0.75, label='experimental spectrum')
        plt.plot(xdata, simulated_ydata, label='simulated spectrum')
    
        plt.xlabel('PPM')  # axis labels
        # plt.yticks([], [])
        plt.title(str(settings.InputFiles[0]) +
                  "\nCarbon NMR of Isomer " + str(isomerindex + 1) + "\n Number of Peaks Found = " + str(len(exppeaks)))
    
        # plot assignments

        for ind1, peak in enumerate(assigned_peaks):

            wh = np.argmin(abs(xdata - peak))
    
            plt.plot([peak, assigned_shifts[ind1]],
                     [ydata[wh], 1.1], linewidth=0.5, color='cyan')

    
        prev = round(exppeaks_ppm[0], 2)
    
        count = 0
    
        # annotate peak locations
    
        for x, txt in enumerate([round(i, 2) for i in exppeaks_ppm]):
    
            if abs(prev - txt) < 5:
    
                count += 1
            else:
                count = 0
                prev = txt
    
            if exppeaks_ppm[x] in assigned_peaks:
                color = 'C1'
            else:
                color = 'grey'
    
            ax1.annotate(txt, (exppeaks_ppm[x], -0.06 - 0.025 * count), color=color, size=10)
    
            plt.plot(exppeaks_ppm[x], ydata[exppeaks[x]], 'o', color=color)
    
        if len(removed) > 0:
            plt.plot(xdata[removed],
                     simulated_ydata[removed], "ro")
    
        # annotate shift positions

        count = 0
    
        ####some quick sorting

        sortshifts = isomer.Cshifts[::-1]

        slabels = isomer.Clabels[::-1]
    
        prev = sortshifts[0]
    
        for x, txt in enumerate(slabels):
    
            if abs(prev - sortshifts[x]) < 4:
                count += 1
            else:
                count = 0
                prev = sortshifts[x]
    
            ax1.annotate(txt, (sortshifts[x], + 2.05 + 0.05 * count), size=18)
    
        simulated_calc_ydata = simulate_calc_data(xdata, isomer.Cshifts, simulated_ydata)
    
        plt.plot(xdata, simulated_calc_ydata + 1.1, label='calculated spectrum')

        plt.legend()

        f_name = "Carbon_" + str(isomerindex + 1 ) + ".svg"

        plt.savefig(gdir / f_name,format = "svg", bbox_inches='tight')

        plt.close()

def simulate_calc_data(spectral_xdata_ppm, calculated_locations, simulated_ydata):
    ###simulate calcutated data

    simulated_calc_ydata = np.zeros(len(spectral_xdata_ppm))

    for peak in calculated_locations:
        y = np.exp(-0.5 * ((spectral_xdata_ppm - peak) / 0.002) ** 2)
        simulated_calc_ydata += y

    scaling_factor = np.amax(simulated_ydata) / np.amax(simulated_calc_ydata)

    simulated_calc_ydata = simulated_calc_ydata * scaling_factor

    return simulated_calc_ydata