import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt
from scipy.optimize import linear_sum_assignment as optimise
from scipy.stats import linregress
import copy
import os
import pickle
from openbabel import *

def AssignProton(NMRData,Isomers,settings):

    #do the assignment

    for isomer in Isomers:

        assigned_shifts, assigned_peaks, assigned_labels,scaled_shifts \
            = iterative_assignment(NMRData.Hshifts,isomer.Hshifts,isomer.Hlabels,NMRData.protondata["integrals"],settings)

    #add to isomers instance

        isomer.Hexp = [''] * len(isomer.Hshifts)

        for label,peak in zip(assigned_labels,assigned_peaks):

            w = isomer.Hlabels.index(label)

            isomer.Hexp[w] = peak

    return Isomers

def iterative_assignment(exp_peaks,calculated_shifts, H_labels,rounded_integrals,settings):

    exp_peaks = np.array(exp_peaks)

    calculated_shifts = np.array(calculated_shifts)

    H_labels = np.array(H_labels)

    lnum = 0

    new_assigned_shifts = []
    old_assigned_shifts = [1]

    #print("calc shifts",calculated_shifts)

    while old_assigned_shifts != new_assigned_shifts:

        if lnum ==0:

            scaled_shifts = external_scale_proton_shifts(calculated_shifts)

        else:
            old_assigned_shifts = copy.copy(new_assigned_shifts)
            old_assigned_peaks = copy.copy(new_assigned_peaks)

            scaled_shifts,slope,intercept = internal_scale_proton_shifts(old_assigned_shifts,old_assigned_peaks,calculated_shifts)

            #scaled_mu,scaled_std = scale_params(slope,intercept)


        ###############assign methyl groups first

        #find methyl groups

        m_protons = methyl_protons(settings.InputFiles[0] + ".sdf")

        m_shifts = np.array([])

        # find the average shifts of these groups

        for m_group in m_protons:

            s = 0

            for proton in m_group:

                w = np.where(H_labels == proton)

                s += scaled_shifts[w]/3

            m_shifts = np.hstack((m_shifts,s))

        #find peaks these can be assigned too

        methyl_peaks = []

        rounded_integrals = np.array(rounded_integrals)

        w = (rounded_integrals - (rounded_integrals % 3)) // 3

        for ind, peak in enumerate(sorted(list(set(exp_peaks)))[::-1]):
            methyl_peaks += [peak] * w[ind]

        #create difference matrix

        diff_matrix = np.zeros((len(m_shifts),len(methyl_peaks)))

        for ind1, i in enumerate(m_shifts):
            for ind2, j in enumerate(methyl_peaks):
                diff_matrix[ind1,ind2] = j-i

        prob_matrix = proton_probabilities(diff_matrix,settings)

        prob_matrix = prob_matrix**2

        prob_matrix = 1 - prob_matrix

        vertical_ind, horizontal_ind = optimise(prob_matrix)

        #unpack this assignment

        opt_labelsm = []

        opt_shiftsm = []

        opt_peaksm = []

        for j in vertical_ind:

            opt_labelsm.extend(m_protons[j])

        for i in horizontal_ind:

            opt_peaksm += 3*[methyl_peaks[i]]

        for label in opt_labelsm:

            w = np.where(H_labels == label)

            opt_shiftsm.append(calculated_shifts[w][0])

        #remove shifts/peaks/labels for the list to assign

        calculated_shiftsp = copy.copy(calculated_shifts)

        exp_peaksp = copy.copy(exp_peaks)

        scaled_shiftsp = copy.copy(scaled_shifts)

        H_labelsp = copy.copy(H_labels)

        #peaks

        for p in opt_peaksm:

            w = np.where(exp_peaksp == p)[0][0]

            exp_peaksp = np.delete(exp_peaksp,w)

        #shifts

        for s in opt_shiftsm:

            w = np.where(calculated_shiftsp == s)[0][0]

            calculated_shiftsp = np.delete(calculated_shiftsp,w)
            scaled_shiftsp = np.delete(scaled_shiftsp, w)

        #labels

        for l in opt_labelsm:

            w = np.where(H_labelsp == l)[0][0]

            H_labelsp = np.delete(H_labelsp,w)

        ###############assigned everything else

        diff_matrix = np.zeros((len(calculated_shiftsp),len(exp_peaksp)))

        for ind1, i in enumerate(scaled_shiftsp):
            for ind2, j in enumerate(exp_peaksp):
                diff_matrix[ind1,ind2] = j-i

        prob_matrix = proton_probabilities(diff_matrix,settings)

        prob_matrix = prob_matrix**2

        prob_matrix = 1 - prob_matrix

        vertical_ind, horizontal_ind = optimise(prob_matrix)

        opt_peaksp =  exp_peaksp[horizontal_ind]

        opt_shiftsp = calculated_shiftsp[vertical_ind]

        opt_labelsp = H_labelsp[vertical_ind]

        ################ combine these assignments

        opt_peaks = np.hstack((opt_peaksm,opt_peaksp))

        opt_shifts = np.hstack((opt_shiftsm,opt_shiftsp))

        opt_labels = np.hstack((opt_labelsm,opt_labelsp))

        #check for any shifts that have not been assigned

        copyshifts = list(copy.copy(calculated_shifts))
        copylabels = list(copy.copy(H_labels))

        for shift,label in zip(opt_shifts,opt_labels):

            copyshifts.remove(shift)
            copylabels.remove(label)

        #assign these to the closest peaks - regardless of integrals

        for shift,label in zip(copyshifts,copylabels):

            mindiff = np.array(exp_peaks - shift).argmin()

            opt_peaks = np.append(opt_peaks,exp_peaks[mindiff])

            opt_labels = np.append(opt_labels,label)

            opt_shifts = np.append(opt_shifts,shift)

        opt_shifts,opt_peaks,opt_labels = removecrossassignments(opt_peaks, opt_shifts, opt_labels)

        #### sort output wrt original H labels

        indv = []

        for label in opt_labels:

            wh = np.where(H_labels == label)

            indv.append(wh[0][0])

        ind = np.argsort(opt_shifts)[::-1]

        assigned_shifts = opt_shifts[indv]

        assigned_peaks = opt_peaks[indv]

        assigned_labels = opt_labels[indv]

        ind = np.argsort(assigned_shifts)

        assigned_shifts = assigned_shifts[ind].tolist()
        assigned_peaks = assigned_peaks[ind].tolist()
        assigned_labels = assigned_labels[ind].tolist()

        lnum += 1

        new_assigned_shifts =copy.copy(assigned_shifts)
        new_assigned_peaks=copy.copy(assigned_peaks)

    return assigned_shifts , assigned_peaks, assigned_labels, scaled_shifts

def external_scale_proton_shifts(calculated_shifts):

    scaled = 0.9770793502768845 * calculated_shifts - 0.019505417520415236

    return scaled

def internal_scale_proton_shifts(assigned_shifts,assigned_peaks,calculated_shifts):

    slope, intercept, r_value, p_value, std_err = linregress(assigned_shifts, assigned_peaks)

    scaled_shifts = calculated_shifts * slope + intercept

    #plt.close()

    #plt.plot(assigned_shifts,assigned_shifts,color= 'grey')

    #plt.plot(assigned_shifts, assigned_peaks,'ro')

    #plt.plot(np.array(assigned_shifts) *slope + intercept,assigned_peaks,'co')

    #plt.plot(assigned_shifts,np.array(assigned_shifts)*slope + intercept)
    #plt.show()

    return scaled_shifts,slope,intercept

def proton_probabilities(diff_matrix,settings):

    # unscaled error data

    #prob_matrix = norm.pdf(diff_matrix,-0.0935539568345324,0.24204738737722833) / norm.pdf(-0.0935539568345324,-0.0935539568345324,0.24204738737722833)

    #scaled error data

    #prob_matrix = norm.pdf(diff_matrix,0,0.2374164273141329) / norm.pdf(0,0,0.2374164273141329)

    if settings.StatsParamFile != 'none':

        Cmeans, Cstdevs, scaled_mu, scaled_std = ReadParamFile(settings.StatsParamFile,'m')

        prob_matrix = np.ones(np.shape(diff_matrix))

        for mu, std in zip(scaled_mu,scaled_std):

            prob_matrix *= norm.pdf(diff_matrix,mu, std) / norm.pdf(mu, mu, std)

        prob_matrix /= len(scaled_mu)

    else:

        prob_matrix = norm.pdf(diff_matrix, 0, 1 ) / norm.pdf(0, 0, 1)

    return prob_matrix

def ReadParamFile(f, t):
    infile = open(f, 'r')
    inp = infile.readlines()
    infile.close()

    if t not in inp[0]:
        print("Wrong parameter file type, exiting...")
        quit()

    if t == 'm':
        Cmeans = [float(x) for x in inp[1].split(',')]
        Cstdevs = [float(x) for x in inp[2].split(',')]
        Hmeans = [float(x) for x in inp[3].split(',')]
        Hstdevs = [float(x) for x in inp[4].split(',')]

        return Cmeans, Cstdevs, Hmeans, Hstdevs

def simulate_spectrum(spectral_xdata_ppm,calc_shifts):

    y = np.zeros(len(spectral_xdata_ppm))

    for shift in calc_shifts:

        y += lorentzian(spectral_xdata_ppm,0.001,shift,0.2)

    return y

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

def scale_params(slope,intercept):

    plinear_regression_data_exp = np.array([1.73,
                                   1.73,
                                   1.73,
                                   1.83,
                                   3.09,
                                   4.5,
                                   4.81,
                                   6.76,
                                   7.11,
                                   7.13,
                                   7.18,
                                   7.22,
                                   7.32,
                                   1.02,
                                   2.67,
                                   2.67,
                                   2.67,
                                   2.67,
                                   2.67,
                                   3.74,
                                   4.46,
                                   6.8,
                                   6.8,
                                   6.8,
                                   6.99,
                                   7.46,
                                   8.11,
                                   1.27,
                                   1.5,
                                   1.65,
                                   2.27,
                                   2.37,
                                   3.48,
                                   3.51,
                                   3.53,
                                   4.01,
                                   4.35,
                                   5.55,
                                   6.62,
                                   7.16,
                                   7.7,
                                   1.41,
                                   1.57,
                                   1.69,
                                   1.83,
                                   2.04,
                                   2.79,
                                   2.79,
                                   2.79,
                                   2.79,
                                   2.79,
                                   2.79,
                                   2.79,
                                   3.31,
                                   3.31,
                                   4.03,
                                   4.82,
                                   6.32,
                                   7.07,
                                   7.2,
                                   7.2,
                                   7.2,
                                   7.2,
                                   7.2,
                                   7.2,
                                   7.2,
                                   7.2,
                                   3.04,
                                   3.21,
                                   3.78,
                                   3.94,
                                   4.1,
                                   4.31,
                                   5.85,
                                   5.88,
                                   6.15,
                                   6.69,
                                   6.73,
                                   6.85,
                                   7.12,
                                   7.28,
                                   7.6,
                                   7.81,
                                   0.8,
                                   0.93,
                                   0.98,
                                   1.01,
                                   1.09,
                                   1.2,
                                   1.31,
                                   1.44,
                                   1.46,
                                   1.58,
                                   1.61,
                                   1.63,
                                   1.7,
                                   1.85,
                                   1.86,
                                   2.04,
                                   2.08,
                                   2.27,
                                   2.34,
                                   2.37,
                                   2.41,
                                   3.65,
                                   5.74,
                                   1.91,
                                   2.69,
                                   3.07,
                                   3.27,
                                   3.52,
                                   4.16,
                                   5.14,
                                   7.02,
                                   7.14,
                                   7.53,
                                   0.89,
                                   1.36,
                                   1.36,
                                   1.36,
                                   1.36,
                                   1.36,
                                   1.61,
                                   1.61,
                                   3.68,
                                   3.9,
                                   4.08,
                                   4.08,
                                   4.08,
                                   5.69,
                                   8,
                                   0.89,
                                   0.9,
                                   1.29,
                                   1.5,
                                   1.5,
                                   1.5,
                                   1.5,
                                   2.04,
                                   2.15,
                                   2.75,
                                   3.24,
                                   3.34,
                                   3.52,
                                   4.18,
                                   4.2,
                                   5.62,
                                   6.78,
                                   1,
                                   1.08,
                                   2.22,
                                   2.28,
                                   2.7,
                                   3.02,
                                   4.1,
                                   4.33,
                                   4.55,
                                   4.6,
                                   4.63,
                                   4.82,
                                   5.09,
                                   6.09,
                                   6.24,
                                   1.03,
                                   1.05,
                                   2.21,
                                   2.69,
                                   2.69,
                                   3.41,
                                   4.03,
                                   4.13,
                                   4.53,
                                   4.63,
                                   4.63,
                                   4.7,
                                   5.33,
                                   6.16,
                                   6.18,
                                   0.97,
                                   1.22,
                                   1.28,
                                   1.42,
                                   1.42,
                                   1.67,
                                   1.67,
                                   1.67,
                                   1.87,
                                   2.14,
                                   2.53,
                                   2.53,
                                   3.41,
                                   3.6,
                                   4.05,
                                   0.93,
                                   1.23,
                                   1.33,
                                   1.41,
                                   1.61,
                                   1.61,
                                   1.61,
                                   1.74,
                                   1.81,
                                   2.13,
                                   2.51,
                                   2.51,
                                   3.43,
                                   3.88,
                                   3.92,
                                   0.91,
                                   1.15,
                                   1.19,
                                   1.44,
                                   1.6,
                                   1.61,
                                   1.89,
                                   1.89,
                                   1.89,
                                   1.89,
                                   2.17,
                                   2.23,
                                   2.5,
                                   2.94,
                                   3.7,
                                   3.87,
                                   6.23,
                                   6.75,
                                   6.81,
                                   7.33,
                                   7.62,
                                   1.24,
                                   1.27,
                                   1.8,
                                   1.8,
                                   1.8,
                                   1.8,
                                   1.8,
                                   1.8,
                                   2.56,
                                   2.8,
                                   3.16,
                                   3.22,
                                   6.44,
                                   6.86,
                                   1.25,
                                   1.27,
                                   1.61,
                                   1.61,
                                   1.61,
                                   1.61,
                                   1.61,
                                   1.61,
                                   2.86,
                                   2.94,
                                   3.24,
                                   3.58,
                                   6.18,
                                   7.14,
                                   0.98,
                                   1.5,
                                   1.96,
                                   2.66,
                                   3.99,
                                   4.15,
                                   5.05,
                                   7.26,
                                   7.34,
                                   7.39,
                                   0.99,
                                   1.56,
                                   1.97,
                                   2.21,
                                   3.91,
                                   4.17,
                                   5.12,
                                   7.3,
                                   1.01,
                                   1.47,
                                   1.56,
                                   1.64,
                                   1.64,
                                   2.08,
                                   2.45,
                                   4.06,
                                   4.36,
                                   5.26,
                                   7.29,
                                   1,
                                   1.52,
                                   1.52,
                                   1.52,
                                   1.72,
                                   1.92,
                                   2.7,
                                   3.76,
                                   4.27,
                                   4.87,
                                   7.26,
                                   7.34,
                                   7.4,
                                   0.98,
                                   1.05,
                                   1.05,
                                   1.18,
                                   1.42,
                                   1.42,
                                   1.42,
                                   1.42,
                                   1.75,
                                   1.75,
                                   1.86,
                                   1.99,
                                   1.99,
                                   2.41,
                                   3.38,
                                   5.84,
                                   2.03,
                                   2.03,
                                   2.4,
                                   2.4,
                                   3.51,
                                   3.63,
                                   4.11,
                                   6.2,
                                   8.06,
                                   8.34,
                                   1.73,
                                   2.22,
                                   2.33,
                                   3.55,
                                   3.61,
                                   3.77,
                                   4.35,
                                   6.05,
                                   7.62,
                                   2.04,
                                   2.22,
                                   2.52,
                                   3.53,
                                   3.53,
                                   4.23,
                                   4.84,
                                   4.87,
                                   5.36,
                                   7.66,
                                   2.27,
                                   3.16,
                                   3.18,
                                   3.22,
                                   3.27,
                                   3.45,
                                   3.71,
                                   3.97,
                                   4.1,
                                   4.15,
                                   6.8,
                                   7.12,
                                   7.15,
                                   7.2,
                                   7.23,
                                   7.28,
                                   7.59,
                                   0.87,
                                   0.93,
                                   0.96,
                                   0.97,
                                   0.97,
                                   1.06,
                                   1.26,
                                   1.29,
                                   1.52,
                                   1.55,
                                   1.64,
                                   1.76,
                                   1.77,
                                   1.94,
                                   2,
                                   2.04,
                                   2.04,
                                   2.08,
                                   2.2,
                                   2.68,
                                   3.15,
                                   3.2,
                                   3.48,
                                   3.7,
                                   4.13,
                                   4.32,
                                   4.57,
                                   4.82,
                                   7.13,
                                   7.19,
                                   7.46,
                                   8.3,
                                   0.83,
                                   1.28,
                                   1.28,
                                   1.53,
                                   1.74,
                                   1.74,
                                   1.83,
                                   2.18,
                                   2.23,
                                   2.68,
                                   2.81,
                                   2.88,
                                   2.97,
                                   2.97,
                                   3.06,
                                   3.06,
                                   3.12,
                                   3.78,
                                   4.19,
                                   4.19,
                                   4.25,
                                   6.97,
                                   6.97,
                                   6.97,
                                   7.21,
                                   7.21,
                                   7.21,
                                   7.21,
                                   7.21,
                                   7.21,
                                   7.21,
                                   7.21,
                                   7.21,
                                   7.21,
                                   0.99,
                                   1.06,
                                   1.3,
                                   1.36,
                                   1.44,
                                   1.61,
                                   1.63,
                                   1.63,
                                   1.63,
                                   1.74,
                                   2.02,
                                   2.06,
                                   2.23,
                                   2.53,
                                   2.87,
                                   3,
                                   3.06,
                                   3.35,
                                   3.65,
                                   4.05,
                                   4.52,
                                   7.15,
                                   7.71,
                                   7.75,
                                   0.81,
                                   0.9,
                                   1.44,
                                   1.61,
                                   1.82,
                                   2.78,
                                   2.78,
                                   2.92,
                                   2.92,
                                   2.92,
                                   3.1,
                                   3.1,
                                   3.67,
                                   3.67,
                                   3.85,
                                   3.85,
                                   3.85,
                                   3.94,
                                   4.99,
                                   5.63,
                                   6.67,
                                   7.23,
                                   7.23,
                                   7.23,
                                   7.53,
                                   0.84,
                                   0.89,
                                   1.2,
                                   1.55,
                                   1.55,
                                   1.55,
                                   1.55,
                                   1.55,
                                   1.55,
                                   1.55,
                                   1.68,
                                   1.68,
                                   1.68,
                                   1.83,
                                   1.83,
                                   1.83,
                                   1.96,
                                   2.04,
                                   2.04,
                                   2.11,
                                   2.13,
                                   2.28,
                                   2.8,
                                   2.8,
                                   3,
                                   1.19,
                                   1.38,
                                   1.49,
                                   1.5,
                                   1.71,
                                   1.81,
                                   1.88,
                                   1.95,
                                   2.1,
                                   2.21,
                                   2.28,
                                   2.39,
                                   2.63,
                                   3.58,
                                   4.21,
                                   4.23,
                                   4.94,
                                   5.69,
                                   5.79,
                                   0.99,
                                   1.07,
                                   1.31,
                                   1.34,
                                   1.37,
                                   1.52,
                                   1.55,
                                   1.61,
                                   1.74,
                                   1.78,
                                   2.04,
                                   2.12,
                                   2.71,
                                   3.22,
                                   4.08,
                                   4.3,
                                   5.09,
                                   5.38,
                                   6.19,
                                   0.84,
                                   0.92,
                                   0.95,
                                   1.15,
                                   1.18,
                                   1.23,
                                   1.64,
                                   1.69,
                                   2.74,
                                   3.48,
                                   0.85,
                                   0.85,
                                   0.93,
                                   1.04,
                                   1.22,
                                   1.49,
                                   1.62,
                                   1.67,
                                   2.74,
                                   3.67,
                                   0.85,
                                   0.89,
                                   0.95,
                                   1.07,
                                   1.07,
                                   1.22,
                                   1.66,
                                   1.66,
                                   2.72,
                                   3.68,
                                   0.85,
                                   0.86,
                                   0.89,
                                   1.16,
                                   1.18,
                                   1.24,
                                   1.66,
                                   1.7,
                                   2.66,
                                   3.61,
                                   0.92,
                                   0.94,
                                   1.08,
                                   1.24,
                                   1.25,
                                   1.3,
                                   1.33,
                                   1.33,
                                   1.36,
                                   1.45,
                                   1.48,
                                   1.49,
                                   1.55,
                                   1.69,
                                   1.79,
                                   1.83,
                                   2.26,
                                   2.28,
                                   2.4,
                                   3.23,
                                   3.3,
                                   3.69,
                                   4.02,
                                   5.19,
                                   5.24,
                                   0.92,
                                   0.94,
                                   1.08,
                                   1.25,
                                   1.28,
                                   1.33,
                                   1.33,
                                   1.36,
                                   1.38,
                                   1.46,
                                   1.48,
                                   1.49,
                                   1.54,
                                   1.64,
                                   1.68,
                                   1.78,
                                   1.83,
                                   2.26,
                                   2.66,
                                   3.23,
                                   3.55,
                                   3.64,
                                   4.04,
                                   5.14,
                                   5.17,
                                   0.98,
                                   1.07,
                                   1.09,
                                   1.11,
                                   1.46,
                                   1.59,
                                   1.63,
                                   1.69,
                                   1.73,
                                   1.84,
                                   1.94,
                                   1.95,
                                   2.17,
                                   2.18,
                                   2.25,
                                   2.52,
                                   2.63,
                                   2.82,
                                   3.69,
                                   6.17,
                                   7.2,
                                   1.14,
                                   1.41,
                                   1.41,
                                   1.49,
                                   1.49,
                                   1.72,
                                   1.72,
                                   3.1,
                                   3.32,
                                   3.71,
                                   3.84,
                                   4.02,
                                   4.97,
                                   6.36,
                                   6.43,
                                   1.12,
                                   1.21,
                                   1.41,
                                   1.41,
                                   1.41,
                                   1.59,
                                   1.71,
                                   2.9,
                                   3.55,
                                   3.71,
                                   3.83,
                                   4.12,
                                   4.85,
                                   6.36,
                                   6.42,
                                   0.79,
                                   1.19,
                                   1.19,
                                   1.39,
                                   1.39,
                                   1.71,
                                   1.71,
                                   2.21,
                                   2.51,
                                   3.22,
                                   3.57,
                                   4.22,
                                   4.68,
                                   5.06,
                                   5.61])

    plinear_regression_data_calc = np.array([2.01,
                                    1.62,
                                    1.83,
                                    2.27,
                                    3.1,
                                    4.95,
                                    4.74,
                                    6.89,
                                    7.03,
                                    7.02,
                                    7.4,
                                    7.49,
                                    7.64,
                                    1.06,
                                    2.59,
                                    2.74,
                                    2.75,
                                    2.75,
                                    2.93,
                                    3.72,
                                    4.79,
                                    6.88,
                                    6.9,
                                    7.4,
                                    7.44,
                                    8.18,
                                    8.45,
                                    1.38,
                                    1.15,
                                    1.18,
                                    2.5,
                                    2.5,
                                    3.39,
                                    3.22,
                                    3.08,
                                    4.22,
                                    4.5,
                                    5.55,
                                    6.91,
                                    7.24,
                                    7.69,
                                    1.32,
                                    1.51,
                                    1.59,
                                    1.87,
                                    1.93,
                                    2.66,
                                    2.73,
                                    2.84,
                                    2.85,
                                    2.88,
                                    2.96,
                                    3.12,
                                    3.29,
                                    3.31,
                                    4.23,
                                    4.69,
                                    6.57,
                                    7.16,
                                    7.38,
                                    7.44,
                                    7.45,
                                    7.48,
                                    7.5,
                                    7.5,
                                    7.57,
                                    7.76,
                                    2.83,
                                    3.48,
                                    3.64,
                                    3.76,
                                    3.82,
                                    4.15,
                                    6,
                                    6.02,
                                    6.06,
                                    6.77,
                                    6.84,
                                    7.21,
                                    7.27,
                                    7.4,
                                    7.45,
                                    7.74,
                                    0.81,
                                    1.18,
                                    1.25,
                                    1.13,
                                    1.22,
                                    1.29,
                                    1.41,
                                    1.51,
                                    1.49,
                                    1.69,
                                    1.54,
                                    1.66,
                                    1.78,
                                    1.86,
                                    1.83,
                                    2.06,
                                    2.16,
                                    2.37,
                                    2.27,
                                    2.44,
                                    2.48,
                                    3.84,
                                    5.86,
                                    2.27,
                                    2.79,
                                    3.06,
                                    3.06,
                                    3.23,
                                    3.76,
                                    5.2,
                                    7.17,
                                    7.26,
                                    7.45,
                                    1.04,
                                    1.25,
                                    1.36,
                                    1.37,
                                    1.43,
                                    1.44,
                                    1.72,
                                    1.74,
                                    4.41,
                                    4.59,
                                    4.34,
                                    4.16,
                                    4.18,
                                    6.05,
                                    7.69,
                                    0.93,
                                    1,
                                    1.31,
                                    1.53,
                                    1.54,
                                    1.56,
                                    1.58,
                                    1.86,
                                    2.05,
                                    2.68,
                                    2.77,
                                    3.38,
                                    4.15,
                                    4.18,
                                    4.2,
                                    4.83,
                                    7.19,
                                    0.98,
                                    1.11,
                                    2.35,
                                    2.34,
                                    2.84,
                                    2.69,
                                    4.1,
                                    4.16,
                                    4.53,
                                    4.64,
                                    4.57,
                                    4.97,
                                    5.21,
                                    5.87,
                                    6.45,
                                    0.97,
                                    0.98,
                                    2.33,
                                    2.48,
                                    2.77,
                                    3.09,
                                    4.92,
                                    4.14,
                                    4.29,
                                    4.59,
                                    4.65,
                                    4.66,
                                    5.21,
                                    5.95,
                                    6.34,
                                    0.95,
                                    1.21,
                                    1.39,
                                    1.46,
                                    1.5,
                                    1.53,
                                    1.7,
                                    1.79,
                                    2.01,
                                    2.15,
                                    2.54,
                                    2.58,
                                    3.51,
                                    3.55,
                                    4.28,
                                    0.91,
                                    1.19,
                                    1.36,
                                    1.5,
                                    1.51,
                                    1.72,
                                    1.82,
                                    1.87,
                                    1.88,
                                    2.14,
                                    2.47,
                                    2.56,
                                    3.5,
                                    3.78,
                                    3.83,
                                    0.99,
                                    1.34,
                                    1.4,
                                    1.62,
                                    1.69,
                                    1.44,
                                    1.91,
                                    1.95,
                                    2.11,
                                    2.01,
                                    2.26,
                                    2.3,
                                    2.49,
                                    3.29,
                                    3.78,
                                    3.89,
                                    6.17,
                                    6.81,
                                    7.05,
                                    7.65,
                                    7.46,
                                    1.2,
                                    1.36,
                                    1.77,
                                    1.81,
                                    1.82,
                                    1.82,
                                    1.94,
                                    1.99,
                                    2.45,
                                    2.72,
                                    3.37,
                                    3.54,
                                    6.81,
                                    6.94,
                                    1.21,
                                    1.36,
                                    1.59,
                                    1.62,
                                    1.63,
                                    1.63,
                                    1.77,
                                    1.82,
                                    2.78,
                                    2.86,
                                    3.44,
                                    3.66,
                                    6.53,
                                    7.06,
                                    1.09,
                                    1.61,
                                    2.01,
                                    2.65,
                                    4.06,
                                    4.15,
                                    5.38,
                                    7.41,
                                    7.56,
                                    7.59,
                                    1.07,
                                    1.56,
                                    1.89,
                                    2.33,
                                    4.08,
                                    4.22,
                                    5.64,
                                    7.51,
                                    1.11,
                                    1.47,
                                    1.53,
                                    1.74,
                                    1.84,
                                    1.86,
                                    2.48,
                                    4.18,
                                    4.22,
                                    5.69,
                                    7.52,
                                    1.1,
                                    1.61,
                                    1.64,
                                    1.68,
                                    1.82,
                                    2.3,
                                    2.7,
                                    4.02,
                                    4.15,
                                    5.33,
                                    7.4,
                                    7.56,
                                    7.6,
                                    0.97,
                                    0.97,
                                    1.11,
                                    1.16,
                                    1.42,
                                    1.49,
                                    1.63,
                                    1.69,
                                    1.71,
                                    1.84,
                                    1.88,
                                    1.95,
                                    2.11,
                                    2.63,
                                    3.26,
                                    5.78,
                                    2.37,
                                    2.42,
                                    2.58,
                                    2.22,
                                    3.52,
                                    4.14,
                                    4.17,
                                    6.1,
                                    7.72,
                                    8.1,
                                    1.84,
                                    2.06,
                                    2.75,
                                    3.63,
                                    4.08,
                                    4.25,
                                    4.44,
                                    5.58,
                                    7.02,
                                    2.4,
                                    2.73,
                                    2.83,
                                    3.94,
                                    4.1,
                                    4.69,
                                    5.27,
                                    5.51,
                                    5.7,
                                    7.53,
                                    2.43,
                                    3.54,
                                    3.75,
                                    4,
                                    4.07,
                                    4.08,
                                    4.16,
                                    4.25,
                                    4.26,
                                    4.46,
                                    6.86,
                                    6.97,
                                    7.16,
                                    7.34,
                                    7.54,
                                    7.6,
                                    7.68,
                                    0.86,
                                    1.05,
                                    0.98,
                                    1.04,
                                    1.08,
                                    1.17,
                                    1.55,
                                    1.57,
                                    2.02,
                                    1.65,
                                    1.67,
                                    1.71,
                                    1.9,
                                    2.11,
                                    2.24,
                                    1.69,
                                    1.48,
                                    2.15,
                                    2.24,
                                    2.38,
                                    3.15,
                                    3.12,
                                    3.75,
                                    4.53,
                                    3.61,
                                    3.96,
                                    4.31,
                                    4.65,
                                    7.2,
                                    7.57,
                                    7.77,
                                    9.01,
                                    0.73,
                                    1.73,
                                    1.81,
                                    1.87,
                                    2.09,
                                    2.31,
                                    2.38,
                                    3,
                                    3.02,
                                    3.04,
                                    3.17,
                                    3.18,
                                    3.28,
                                    3.29,
                                    3.52,
                                    3.62,
                                    3.64,
                                    3.88,
                                    3.97,
                                    4.04,
                                    4.11,
                                    7.14,
                                    7.22,
                                    7.26,
                                    7.4,
                                    7.41,
                                    7.48,
                                    7.52,
                                    7.55,
                                    7.56,
                                    7.57,
                                    7.62,
                                    7.68,
                                    7.96,
                                    0.92,
                                    1.11,
                                    1.27,
                                    1.35,
                                    1.38,
                                    1.45,
                                    1.47,
                                    1.58,
                                    1.59,
                                    1.61,
                                    1.63,
                                    1.72,
                                    1.91,
                                    1.99,
                                    2.06,
                                    2.27,
                                    2.67,
                                    3.12,
                                    3.13,
                                    3.19,
                                    3.28,
                                    7.38,
                                    7.99,
                                    8.16,
                                    0.91,
                                    1.11,
                                    1.76,
                                    2.1,
                                    2.12,
                                    2.2,
                                    2.42,
                                    2.67,
                                    2.76,
                                    2.8,
                                    2.86,
                                    3.66,
                                    3.69,
                                    3.75,
                                    3.81,
                                    3.92,
                                    4.01,
                                    4.09,
                                    4.69,
                                    5.65,
                                    6.88,
                                    7.36,
                                    7.53,
                                    7.64,
                                    7.69,
                                    0.87,
                                    0.99,
                                    1.3,
                                    1.51,
                                    1.51,
                                    1.53,
                                    1.53,
                                    1.53,
                                    1.54,
                                    1.56,
                                    1.59,
                                    1.59,
                                    1.8,
                                    1.89,
                                    1.94,
                                    1.96,
                                    1.97,
                                    1.97,
                                    2,
                                    2.04,
                                    2.05,
                                    2.26,
                                    2.79,
                                    2.86,
                                    2.95,
                                    1.26,
                                    1.41,
                                    1.45,
                                    1.68,
                                    1.79,
                                    1.92,
                                    2,
                                    1.79,
                                    2.41,
                                    2.28,
                                    2.47,
                                    2.27,
                                    2.72,
                                    3.56,
                                    4.3,
                                    4.36,
                                    4.92,
                                    6.12,
                                    6.15,
                                    1.35,
                                    1.24,
                                    1.7,
                                    1.53,
                                    1.69,
                                    1.97,
                                    1.81,
                                    1.93,
                                    2.06,
                                    2.15,
                                    2.34,
                                    2.37,
                                    2.59,
                                    3.46,
                                    4.38,
                                    4.57,
                                    4.85,
                                    5.74,
                                    6.72,
                                    0.83,
                                    0.85,
                                    0.99,
                                    1.28,
                                    1.44,
                                    1.68,
                                    1.77,
                                    1.81,
                                    2.84,
                                    3.39,
                                    0.77,
                                    0.89,
                                    1,
                                    1.18,
                                    1.31,
                                    1.8,
                                    1.84,
                                    1.87,
                                    2.73,
                                    3.78,
                                    0.9,
                                    0.97,
                                    1,
                                    1.04,
                                    1.04,
                                    1.3,
                                    1.79,
                                    1.83,
                                    2.72,
                                    3.76,
                                    0.87,
                                    0.95,
                                    1.01,
                                    1.19,
                                    1.22,
                                    1.37,
                                    1.77,
                                    1.78,
                                    2.82,
                                    3.5,
                                    0.98,
                                    1.04,
                                    1.22,
                                    1.63,
                                    1.27,
                                    1.59,
                                    1.36,
                                    1.41,
                                    1.49,
                                    1.65,
                                    1.54,
                                    1.78,
                                    1.59,
                                    1.64,
                                    1.77,
                                    1.87,
                                    2.31,
                                    2.46,
                                    2.52,
                                    3.14,
                                    3.29,
                                    4,
                                    4.19,
                                    4.74,
                                    5.18,
                                    1.04,
                                    1.06,
                                    1.4,
                                    1.35,
                                    1.3,
                                    1.46,
                                    1.51,
                                    1.61,
                                    1.56,
                                    1.55,
                                    1.67,
                                    1.57,
                                    1.51,
                                    1.66,
                                    2.13,
                                    1.72,
                                    2.2,
                                    2.48,
                                    2.54,
                                    3.22,
                                    3.8,
                                    3.64,
                                    4.36,
                                    4.95,
                                    4.73,
                                    1.07,
                                    1.15,
                                    1.31,
                                    1.14,
                                    1.68,
                                    1.58,
                                    1.63,
                                    2.1,
                                    1.78,
                                    2.03,
                                    1.99,
                                    2.36,
                                    2.2,
                                    2.24,
                                    2.56,
                                    2.66,
                                    2.76,
                                    2.88,
                                    3.67,
                                    6.25,
                                    7.26,
                                    1.18,
                                    1.83,
                                    1.51,
                                    1.56,
                                    0.97,
                                    2.27,
                                    1.83,
                                    3.15,
                                    3.63,
                                    3.3,
                                    4.2,
                                    3.68,
                                    4.75,
                                    6.42,
                                    6.17,
                                    1.1,
                                    1.2,
                                    1.74,
                                    1.36,
                                    1.62,
                                    1.66,
                                    1.31,
                                    3.15,
                                    3.53,
                                    3.56,
                                    4.83,
                                    3.89,
                                    4.9,
                                    6.3,
                                    6.14,
                                    0.96,
                                    1.38,
                                    1.37,
                                    1.63,
                                    1.55,
                                    1.91,
                                    2.02,
                                    2.51,
                                    2.61,
                                    3.84,
                                    4.16,
                                    4.73,
                                    4.57,
                                    5.64,
                                    6.17])

    scaled = plinear_regression_data_calc * slope +intercept

    errors = plinear_regression_data_exp - scaled

    mu, std = norm.fit(errors)

    return mu,std

def remove_labile_protons(sdffile,lbls,shifts):

    f = sdffile + '.sdf'

    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()
    obconversion.ReadFile(obmol, f)

    CI = []

    for atom in OBMolAtomIter(obmol):

        if atom.GetAtomicNum() == 1:

            for NbrAtom in OBAtomAtomIter(atom):

                if (NbrAtom.GetAtomicNum() == 8):
                    CI.append( 'H' + str(atom.GetIndex() + 1))

    #remove these carbons

    for C in CI:

        ind = lbls.index(C)

        lbls.remove(C)

        for l in shifts:
            l.pop(ind)

    return lbls,shifts

def removecrossassignments(exp,calc,labels):

    #sort these in decending order

    s  = np.argsort(calc)[::-1]

    calc = calc[s]

    exp = exp[s]

    labels = labels[s]

    #generate difference matrix
    switch = 0

    expcopy = np.array(exp)

    while switch == 0:

        swapm = np.zeros([len(calc), len(calc)])

        for i,Hi in enumerate(expcopy):
            for j,Hj in enumerate(expcopy):

                if i>j:

                    swapm[i, j] = 0
                else:
                    swapm[i,j] = round(Hi - Hj,1)

        w = np.argwhere(swapm < 0)

        if len(w > 0):
            expcopy[w[0]] = expcopy[w[0][::-1]]

        else:
            switch =1

    return calc, expcopy,labels

def methyl_protons(file):

    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()
    obconversion.ReadFile(obmol, file)

    methyl_protons = []

    for atom in OBMolAtomIter(obmol):

        count = 0

        nbrprotons = []

        for NbrAtom in OBAtomAtomIter(atom):

            if (atom.GetAtomicNum() == 6) & (NbrAtom.GetAtomicNum() == 1):

                l = NbrAtom.GetIndex()

                count += 1

                nbrprotons.append('H' + str(l + 1))

        if count == 3:
            methyl_protons.append(nbrprotons)

    return methyl_protons
