import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt
from scipy.optimize import linear_sum_assignment as optimise
from scipy.stats import linregress
import copy
from openbabel import *
from scipy.stats import gaussian_kde


def AssignCarbon(NMRData, Isomers, settings):
    for isomer in Isomers:

        assigned_shifts, assigned_peaks, assigned_labels, scaled_shifts = iterative_assignment(
            NMRData.carbondata["exppeaks"],
            NMRData.carbondata["xdata"],
            NMRData.carbondata["ydata"],
            isomer.Cshifts, isomer.Clabels)

        # add to isomers instance

        isomer.Cexp = [''] * len(isomer.Cshifts)

        for label, peak in zip(assigned_labels, assigned_peaks):
            w = isomer.Clabels.index(label)

            isomer.Cexp[w] = peak

    return Isomers


def iterative_assignment(picked_peaks, spectral_xdata_ppm, total_spectral_ydata, calculated_shifts, C_labels):
    calculated_shifts = np.array(calculated_shifts)

    original_C_labels = np.array(C_labels)

    s = np.argsort(np.array(calculated_shifts))

    calculated_shifts = calculated_shifts[s]

    scaled_shifts = copy.copy(calculated_shifts)

    C_labels = original_C_labels[s]

    exp_peaks = spectral_xdata_ppm[picked_peaks]

    new_assigned_peaks = []

    new_assigned_shifts = []

    for lnum in range(0, 2):

        if lnum == 0:

            scaled_shifts = external_scale_carbon_shifts(calculated_shifts)

            scaled_mu = 0

            scaled_std = 2.486068603518297

            copy_calc_shifts = copy.copy(calculated_shifts)

        elif lnum == 1:

            old_assigned_shifts = copy.copy(new_assigned_shifts)

            old_assigned_peaks = copy.copy(new_assigned_peaks)

            scaled_shifts, slope, intercept = internal_scale_carbon_shifts(old_assigned_shifts, old_assigned_peaks,
                                                                           calculated_shifts)
            scaled_mu = 0

            scaled_std = 10

            copy_calc_shifts = copy.copy(calculated_shifts)

        ####calculate difference matrix

        diff_matrix = np.zeros((len(calculated_shifts), len(exp_peaks)))

        for ind1, i in enumerate(scaled_shifts):
            for ind2, j in enumerate(exp_peaks):
                diff_matrix[ind1, ind2] = j - i

        ####find any errors larger than 10 ppm and nans
        ####calculate pos matirx

        pos_matrix = carbon_probabilities(diff_matrix, scaled_mu, scaled_std)

        pos_matrix[abs(diff_matrix) >= 10] = 0

        pos_matrix[np.isnan(pos_matrix)] = 0

        ####calculate amp matrix

        amp_matrix = amp_kde(total_spectral_ydata, picked_peaks, pos_matrix, calculated_shifts)

        ####duplicate the pos matrix along the horizontal to allow multiple assignment weighting

        pos_matrixc = copy.copy(pos_matrix)

        for d in range(0, len(calculated_shifts) - 1):
            pos_matrix = np.hstack((pos_matrix, pos_matrixc))

        ####calculate the probability matrix

        prob_matrix = (pos_matrix * amp_matrix) ** 0.5

        ####check for any shifts that have zero probabilites for all peaks

        b = np.where(np.sum(prob_matrix, 1) == 0)

        prob_matrix = np.delete(prob_matrix, b, 0)

        unassignable_shifts = calculated_shifts[b]

        copy_calc_shifts = np.delete(copy_calc_shifts, b)

        copy_labels = np.delete(C_labels, b)

        ####do the assignment

        vertical_ind, horizontal_ind = optimise(1 - prob_matrix)

        horizontal_ind = horizontal_ind % len(picked_peaks)

        opt_peaks = exp_peaks[horizontal_ind]

        opt_shifts = copy_calc_shifts[vertical_ind]

        opt_labels = copy_labels[vertical_ind]

        ####do some sorting

        so = np.argsort(opt_shifts)

        new_assigned_peaks = opt_peaks[so]

        new_assigned_shifts = opt_shifts[so]

        new_assigned_labels = opt_labels[so]

    ############################
    # in the third round only reassign shifts that have had a change of bias

    old_assigned_shifts = copy.copy(new_assigned_shifts)

    old_assigned_peaks = copy.copy(new_assigned_peaks)

    new_assigned_shifts = copy.copy(new_assigned_shifts)

    new_assigned_peaks = copy.copy(new_assigned_peaks)

    new_assigned_labels = copy.copy(new_assigned_labels)

    bias_weights = []

    # find unassigned peaks

    ampdivide = np.zeros(len(picked_peaks))

    peak_amps = total_spectral_ydata[picked_peaks]

    reassign_shifts_ind = []

    for i in old_assigned_peaks:
        w = np.where(exp_peaks == i)

        ampdivide[w] += 1

    c = 0

    for shift, peak in zip(old_assigned_shifts, old_assigned_peaks):

        # find where peaks are within 20 ppm window

        w = np.where((exp_peaks < peak + 10) & (exp_peaks > peak - 10))[0]

        if len(w) > 0:

            # find maximum peak height within this window - when taking into account how many times the peak has already been assigned

            # find amplitude of peak given how many times it has been assigned

            assigned_amp = (peak_amps[exp_peaks == peak] / ampdivide[exp_peaks == peak])[0]

            # find amplitude of max peak in the 20 ppm window given how many times it would be assigned if the current shift was assigned to it as well

            div_amps = peak_amps / (ampdivide + 1)

            pi = np.where(exp_peaks == peak)

            div_amps[pi] = peak_amps[pi] / ampdivide[pi]

            max_window_amp = np.max(div_amps[w])

            ratio = max_window_amp / assigned_amp

            if ratio > 1:
                bias_weights.append(ratio)

                reassign_shifts_ind.append(c)

        c += 1

    ####reassign the shifts with a bias above zero in order of bias to peak within ten ppm with largest unassigned amplitude

    bias_weights = np.array(bias_weights)

    reassign_shifts = np.array(old_assigned_shifts)[reassign_shifts_ind]

    s = np.argsort(bias_weights)

    reassign_shifts = reassign_shifts[s]

    reassign_shifts_ind = np.array(reassign_shifts_ind)[s]

    for shift, ind in zip(reassign_shifts, reassign_shifts_ind):

        # find peak this shift is assigned to

        p = new_assigned_peaks[ind]

        pi = np.where(exp_peaks == p)

        new_peak_amps = peak_amps / (ampdivide + 1)

        new_peak_amps[pi] = peak_amps[pi] / (ampdivide[pi])

        # find peaks within 10 ppm

        w = np.where((exp_peaks < p + 10) & (
                exp_peaks > p - 10))[0]

        if len(w) > 0:
            assigned_peak = exp_peaks[w[np.argmax(new_peak_amps[w])]]

            new_assigned_peaks[ind] = assigned_peak

        # recalculate estimated peak heights

        ampdivide = np.zeros(len(picked_peaks))

        for i in new_assigned_peaks:
            w = np.where(exp_peaks == i)

            ampdivide[w] += 1

    #############################

    # remove cross assignments

    new_assigned_shifts, new_assigned_peaks, new_assigned_labels = removecrossassignments(new_assigned_peaks,
                                                                                          new_assigned_shifts,
                                                                                          new_assigned_labels)

    #### sortoutput wrt original H labels

    assigned_labels = []

    assigned_shifts = []

    assigned_peaks = []

    for label in original_C_labels:

        wh = np.where(new_assigned_labels == label)[0]

        assigned_labels.append(label)

        if len(wh) > 0:

            assigned_shifts.append(new_assigned_shifts[wh[0]])

            assigned_peaks.append(new_assigned_peaks[wh[0]])

        else:
            assigned_shifts.append('')

            assigned_peaks.append('')

    '''

    assigned_shifts = np.array(new_assigned_shifts)[indv]

    assigned_peaks = np.array(new_assigned_peaks)[indv]

    assigned_labels = new_assigned_labels[indv]

    '''

    return assigned_shifts, assigned_peaks, assigned_labels, scaled_shifts


def external_scale_carbon_shifts(calculated_shifts):
    scaled = calculated_shifts * 0.9601578792266342 - 1.2625604390657088

    return scaled


def internal_scale_carbon_shifts(assigned_shifts, assigned_peaks, calculated_shifts):
    slope, intercept, r_value, p_value, std_err = linregress(assigned_shifts, assigned_peaks)

    scaled_shifts = calculated_shifts * slope + intercept

    # plt.close()

    # plt.plot(assigned_shifts,assigned_shifts,color= 'grey')

    # plt.plot(assigned_shifts, assigned_peaks,'ro')

    # plt.plot(np.array(assigned_shifts) *slope + intercept,assigned_peaks,'co')

    # plt.plot(assigned_shifts,np.array(assigned_shifts)*slope + intercept)

    # plt.show()

    return scaled_shifts, slope, intercept


'''
def amp_weighting(total_spectral_ydata, picked_peaks, prob_matrix,shifts,steep_weights):

    print(np.round(steep_weights,2))

    peak_amps = total_spectral_ydata[picked_peaks]

    peak_amps = peak_amps / np.max(peak_amps)

    samps = np.sort(peak_amps)

    thresh = samps[max([len(samps) - len(shifts),0])]

    #plt.plot(total_spectral_ydata)
    #plt.plot(picked_peaks,total_spectral_ydata[picked_peaks],'co')
    #plt.show()

    def weight_curve(x, thresh, steep):

        y = 1 / (1 + np.exp(-steep * (x - thresh)))

        return y

    steep =100* np.std(peak_amps)

    #x = np.linspace(0,1,1000)

    #plt.plot(x,weight_curve(x,thresh,steep))

    amp_matrix = np.zeros((len(shifts),len(picked_peaks)))

    for  i in range(0,len(amp_matrix[:,0])):
        amp_matrix[i,:] = weight_curve(peak_amps, thresh, steep * steep_weights[i])

    #plt.plot(peak_amps,weight_curve(peak_amps,thresh,steep),'co')
    #plt.show()

    return amp_matrix
'''


def amp_weighting(total_spectral_ydata, picked_peaks, prob_matrix, shifts, steep_weights):
    peak_amps = total_spectral_ydata[picked_peaks]

    peak_amps = peak_amps

    peak_amps = peak_amps / np.max(peak_amps)

    duplicated_amps = copy.copy(peak_amps)

    for d in range(0, len(shifts) - 1):
        duplicated_amps = np.hstack((duplicated_amps, peak_amps * (0.25) ** (d + 1)))

    samps = np.sort(peak_amps)

    thresh = samps[max([len(samps) - len(shifts), 0])]

    # plt.plot(total_spectral_ydata)
    # plt.plot(picked_peaks,total_spectral_ydata[picked_peaks],'co')
    # plt.show()

    def weight_curve(x, thresh, steep):

        y = 1 / (1 + np.exp(-steep * (x - thresh)))

        return y

    steep = np.std(peak_amps)

    x = np.linspace(0, 1, 1000)

    plt.plot(x, weight_curve(x, thresh, 100 * steep))
    plt.plot(x, weight_curve(x, thresh, 400 * steep))

    amp_matrix = np.zeros((len(shifts), len(duplicated_amps)))

    for i in range(0, len(amp_matrix[:, 0])):
        amp_matrix[i, :] = weight_curve(duplicated_amps, thresh, steep * 100 ** steep_weights[i])

    plt.plot(peak_amps, weight_curve(peak_amps, thresh, 100 * steep), 'o', color='C0')
    plt.plot(peak_amps, weight_curve(peak_amps, thresh, 400 * steep), 'co', color='C1')
    plt.show()

    return amp_matrix


def amp_kde(total_spectral_ydata, picked_peaks, prob_matrix, shifts):
    peak_amps = total_spectral_ydata[picked_peaks]

    peak_amps = peak_amps

    peak_amps = peak_amps / np.max(peak_amps)

    x = np.linspace(0, 1, 1000)

    kde = gaussian_kde(peak_amps)

    y = kde.evaluate(x)

    # plt.plot(x,y)
    # plt.show()

    # find minima

    # find maxima in second derivative,

    ddy = np.diff(y, 2)

    ddy1 = np.roll(ddy, 1)

    ddyn1 = np.roll(ddy, -1)

    w = np.where((ddy[1:-1] > ddy1[1:-1]) & (ddy[1:-1] > ddyn1[1:-1]))[0] + 2

    # add zero and one values

    if w[0] != 0:
        w = np.hstack((0, w))

    if w[-1] != len(ddy) - 2:
        w = np.hstack((w, len(y) - 1))

    minima = x[w]

    i = 0

    groups = np.zeros(len(peak_amps))

    number_in_group = []

    for m, m1 in zip(minima[:-1], minima[1:]):
        w = np.where((peak_amps > m) & (peak_amps <= m1))[0]

        groups[w] = i

        number_in_group.append(len(w))

        i += 1

    groups = groups.astype(int)

    # convert group numbers to weights

    cumsum = np.cumsum(number_in_group[::-1])[::-1]

    weight_values = len(shifts) / cumsum

    weight_values /= np.max(weight_values)

    peak_weights = weight_values[groups]

    duplicated_weights = copy.copy(peak_weights)

    # do multiple assignment weights

    for d in range(0, len(shifts) - 1):
        duplicated_weights = np.hstack(
            (duplicated_weights, peak_weights * (0.125 ** (np.max(groups) - groups + 1)) ** (d + 1)))

    duplicated_weightsc = copy.copy(duplicated_weights)

    # duplicate along vertical

    for d in range(0, len(shifts) - 1):
        duplicated_weights = np.vstack((duplicated_weights, duplicated_weightsc))

    # renormalise

    for i in range(duplicated_weights.shape[0]):
        duplicated_weights[i, :] = duplicated_weights[i, :] / np.sum(duplicated_weights[i, :])

    return duplicated_weights


def multiple_assignment_weighting(prob_matrix):
    # shifts are columns
    # peaks are rows
    # duplicate matrix along the horizontal but multiply by 0.5 each time

    pmcopy = copy.copy(prob_matrix)

    for i, shift in enumerate(pmcopy[:, 0]):
        prob_matrix = np.hstack((prob_matrix, pmcopy * (1 / (i + 1))))

    return prob_matrix


def carbon_probabilities(diff_matrix, scaled_mu, scaled_std):
    # unscaled error data
    # prob_matrix = norm.pdf(diff_matrix,-4.585221745350501,3.2910538642479277) / norm.pdf(0,-4.585221745350501,3.2910538642479277)

    # scaled error data
    # prob_matrix = norm.pdf(diff_matrix,0,2.486068603518297) / norm.pdf(0,0,2.486068603518297)

    prob_matrix = norm.pdf(diff_matrix, scaled_mu, scaled_std) / norm.pdf(scaled_mu, scaled_mu, scaled_std)

    # prob_matrix = norm.pdf(diff_matrix, 0, 10) / norm.pdf(scaled_mu, 0, 10)

    # normalise rows

    for i in range(prob_matrix.shape[0]):
        prob_matrix[i, :] = prob_matrix[i, :] / np.sum(prob_matrix[i, :])

    return prob_matrix


def simulate_spectrum(spectral_xdata_ppm, calc_shifts):
    y = np.zeros(len(spectral_xdata_ppm))

    for shift in calc_shifts:
        y += lorentzian(spectral_xdata_ppm, 0.001, shift, 0.2)

    return y


def simulate_spectrum(spectral_xdata_ppm, calc_shifts, assigned_peaks, set_exp):
    for ind, shift in enumerate(calc_shifts):
        exp_p = assigned_peaks[ind]

        ind2 = set_exp.index(exp_p)

        y = lorentzian(spectral_xdata_ppm, 0.001, shift, 0.2)

        plt.plot(spectral_xdata_ppm, y + 1.05, color='C' + str(ind2))


def simulate_calc_data(spectral_xdata_ppm, calculated_locations, simulated_ydata):
    ###simulate calcutated data

    simulated_calc_ydata = np.zeros(len(spectral_xdata_ppm))

    for peak in calculated_locations:
        y = np.exp(-0.5 * ((spectral_xdata_ppm - peak) / 0.002) ** 2)
        simulated_calc_ydata += y

    scaling_factor = np.amax(simulated_ydata) / np.amax(simulated_calc_ydata)

    simulated_calc_ydata = simulated_calc_ydata * scaling_factor

    return simulated_calc_ydata


def lorentzian(p, w, p0, A):
    x = (p0 - p) / (w / 2)
    L = A / (1 + x ** 2)

    return L


def remove_iodine(sdffile, lbls, shifts):
    f = sdffile + '.sdf'

    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()
    obconversion.ReadFile(obmol, f)

    CI = []

    for atom in OBMolAtomIter(obmol):

        if atom.GetAtomicNum() == 6:

            for NbrAtom in OBAtomAtomIter(atom):

                if (NbrAtom.GetAtomicNum() == 53):
                    CI.append('C' + str(atom.GetIndex() + 1))

    # remove these carbons

    for C in CI:

        ind = lbls.index(C)

        lbls.remove(C)

        for l in shifts:
            l.pop(ind)

    return lbls, shifts


def removecrossassignments(exp, calc, labels):
    # sort these in decending order

    s = np.argsort(calc)[::-1]

    calc = calc[s]

    exp = exp[s]

    labels = labels[s]

    # generate difference matrix
    switch = 0

    expcopy = np.array(exp)

    while switch == 0:

        swapm = np.zeros([len(calc), len(calc)])

        for i, Hi in enumerate(expcopy):
            for j, Hj in enumerate(expcopy):

                if i > j:

                    swapm[i, j] = 0
                else:
                    swapm[i, j] = round(Hi - Hj, 1)

        w = np.argwhere(swapm < 0)

        if len(w > 0):
            expcopy[w[0]] = expcopy[w[0][::-1]]

        else:
            switch = 1

    return calc, expcopy, labels