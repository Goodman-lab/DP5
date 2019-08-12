import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt
from scipy.optimize import linear_sum_assignment as optimise
from scipy.stats import linregress
import copy
from openbabel import *
from scipy.stats import gaussian_kde

def AssignCarbon(NMRData,Isomers,settings):

    for isomer in Isomers:

        assigned_shifts, assigned_peaks, assigned_labels, scaled_shifts = iterative_assignment(NMRData.carbondata["exppeaks"],
                                                                                               NMRData.carbondata["xdata"],
                                                                                               NMRData.carbondata["ydata"],
                                                                                               isomer.Cshifts, isomer.Clabels)

    #add to isomers instance

        isomer.Cexp = [''] * len(isomer.Cshifts)

        for label,peak in zip(assigned_labels,assigned_peaks):

            w = isomer.Clabels.index(label)

            isomer.Cexp[w] = peak

    return Isomers

def iterative_assignment(picked_peaks, spectral_xdata_ppm, total_spectral_ydata, calculated_shifts, C_labels):

    calculated_shifts = np.array(calculated_shifts)

    original_C_labels = np.array(C_labels)

    s = np.argsort(np.array(calculated_shifts))

    calculated_shifts = calculated_shifts[s]

    C_labels = original_C_labels[s]

    # print("Clabels = ",C_labels)
    # print("calc shifts = " ,calculated_shifts)

    exp_peaks = spectral_xdata_ppm[picked_peaks]

    lnum = 0

    new_assigned_peaks = []
    old_assigned_peaks = [1]

    steep_weights = np.ones(len(calculated_shifts))

    while lnum < 2:

        diff_matrix = np.zeros((len(calculated_shifts), len(exp_peaks)))

        if lnum == 0:

            old_assigned_peaks = copy.copy(new_assigned_peaks)

            scaled_shifts = external_scale_carbon_shifts(calculated_shifts)

            scaled_mu = 0

            scaled_std = 2.486068603518297

        elif lnum == 1:

            old_assigned_shifts = copy.copy(new_assigned_shifts)
            old_assigned_peaks = copy.copy(new_assigned_peaks)

            scaled_shifts, slope, intercept = internal_scale_carbon_shifts(old_assigned_shifts, old_assigned_peaks,
                                                                           calculated_shifts)

            # scaled_mu,scaled_std = scale_params(slope,intercept)

            scaled_mu = 0

            scaled_std = 10

        else:
            old_assigned_shifts = copy.copy(new_assigned_shifts)
            old_assigned_peaks = copy.copy(new_assigned_peaks)

        for ind1, i in enumerate(scaled_shifts):
            for ind2, j in enumerate(exp_peaks):
                diff_matrix[ind1, ind2] = j - i

        pos_matrix = carbon_probabilities(diff_matrix, scaled_mu, scaled_std)

        # change nans to zeros

        pos_matrix[np.isnan(pos_matrix)] = 0

        # amp_matrix = amp_weighting(total_spectral_ydata, picked_peaks, pos_matrix,calculated_shifts,steep_weights)

        amp_matrix = amp_kde(total_spectral_ydata, picked_peaks, pos_matrix, calculated_shifts, steep_weights)

        # print("exp_peaks",np.round(exp_peaks,2))

        # apply multiple assignment penalty to amp weighting and duplicate to pos matrix along the horizontal

        pos_matrixc = copy.copy(pos_matrix)

        for d in range(0, len(calculated_shifts) - 1):
            pos_matrix = np.hstack((pos_matrix, pos_matrixc))

        prob_matrix = (pos_matrix * amp_matrix) ** 0.5

        # do the assignment

        prob_matrix[np.isnan(prob_matrix)] = 0

        ####################################################

        # for each shift work out if assignment is based on position or amplitude

        # check if a row is all zeros

        vertical_ind, horizontal_ind = optimise(1 - prob_matrix)

        for v, h in zip(vertical_ind, horizontal_ind):

            # find the peaks that have not been assigned

            # find the shift corresponding to the vertical index

            s = exp_peaks[h % len(picked_peaks)]

            # find all unassigned shifts within 20 ppm

            unassigned_ind = np.array(list(set(np.arange(0, np.size(prob_matrix, 1) - 1)) - set(horizontal_ind)))

            w = np.where((exp_peaks[unassigned_ind % len(picked_peaks)] < s + 20) & (
                        exp_peaks[unassigned_ind % len(picked_peaks)] > s - 20))[0]

            if len(w) > 0:
                # print("shift ",s,"assigned score " ,amp_matrix[v, h] / np.sum(amp_matrix[v, :]))

                window_unassigned = unassigned_ind[w]

                # print("unassigned ", np.round(exp_peaks[window_unassigned % len(picked_peaks)]))

                # print("amp scores ", amp_matrix[v,window_unassigned] / amp_matrix[v, h])

                # bayesamp = amp_matrix[v, h] / np.sum(amp_matrix[v, :])

                assigned_amp = amp_matrix[v, h]

                max_window_amp = np.max(amp_matrix[v, window_unassigned])

                ratio = (max_window_amp / assigned_amp)

                steep_weights[v] = max(1, ratio)

        ####################################################

        # matrix has been duplicated along the horizontal add this line to get original indicies

        horizontal_ind = horizontal_ind % len(picked_peaks)

        opt_peaks = exp_peaks[horizontal_ind]

        opt_shifts = calculated_shifts[vertical_ind]

        opt_labels = C_labels[vertical_ind]

        # sort

        so = np.argsort(opt_shifts)

        opt_peaks = opt_peaks[so]

        opt_shifts = opt_shifts[so]

        opt_labels = opt_labels[so]

        # now must check - if two shifts are assigned to the same peak - v likely one of them will have had their amplitude steepness increased
        # but now must make sure this increase is performed for all the shifts assigned to the same peak, this prevents the next assignement just reversing the order
        # i.e if shifts 1 and 2 assigned to a peak the scores will be [f(x), f(0.25x)] the steepness scores might look like [ 0 , 4]
        # hence in the next round of assignement, instead of assigning a peak else where the system will prioritise peak 2 and give the scores [f(0.25x), f(x)]
        # the steepness scores wil change to [4,0] but this is not helpful, the cycle will just continue.

        for s, p in zip(opt_shifts, set(opt_peaks)):
            # find where peak have been assigned multiple times

            wh = np.where(opt_peaks == p)

            # change all of the weights to this value

            steep_weights[wh] = np.max(steep_weights[wh])

        if lnum == 1:
            opt_shifts, opt_peaks, opt_labels = removecrossassignments(opt_peaks, opt_shifts, opt_labels)

        #### sort output wrt original H labels

        indv = []

        for label in original_C_labels:
            wh = np.where(opt_labels == label)

            indv.append(wh[0][0])

        assigned_shifts = opt_shifts[indv]

        assigned_peaks = opt_peaks[indv]

        assigned_labels = opt_labels[indv]

        so = np.argsort(assigned_shifts)

        # print("assigned shifts", np.round(assigned_shifts[so],2))
        # print("assigned labels" ,assigned_labels[so])
        # print("assigned peaks", np.round(assigned_peaks[so],2))
        # print("steep weights",np.round(steep_weights,2))

        # ind = np.argsort(assigned_shifts)

        # assigned_shifts = assigned_shifts[ind].tolist()
        # assigned_peaks = assigned_peaks[ind].tolist()
        # assigned_labels = assigned_labels[ind].tolist()

        new_assigned_shifts = copy.copy(assigned_shifts)
        new_assigned_peaks = copy.copy(assigned_peaks)

        lnum += 1

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

    print(peak_amps)

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


def amp_kde(total_spectral_ydata, picked_peaks, prob_matrix, shifts, steep_weights):
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

    # maxima = x[np.r_[True, y[1:] > y[:-1]] & np.r_[y[:-1] > y[1:], True]]

    # add zero and one values

    if w[0] != 0:
        w = np.hstack((0, w))

    if w[-1] != len(ddy) - 2:
        w = np.hstack((w, len(y) - 1))

    minima = x[w]

    i = 0

    groups = np.zeros(len(peak_amps))

    number_in_group = []

    ##plt.plot(x,y)
    ##plt.plot(x[:-2],ddy/np.max(ddy))

    ##plt.plot(minima,y[w],'co')

    ##plt.show()

    for m, m1 in zip(minima[:-1], minima[1:]):
        w = np.where((peak_amps > m) & (peak_amps <= m1))[0]

        groups[w] = i

        number_in_group.append(len(w))

        i += 1

    groups = groups.astype(int)

    # print("number in group",number_in_group)

    # print("groups",groups)

    # convert group numbers to weights

    cumsum = np.cumsum(number_in_group[::-1])[::-1]

    weight_values = len(shifts) / cumsum

    weight_values /= np.max(weight_values)

    peak_weights = weight_values[groups]

    # print("peak weights",np.round( peak_weights,2))

    duplicated_weights = copy.copy(peak_weights)

    # do multiple assignment weights

    for d in range(0, len(shifts) - 1):
        duplicated_weights = np.hstack(
            (duplicated_weights, peak_weights * (0.125 ** (np.max(groups) - groups + 1)) ** (d + 1)))
        # duplicated_weights = np.hstack(
        #   (duplicated_weights, peak_weights * (0.125) ** (d + 1)))

    duplicated_weightsc = copy.copy(duplicated_weights)

    # duplicate along vertical

    for d in range(0, len(shifts) - 1):
        duplicated_weights = np.vstack((duplicated_weights, duplicated_weightsc))

    # print(duplicated_weights)

    duplicated_weights = duplicated_weights ** np.transpose([steep_weights])

    # print(duplicated_weights)

    # renormalise

    for i in range(duplicated_weights.shape[0]):
        duplicated_weights[i, :] = duplicated_weights[i, :] / np.sum(duplicated_weights[i, :])

    return duplicated_weights


def multiple_assignment_weighting(prob_matrix):
    # shifts are columns
    # peaks are rows
    # duplicate matrix along the horizontal but multiply by 0.5 each time

    pmcopy = copy.copy(prob_matrix)

    for shift in pmcopy[:, 0]:
        prob_matrix = np.hstack((prob_matrix, pmcopy * 0.25))

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


def scale_params(slope, intercept):
    clinear_regression_data_exp = np.array([24.5,
                                            36.36,
                                            59.4,
                                            59.58,
                                            71.05,
                                            114.64,
                                            115.71,
                                            115.79,
                                            118.22,
                                            127.5,
                                            127.52,
                                            127.88,
                                            133.98,
                                            142.14,
                                            157.42,
                                            158,
                                            161.02,
                                            167.33,
                                            19.4,
                                            42,
                                            53.6,
                                            54.7,
                                            55.3,
                                            71.3,
                                            114,
                                            117.8,
                                            119.5,
                                            123.7,
                                            125.7,
                                            130.2,
                                            130.6,
                                            133.4,
                                            147.7,
                                            158.3,
                                            160.2,
                                            22.6,
                                            33.7,
                                            34.3,
                                            42.8,
                                            44.7,
                                            45.4,
                                            69.1,
                                            72.1,
                                            116.6,
                                            123.6,
                                            124.4,
                                            134,
                                            136.6,
                                            142.3,
                                            159.3,
                                            163.6,
                                            165.2,
                                            176.9,
                                            182.4,
                                            20.1,
                                            24.8,
                                            25.8,
                                            28.6,
                                            38.9,
                                            46.7,
                                            47.7,
                                            56.1,
                                            58.3,
                                            72.4,
                                            126.3,
                                            127.2,
                                            127.5,
                                            128.2,
                                            128.4,
                                            128.6,
                                            129,
                                            135.1,
                                            135.7,
                                            143,
                                            155.4,
                                            23.9,
                                            33.6,
                                            52.1,
                                            56.2,
                                            56.7,
                                            101.2,
                                            106.6,
                                            107.5,
                                            108.2,
                                            111.2,
                                            118.6,
                                            120.1,
                                            120.7,
                                            122.5,
                                            126.2,
                                            132.8,
                                            135.4,
                                            136.6,
                                            147.1,
                                            147.9,
                                            166.4,
                                            166.8,
                                            10.92,
                                            17.14,
                                            20.39,
                                            23.08,
                                            29.83,
                                            31.31,
                                            32.59,
                                            33.62,
                                            35.34,
                                            35.38,
                                            36.18,
                                            38.43,
                                            42.53,
                                            50.21,
                                            53.68,
                                            80.93,
                                            123.41,
                                            171.8,
                                            199.7,
                                            28.2,
                                            47.1,
                                            53.6,
                                            55.5,
                                            63.6,
                                            64.4,
                                            76.2,
                                            125.7,
                                            126.1,
                                            126.7,
                                            146.6,
                                            169.7,
                                            14.3,
                                            18.7,
                                            22.2,
                                            27.9,
                                            28.3,
                                            65.9,
                                            73.9,
                                            74.6,
                                            79.7,
                                            91.1,
                                            129.6,
                                            138,
                                            152.2,
                                            152.2,
                                            152.2,
                                            9.3,
                                            9.5,
                                            14.2,
                                            23.7,
                                            25.8,
                                            26.3,
                                            33.6,
                                            49.2,
                                            59,
                                            60.8,
                                            74.8,
                                            81.7,
                                            129.6,
                                            137.5,
                                            166.3,
                                            170.9,
                                            17.9,
                                            20.4,
                                            38.3,
                                            41.3,
                                            46.4,
                                            50.6,
                                            55.2,
                                            69.2,
                                            72.6,
                                            73.4,
                                            75.9,
                                            77,
                                            115.8,
                                            127.3,
                                            133.8,
                                            161.3,
                                            165,
                                            172.5,
                                            17.4,
                                            20.7,
                                            30.7,
                                            37.3,
                                            47.3,
                                            49.2,
                                            53.8,
                                            67.6,
                                            71.5,
                                            73.1,
                                            73.7,
                                            75.8,
                                            118.2,
                                            124.8,
                                            139,
                                            160.8,
                                            164.7,
                                            170.5,
                                            18.3,
                                            21.4,
                                            25.2,
                                            29.8,
                                            39.5,
                                            41.3,
                                            42.75,
                                            42.8,
                                            53.4,
                                            65.6,
                                            76.9,
                                            83.4,
                                            208.2,
                                            18.6,
                                            21.1,
                                            25.9,
                                            30,
                                            39.2,
                                            39.5,
                                            42.3,
                                            42.9,
                                            53.1,
                                            66.3,
                                            77,
                                            84.7,
                                            208,
                                            16,
                                            20.2,
                                            20.2,
                                            22,
                                            22.4,
                                            24.4,
                                            31.9,
                                            31.9,
                                            34.8,
                                            40.7,
                                            42.6,
                                            71.6,
                                            101.1,
                                            112.2,
                                            112.8,
                                            113.1,
                                            126.2,
                                            128.6,
                                            129.5,
                                            143.5,
                                            155.9,
                                            161.3,
                                            162.8,
                                            179.5,
                                            18.4,
                                            19.7,
                                            20.1,
                                            20.1,
                                            24.6,
                                            29.9,
                                            32.9,
                                            36.1,
                                            47.9,
                                            57.9,
                                            92.5,
                                            125.3,
                                            125.5,
                                            133.8,
                                            138.9,
                                            140.7,
                                            151.1,
                                            179,
                                            183,
                                            185.8,
                                            18.9,
                                            19.3,
                                            20.1,
                                            20.1,
                                            25,
                                            27.5,
                                            27.8,
                                            38.8,
                                            46.9,
                                            54.7,
                                            82.1,
                                            124.1,
                                            125.5,
                                            133,
                                            133.7,
                                            137.8,
                                            150.8,
                                            179,
                                            184,
                                            185.8,
                                            14.1,
                                            19.1,
                                            35.3,
                                            43.4,
                                            76.8,
                                            78.3,
                                            85.6,
                                            125.6,
                                            127.3,
                                            128.4,
                                            143.3,
                                            14.1,
                                            19.1,
                                            36.5,
                                            44,
                                            76.8,
                                            79.5,
                                            87.1,
                                            125.9,
                                            127.4,
                                            128.3,
                                            142.1,
                                            14.3,
                                            19.7,
                                            31.3,
                                            44.7,
                                            73.7,
                                            78.2,
                                            83.1,
                                            125.3,
                                            127.2,
                                            128.4,
                                            143.8,
                                            14.3,
                                            19.6,
                                            30.9,
                                            44.2,
                                            73.1,
                                            78.7,
                                            83.6,
                                            125.9,
                                            127.3,
                                            128.4,
                                            143,
                                            13.2,
                                            20.5,
                                            24,
                                            25.5,
                                            25.8,
                                            33.5,
                                            34.2,
                                            36.5,
                                            38.2,
                                            45.6,
                                            50.7,
                                            80.1,
                                            94.3,
                                            106,
                                            172.7,
                                            25.52,
                                            32.21,
                                            62.71,
                                            82.11,
                                            84.52,
                                            124.38,
                                            138.25,
                                            145.68,
                                            147.66,
                                            156.68,
                                            10.24,
                                            34.33,
                                            58.27,
                                            58.9,
                                            81.56,
                                            82.11,
                                            107.59,
                                            134.1,
                                            148.48,
                                            161.78,
                                            39.2,
                                            54.1,
                                            55.2,
                                            63.1,
                                            70.4,
                                            109.3,
                                            116.2,
                                            136,
                                            151.3,
                                            151.5,
                                            153.5,
                                            156.9,
                                            18.81,
                                            33.44,
                                            61.43,
                                            70.42,
                                            74.67,
                                            78.48,
                                            81.21,
                                            81.32,
                                            115.88,
                                            123.39,
                                            126.25,
                                            126.36,
                                            126.95,
                                            129.06,
                                            129.65,
                                            130.52,
                                            134.93,
                                            137.36,
                                            138.24,
                                            140.22,
                                            143.63,
                                            161.37,
                                            16.1,
                                            18.4,
                                            19.7,
                                            21.7,
                                            22.1,
                                            23.1,
                                            23.3,
                                            24.4,
                                            24.5,
                                            25.6,
                                            27.3,
                                            28.1,
                                            29.9,
                                            36.2,
                                            37.77,
                                            37.82,
                                            47.8,
                                            50.9,
                                            52.5,
                                            59.2,
                                            61.5,
                                            65.1,
                                            122.7,
                                            123.3,
                                            123.9,
                                            127.1,
                                            131.6,
                                            137.1,
                                            168.8,
                                            169.2,
                                            170.1,
                                            170.6,
                                            174.1,
                                            174.3,
                                            16.1,
                                            18.7,
                                            19.6,
                                            21.7,
                                            25.4,
                                            38.2,
                                            40,
                                            40.8,
                                            41.1,
                                            41.8,
                                            48.7,
                                            54.4,
                                            63.1,
                                            69.7,
                                            70.2,
                                            124.6,
                                            126,
                                            126.3,
                                            128.2,
                                            128.4,
                                            129,
                                            129.2,
                                            129.3,
                                            130.3,
                                            138,
                                            138.1,
                                            154.2,
                                            156.5,
                                            168.8,
                                            170.7,
                                            10.4,
                                            18.6,
                                            18.9,
                                            22.1,
                                            24.3,
                                            28.2,
                                            29.3,
                                            32.7,
                                            34.7,
                                            35.4,
                                            40.8,
                                            46.9,
                                            50.3,
                                            52,
                                            55.2,
                                            60.8,
                                            61.5,
                                            64.1,
                                            68.7,
                                            84,
                                            105,
                                            121.8,
                                            127,
                                            133.3,
                                            134.6,
                                            142.6,
                                            185.8,
                                            194.7,
                                            19.9,
                                            20.2,
                                            25.8,
                                            27.3,
                                            35.7,
                                            45.4,
                                            53.7,
                                            55.1,
                                            58.9,
                                            69.6,
                                            70.8,
                                            72.8,
                                            73.4,
                                            109.3,
                                            114.1,
                                            126,
                                            126.5,
                                            126.5,
                                            128.5,
                                            129.4,
                                            137.7,
                                            150.8,
                                            155.4,
                                            21.1,
                                            22.2,
                                            23.2,
                                            26.5,
                                            32.7,
                                            34.6,
                                            34.8,
                                            37.6,
                                            37.6,
                                            40.1,
                                            40.1,
                                            41.4,
                                            42,
                                            43.6,
                                            56,
                                            58.7,
                                            65.2,
                                            18.5,
                                            21.7,
                                            21.9,
                                            23.6,
                                            31,
                                            32.1,
                                            40.5,
                                            66.7,
                                            71.3,
                                            71.5,
                                            74,
                                            132.1,
                                            137,
                                            170.1,
                                            19.2,
                                            22.6,
                                            25.3,
                                            27.8,
                                            30.7,
                                            32.1,
                                            39.9,
                                            67.2,
                                            69.6,
                                            69.8,
                                            71.5,
                                            129,
                                            138.2,
                                            169.8,
                                            14.6,
                                            16.6,
                                            21.3,
                                            24.3,
                                            25.2,
                                            33.4,
                                            39.4,
                                            42.3,
                                            78.2,
                                            181.4,
                                            9.9,
                                            15.7,
                                            21.3,
                                            24.3,
                                            25.3,
                                            33.5,
                                            41.6,
                                            41.8,
                                            76.3,
                                            181.9,
                                            11.1,
                                            14.5,
                                            21.6,
                                            23.6,
                                            25,
                                            33.3,
                                            42.1,
                                            42.6,
                                            75.6,
                                            181.6,
                                            12.5,
                                            14.2,
                                            22.2,
                                            23.1,
                                            25,
                                            32.1,
                                            43.1,
                                            43.2,
                                            75.6,
                                            181.4,
                                            14.2,
                                            19.7,
                                            25.2,
                                            35.6,
                                            36.2,
                                            37.1,
                                            38.8,
                                            41.3,
                                            43.2,
                                            45.6,
                                            46,
                                            56.6,
                                            69.1,
                                            72.8,
                                            75.7,
                                            76.5,
                                            84.5,
                                            174.7,
                                            14.1,
                                            20,
                                            26,
                                            32.6,
                                            36.2,
                                            37.4,
                                            37.9,
                                            41,
                                            43.2,
                                            43.5,
                                            45.2,
                                            56.4,
                                            69.2,
                                            71.3,
                                            73.9,
                                            77,
                                            77.1,
                                            173,
                                            15.1,
                                            17.2,
                                            19.2,
                                            23.3,
                                            24.2,
                                            30,
                                            31.3,
                                            35.5,
                                            35.5,
                                            36.7,
                                            43.7,
                                            47.9,
                                            48.6,
                                            50.3,
                                            52,
                                            109.6,
                                            122.7,
                                            140.6,
                                            148.1,
                                            178.4,
                                            181,
                                            19.3,
                                            19.4,
                                            32,
                                            35.1,
                                            40,
                                            53.9,
                                            67.8,
                                            71.5,
                                            102.8,
                                            112.5,
                                            120.6,
                                            137.6,
                                            159.3,
                                            160.7,
                                            170.7,
                                            204.4,
                                            21.5,
                                            22.7,
                                            32.5,
                                            35.6,
                                            39.5,
                                            54.5,
                                            67.1,
                                            73.5,
                                            102.7,
                                            112.1,
                                            121.2,
                                            137,
                                            158.5,
                                            160.4,
                                            170.9,
                                            204.9,
                                            13.2,
                                            18.1,
                                            39,
                                            41.1,
                                            45.2,
                                            71.6,
                                            72.5,
                                            76,
                                            78.1,
                                            128.1,
                                            136.4,
                                            169.8])

    clinear_regression_data_calc = np.array([28.25,
                                             42.92,
                                             64.52,
                                             66.83,
                                             73.75,
                                             119.92,
                                             117.39,
                                             121.44,
                                             123.62,
                                             133.35,
                                             134.76,
                                             140.62,
                                             138.19,
                                             150.68,
                                             160.63,
                                             165.72,
                                             167.88,
                                             172.62,
                                             21.15,
                                             46.68,
                                             54.97,
                                             60.72,
                                             61.22,
                                             75.95,
                                             118.25,
                                             118.47,
                                             121.54,
                                             127.32,
                                             128.31,
                                             135.8,
                                             137.73,
                                             143.05,
                                             147.87,
                                             161.9,
                                             163.21,
                                             22.99,
                                             37.71,
                                             34.52,
                                             49.42,
                                             42.87,
                                             41.08,
                                             74.44,
                                             75.45,
                                             120.78,
                                             131.39,
                                             128.14,
                                             138.96,
                                             144.05,
                                             151.8,
                                             167.93,
                                             173.7,
                                             169.37,
                                             184.5,
                                             184.76,
                                             23.19,
                                             27.6,
                                             29.5,
                                             32.98,
                                             38.82,
                                             48.89,
                                             49.29,
                                             60.03,
                                             60.22,
                                             75.11,
                                             130.98,
                                             132.06,
                                             132.07,
                                             132.98,
                                             133.85,
                                             134.31,
                                             135.27,
                                             138.35,
                                             143.19,
                                             151.23,
                                             162.04,
                                             26.39,
                                             34.55,
                                             55.82,
                                             57.73,
                                             61.32,
                                             107.4,
                                             109.9,
                                             113.97,
                                             114.71,
                                             116.31,
                                             123.46,
                                             124.45,
                                             126.12,
                                             126.35,
                                             130.98,
                                             138.37,
                                             139.12,
                                             144.09,
                                             149.42,
                                             149.75,
                                             179.68,
                                             181.15,
                                             12.94,
                                             19.35,
                                             24.44,
                                             26.57,
                                             33.4,
                                             35.17,
                                             36.94,
                                             35.95,
                                             39.84,
                                             40.07,
                                             39.52,
                                             44.56,
                                             47.89,
                                             52.68,
                                             58.06,
                                             84.8,
                                             131.23,
                                             179.91,
                                             205.61,
                                             34.33,
                                             50.45,
                                             55.11,
                                             59.25,
                                             68.63,
                                             73.22,
                                             83.78,
                                             132.42,
                                             135.27,
                                             135.79,
                                             151.68,
                                             186.57,
                                             16.97,
                                             22.59,
                                             28.21,
                                             32.54,
                                             34.01,
                                             68.87,
                                             80.17,
                                             81.55,
                                             84.17,
                                             98.84,
                                             131.05,
                                             144.94,
                                             156.61,
                                             158.51,
                                             159.98,
                                             11.03,
                                             12.23,
                                             16.24,
                                             24.25,
                                             28.57,
                                             29.24,
                                             36.72,
                                             51.87,
                                             62.23,
                                             68.2,
                                             76.61,
                                             83.19,
                                             137.39,
                                             149.72,
                                             174.77,
                                             177.37,
                                             19.44,
                                             22.99,
                                             42.64,
                                             46.21,
                                             50.93,
                                             54.83,
                                             59.75,
                                             72.75,
                                             79.96,
                                             76.77,
                                             79.97,
                                             84.3,
                                             124.95,
                                             137.75,
                                             145.14,
                                             169.85,
                                             172.75,
                                             185.3,
                                             20.65,
                                             23.08,
                                             36.36,
                                             41.88,
                                             51.69,
                                             53.34,
                                             59.11,
                                             71.31,
                                             78.59,
                                             76.54,
                                             81.38,
                                             80.09,
                                             127.21,
                                             135.49,
                                             150.25,
                                             169.81,
                                             172.37,
                                             183.85,
                                             22.59,
                                             23.18,
                                             27.03,
                                             32.23,
                                             43.2,
                                             43.52,
                                             44.32,
                                             46.54,
                                             55.89,
                                             69.35,
                                             80.09,
                                             87.17,
                                             221.39,
                                             22.74,
                                             22.95,
                                             27.35,
                                             32.3,
                                             42.99,
                                             43.75,
                                             44.23,
                                             46.02,
                                             55.32,
                                             69.29,
                                             79.74,
                                             88.94,
                                             221.49,
                                             18.08,
                                             22.2,
                                             21.05,
                                             25.83,
                                             24.62,
                                             28.28,
                                             31.73,
                                             35.84,
                                             39.94,
                                             49.33,
                                             46.86,
                                             73.81,
                                             105.61,
                                             117.61,
                                             119.72,
                                             122.88,
                                             138.45,
                                             135.4,
                                             141.75,
                                             146.67,
                                             163.13,
                                             169.71,
                                             167.48,
                                             182.75,
                                             21.14,
                                             21.61,
                                             21.79,
                                             22.94,
                                             30.44,
                                             34.29,
                                             38.14,
                                             39.94,
                                             49.31,
                                             62.7,
                                             89.13,
                                             129.17,
                                             132.2,
                                             138.79,
                                             141.65,
                                             149.79,
                                             157.24,
                                             188.93,
                                             190.54,
                                             193.58,
                                             21.15,
                                             21.31,
                                             21.35,
                                             21.95,
                                             28.88,
                                             29.65,
                                             30.47,
                                             42.3,
                                             49.33,
                                             58.65,
                                             82.27,
                                             127.78,
                                             131.6,
                                             139.29,
                                             139.35,
                                             145.2,
                                             157,
                                             189.89,
                                             190.3,
                                             192.94,
                                             16.32,
                                             23.69,
                                             40.15,
                                             46.15,
                                             79.36,
                                             81.68,
                                             88.08,
                                             129.27,
                                             130.98,
                                             133.42,
                                             152.43,
                                             16.23,
                                             23.69,
                                             41.3,
                                             45.13,
                                             79.62,
                                             83.94,
                                             89.9,
                                             129.59,
                                             130.72,
                                             133.02,
                                             152.16,
                                             16.51,
                                             24.3,
                                             35.09,
                                             45.73,
                                             77.22,
                                             83.86,
                                             86.16,
                                             129.04,
                                             130.61,
                                             133.15,
                                             152.32,
                                             16.52,
                                             23.61,
                                             36.19,
                                             45.92,
                                             76.94,
                                             82.03,
                                             86.25,
                                             129.57,
                                             131.26,
                                             133.73,
                                             153.9,
                                             16.2,
                                             21.35,
                                             27.58,
                                             27.6,
                                             29.42,
                                             36.67,
                                             38.75,
                                             40.34,
                                             41.54,
                                             49.39,
                                             53.93,
                                             85.89,
                                             101.56,
                                             110.08,
                                             184.49,
                                             29.75,
                                             34.91,
                                             70,
                                             82.95,
                                             93.75,
                                             132.37,
                                             144.67,
                                             150.2,
                                             148.5,
                                             160.95,
                                             15.76,
                                             38.65,
                                             68.01,
                                             68.72,
                                             91.9,
                                             98.9,
                                             119.82,
                                             142.91,
                                             155.03,
                                             167.26,
                                             44.13,
                                             60.46,
                                             62.2,
                                             70.31,
                                             78.74,
                                             121.01,
                                             125.15,
                                             143.79,
                                             148.85,
                                             156.05,
                                             160.86,
                                             166.48,
                                             22.28,
                                             35.95,
                                             67.84,
                                             76.56,
                                             80.18,
                                             81.35,
                                             81.76,
                                             85.74,
                                             121.51,
                                             131.3,
                                             132,
                                             132.07,
                                             135.8,
                                             136.41,
                                             137.04,
                                             138.19,
                                             143.25,
                                             144.57,
                                             147.99,
                                             148.01,
                                             152.72,
                                             169.23,
                                             18.88,
                                             19.39,
                                             21.92,
                                             20.41,
                                             22.53,
                                             25.17,
                                             24.76,
                                             27.42,
                                             28.96,
                                             29.27,
                                             26.99,
                                             23.98,
                                             36.25,
                                             39.32,
                                             39.57,
                                             39.64,
                                             54.3,
                                             57.13,
                                             46.58,
                                             67.39,
                                             58.26,
                                             69.27,
                                             126.45,
                                             127.8,
                                             124.44,
                                             135.05,
                                             139.29,
                                             144.26,
                                             185.87,
                                             179.44,
                                             186.57,
                                             178.26,
                                             188.39,
                                             183.95,
                                             19.22,
                                             20.43,
                                             21.01,
                                             23.68,
                                             29.63,
                                             39.4,
                                             40.6,
                                             41.88,
                                             43.08,
                                             43.61,
                                             55.42,
                                             63.1,
                                             68.06,
                                             74,
                                             75.75,
                                             130.71,
                                             130.72,
                                             131.43,
                                             133.67,
                                             133.94,
                                             134.42,
                                             134.58,
                                             134.98,
                                             138.55,
                                             146.85,
                                             147.78,
                                             160,
                                             160.23,
                                             179.17,
                                             180.48,
                                             11.18,
                                             22.73,
                                             23.17,
                                             23.59,
                                             23.7,
                                             29.1,
                                             30.11,
                                             35.73,
                                             36.77,
                                             40.12,
                                             43.34,
                                             50.73,
                                             52.33,
                                             56.29,
                                             57.36,
                                             60.79,
                                             63.17,
                                             65.71,
                                             75.15,
                                             92.11,
                                             120.92,
                                             128.02,
                                             136.17,
                                             139.23,
                                             141.48,
                                             147.33,
                                             193.82,
                                             205.47,
                                             21.51,
                                             22.04,
                                             28.04,
                                             31.38,
                                             39.18,
                                             47.37,
                                             56.97,
                                             62.22,
                                             62.62,
                                             71.12,
                                             73.32,
                                             77.81,
                                             78.27,
                                             116.67,
                                             119.93,
                                             131.4,
                                             131.69,
                                             133.78,
                                             133.9,
                                             135.58,
                                             142.31,
                                             158.37,
                                             161.36,
                                             23.33,
                                             24.25,
                                             25.61,
                                             28.78,
                                             34.95,
                                             36.75,
                                             37.57,
                                             38.35,
                                             40.56,
                                             40.63,
                                             42.03,
                                             42.09,
                                             42.55,
                                             45.7,
                                             58.65,
                                             59.83,
                                             65.85,
                                             19.42,
                                             24.33,
                                             23.88,
                                             27.02,
                                             34.88,
                                             34.09,
                                             42.93,
                                             69.84,
                                             70.72,
                                             71.43,
                                             77.67,
                                             143.19,
                                             146.62,
                                             182.58,
                                             21.01,
                                             26.32,
                                             29.3,
                                             30.78,
                                             34.77,
                                             35.67,
                                             43.27,
                                             71.18,
                                             71.48,
                                             73.97,
                                             73.97,
                                             138.37,
                                             148.51,
                                             183.49,
                                             15.7,
                                             17.85,
                                             20.99,
                                             25.54,
                                             28.88,
                                             39.1,
                                             42.5,
                                             43.83,
                                             82.72,
                                             186.81,
                                             11.91,
                                             14.45,
                                             20.89,
                                             25.63,
                                             28.91,
                                             36.25,
                                             41.51,
                                             44.2,
                                             79.94,
                                             188.17,
                                             12.27,
                                             15.88,
                                             20.87,
                                             25.47,
                                             29.05,
                                             36.52,
                                             41.29,
                                             42.05,
                                             79.37,
                                             187.85,
                                             15.4,
                                             16.93,
                                             21.57,
                                             24.99,
                                             28.85,
                                             38.28,
                                             43.19,
                                             43.52,
                                             80.73,
                                             186.7,
                                             15.91,
                                             22.73,
                                             27.59,
                                             39.42,
                                             40.4,
                                             41.73,
                                             40.33,
                                             44.96,
                                             45.02,
                                             47.99,
                                             48.83,
                                             57.97,
                                             76.91,
                                             75.69,
                                             76.48,
                                             78.82,
                                             87.18,
                                             185.54,
                                             16.04,
                                             23.59,
                                             27.65,
                                             35.3,
                                             39.56,
                                             40.69,
                                             40.96,
                                             42.71,
                                             45.56,
                                             45.91,
                                             47.53,
                                             57.25,
                                             76.91,
                                             75.24,
                                             77.33,
                                             78.6,
                                             78.46,
                                             185.08,
                                             18.86,
                                             19.59,
                                             23.52,
                                             27.38,
                                             27.09,
                                             32.78,
                                             35.91,
                                             37.79,
                                             40.74,
                                             38.74,
                                             47.78,
                                             53.46,
                                             53.68,
                                             53.48,
                                             53.38,
                                             115.1,
                                             129.17,
                                             146.22,
                                             155.33,
                                             191.29,
                                             184.73,
                                             22.84,
                                             25.79,
                                             38.1,
                                             40.77,
                                             39.33,
                                             49.54,
                                             74.73,
                                             74.98,
                                             104.22,
                                             116.54,
                                             126.14,
                                             144.53,
                                             161.62,
                                             164.21,
                                             183.41,
                                             224.46,
                                             22.16,
                                             24.65,
                                             36.42,
                                             36.75,
                                             43.16,
                                             50.72,
                                             75.09,
                                             70.41,
                                             105.07,
                                             116.75,
                                             124.6,
                                             148.33,
                                             165.63,
                                             164.99,
                                             182.45,
                                             219.26,
                                             16.07,
                                             21.86,
                                             42.7,
                                             42.45,
                                             45.92,
                                             75.93,
                                             75.26,
                                             81.15,
                                             81.93,
                                             138.31,
                                             147.38,
                                             182.28])

    scaled = clinear_regression_data_calc * slope + intercept

    errors = clinear_regression_data_exp - scaled

    mu, std = norm.fit(errors)

    return mu, std


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