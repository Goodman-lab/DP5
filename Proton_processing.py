from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import gmean
from lmfit import Minimizer, Parameters, report_fit
import nmrglue as ng
from scipy.stats import norm
import pickle

try:
    from openbabel.openbabel import OBConversion, OBMol, OBAtomAtomIter, OBMolAtomIter
except ImportError:
    from openbabel import *

from scipy.optimize import linear_sum_assignment as optimise
import copy
from scipy.interpolate import InterpolatedUnivariateSpline
import pathos.multiprocessing as mp
import time
import itertools
import os

def process_proton(NMR_file, settings, datatype):

    total_spectral_ydata, spectral_xdata_ppm, corr_distance, uc, noise_std, peak_regions = spectral_processing(NMR_file,
                                                                                                               datatype)

    print('peak picking')

    gradient_peaks, gradient_regions, gradient_groups, std = gradient_peak_picking(total_spectral_ydata, corr_distance,
                                                                                   uc, noise_std, peak_regions)

    import time

    start = time.time()

    picked_peaks, grouped_peaks, peak_regions, sim_y, total_params = multiproc_BIC_minimisation(gradient_regions,
                                                                                                gradient_groups,
                                                                                                total_spectral_ydata,
                                                                                                corr_distance,
                                                                                                uc, noise_std)

    end = time.time()

    print("minimisation time = " + str((end - start) / 60) + " mins")

    peak_regions, picked_peaks, grouped_peaks, spectral_xdata_ppm, solvent_region_ind = editsolvent_removal2(
        settings.Solvent, total_spectral_ydata, spectral_xdata_ppm, picked_peaks, peak_regions, grouped_peaks,
        total_params,
        uc)

    print('simulating spectrum')

    sim_regions, full_sim_data = simulate_regions(total_params, peak_regions, grouped_peaks, total_spectral_ydata,
                                                  spectral_xdata_ppm)

    peak_regions, grouped_peaks, sim_regions, integral_sum, cummulative_vectors, integrals, number_of_protons_structure, optimum_proton_number, total_integral = find_integrals(

        settings.InputFiles[0], peak_regions, grouped_peaks, sim_regions, picked_peaks, total_params,
        total_spectral_ydata, solvent_region_ind)

    # find region centres

    centres = weighted_region_centres(peak_regions, total_spectral_ydata)

    ################

    exp_peaks = []

    integrals = np.array([int(i) for i in integrals])

    for ind, peak in enumerate(centres):
        exp_peaks += [peak] * integrals[ind]

    integrals = integrals[integrals > 0.5]

    exp_peaks = spectral_xdata_ppm[exp_peaks]

    exp_peaks = np.array([ round(i,4) for i in exp_peaks])

    return exp_peaks, spectral_xdata_ppm, total_spectral_ydata, integrals, peak_regions, centres, cummulative_vectors, integral_sum, picked_peaks, total_params, sim_regions


def guess_udic(dic, data):
    """
    Guess parameters of universal dictionary from dic, data pair.
    Parameters
    ----------
    dic : dict
        Dictionary of JCAMP-DX parameters.
    data : ndarray
        Array of NMR data.
    Returns
    -------
    udic : dict
        Universal dictionary of spectral parameters.
    """

    # create an empty universal dictionary
    udic = ng.fileiobase.create_blank_udic(1)

    # update default values (currently only 1D possible)
    # "label"
    try:
        label_value = dic[".OBSERVENUCLEUS"][0].replace("^", "")
        udic[0]["label"] = label_value
    except KeyError:
        # sometimes INSTRUMENTAL PARAMETERS is used:
        try:
            label_value = dic["INSTRUMENTALPARAMETERS"][0].replace("^", "")
            udic[0]["label"] = label_value
        except KeyError:
            pass

    # "obs"
    obs_freq = None
    try:
        obs_freq = float(dic[".OBSERVEFREQUENCY"][0])
        udic[0]["obs"] = obs_freq

    except KeyError:
        pass

    # "size"
    if isinstance(data, list):
        data = data[0]  # if list [R,I]
    if data is not None:
        udic[0]["size"] = len(data)

    # "sw"
    # get firstx, lastx and unit
    firstx, lastx, isppm = ng.jcampdx._find_firstx_lastx(dic)

    # ppm data: convert to Hz
    if isppm:
        if obs_freq:
            firstx = firstx * obs_freq
            lastx = lastx * obs_freq
        else:
            firstx, lastx = (None, None)

    if firstx is not None and lastx is not None:
        udic[0]["sw"] = abs(lastx - firstx)

        # keys not found in standard&required JCAMP-DX keys and thus left default:
        # car, complex, encoding

        udic[0]['car'] = firstx - abs(lastx - firstx) / 2

    return udic


def spectral_processing(file, datatype):
    print('Processing Spectrum')

    if datatype == 'jcamp':

        dic, total_spectral_ydata = ng.jcampdx.read(file)  # read file

        total_spectral_ydata = total_spectral_ydata[0] + 1j * total_spectral_ydata[1]

        total_spectral_ydata = ng.proc_base.ifft_positive(total_spectral_ydata)

    else:

        dic, total_spectral_ydata = ng.bruker.read(file)  # read file

        total_spectral_ydata = ng.bruker.remove_digital_filter(dic, total_spectral_ydata)  # remove the digital filter

    total_spectral_ydata = ng.proc_base.zf_double(total_spectral_ydata, 4)

    total_spectral_ydata = ng.proc_base.fft_positive(total_spectral_ydata)  # Fourier transform

    corr_distance = estimate_autocorrelation(total_spectral_ydata)

    # normalise the data

    m = max(np.max(abs(np.real(total_spectral_ydata))), np.max(abs(np.imag(total_spectral_ydata))))

    total_spectral_ydata = np.real(total_spectral_ydata / m) + 1j * np.imag(total_spectral_ydata / m)

    if datatype == 'jcamp':

        udic = guess_udic(dic, total_spectral_ydata)

    else:

        udic = ng.bruker.guess_udic(dic, total_spectral_ydata)  # sorting units

    uc = ng.fileiobase.uc_from_udic(udic)  # unit conversion element

    spectral_xdata_ppm = uc.ppm_scale()  # ppmscale creation

    # baseline and phasing

    tydata = ACMEWLRhybrid(total_spectral_ydata, corr_distance)

    # find final noise distribution
    classification, sigma = baseline_find_signal(tydata, corr_distance, True, 1)

    # fall back phasing if fit doesnt converge

    # calculate negative area

    # draw regions

    peak_regions = []

    c1 = np.roll(classification, 1)

    diff = classification - c1

    s_start = np.where(diff == 1)[0]

    s_end = np.where(diff == -1)[0] - 1

    for r in range(len(s_start)):
        peak_regions.append(np.arange(s_start[r], s_end[r]))

    tydata = tydata / np.max(abs(tydata))

    print("corrected")

    return tydata, spectral_xdata_ppm, corr_distance, uc, sigma, peak_regions


def estimate_autocorrelation(total_spectral_ydata):
    # note this region may have a baseline distortion

    y = np.real(total_spectral_ydata[0:10000])

    params = Parameters()

    # define a basleine polnomial

    order = 6

    for p in range(order + 1):
        params.add('p' + str(p), value=0)

    def poly(params, order, y):

        bl = np.zeros(len(y))
        x = np.arange(len(y))

        for p in range(order + 1):
            bl += params['p' + str(p)] * x ** (p)

        return bl

    def res(params, order, y):

        bl = poly(params, order, y)

        r = abs(y - bl)

        return r

    out = Minimizer(res, params,
                    fcn_args=(order, y))

    results = out.minimize()

    bl = poly(results.params, order, y)

    y = y - bl

    t0 = np.sum(y * y)

    c = 1

    tc = 1

    t = []

    while tc > 0.36:
        tc = np.sum(np.roll(y, c) * y) / t0

        t.append(tc)

        c += 1

    print('autocorrelation distance = ' + str(c))

    return c


def acme(y, corr_distance):
    params = Parameters()

    phase_order = 3

    for p in range(phase_order + 1):
        params.add('p' + str(p), value=0, min=-np.pi, max=np.pi)

    def acmescore(params, im, real, phase_order):

        """
        Phase correction using ACME algorithm by Chen Li et al.
        Journal of Magnetic Resonance 158 (2002) 164-168

        Parameters
        ----------
        pd : tuple
            Current p0 and p1 values
        data : ndarray
            Array of NMR data.

        Returns
        -------
        score : float
            Value of the objective function (phase score)

        """

        data = ps(params, im, real, phase_order)

        ##########

        # calculate entropy of non corrected data, calculate penalty for baseline corrected data
        #  - keep as vector to use the default lmfit method

        # Calculation of first derivatives of signal regions

        ds1 = np.abs((data[1:] - data[:-1]))

        p1 = ds1 / np.sum(ds1)

        # Calculation of entropy
        p1[p1 == 0] = 1

        h1 = -p1 * np.log(p1)
        # h1s = np.sum(h1)

        # Calculation of penalty
        pfun = 0.0
        as_ = data - np.abs(data)
        # as_ = databl - np.abs(databl)
        sumas = np.sum(as_)

        if sumas < 0:
            # pfun = pfun + np.sum((as_ / 2) ** 2)
            pfun = (as_[1:] / 2) ** 2

        p = 1000 * pfun

        return h1 + p

    out = Minimizer(acmescore, params,
                    fcn_args=(np.imag(y), np.real(y), phase_order))

    results = out.minimize()

    p = results.params

    p.pretty_print()

    y = ps(p, np.imag(y), np.real(y), phase_order)

    classification, sigma = baseline_find_signal(y, corr_distance, True, 1)
    r = gen_baseline(np.real(y), classification, corr_distance)
    y -= r

    return y


def ACMEWLRhybrid(y, corr_distance):

    def residual_function(params, im, real):

        # phase the region

        data = ps(params, im, real, 0)

        # make new baseline for this region

        r = np.linspace(data[0], data[-1], len(real))

        # find negative area

        data -= r

        ds1 = np.abs((data[1:] - data[:-1]))

        p1 = ds1 / np.sum(ds1)

        # Calculation of entropy
        p1[p1 == 0] = 1

        h1 = -p1 * np.log(p1)
        h1s = np.sum(h1)

        # Calculation of penalty
        pfun = 0.0

        as_ = data - np.abs(data)

        sumas = np.sum(as_)

        if sumas < 0:
            pfun = (as_[1:] / 2) ** 2

        p = np.sum(pfun)

        return h1s + 1000 * p

    # find regions

    classification, sigma = baseline_find_signal(y, corr_distance, True, 1)

    c1 = np.roll(classification, 1)

    diff = classification - c1

    s_start = np.where(diff == 1)[0]

    s_end = np.where(diff == -1)[0] - 1

    peak_regions = []

    for r in range(len(s_start)):
        peak_regions.append(np.arange(s_start[r], s_end[r]))

    # for region in peak_regions:
    #    plt.plot(region,y[region],color = 'C1')

    # phase each region independently

    phase_angles = []

    weights = []

    centres = []

    for region in peak_regions:
        params = Parameters()

        params.add('p0', value=0, min=-np.pi, max=np.pi)

        out = Minimizer(residual_function, params,
                        fcn_args=(np.imag(y[region]), np.real(y[region])))

        results = out.minimize('brute')

        p = results.params

        phase_angles.append(p['p0'] * 1)

        # find weight

        data = ps(p, np.imag(y[region]), np.real(y[region]), 0)

        # make new baseline for this region

        r = np.linspace(data[0], data[-1], len(data))

        # find negative area

        res = data - r

        weights.append(abs(np.sum(res[res > 0] / np.sum(y[y > 0]))))

        centres.append(np.median(region) / len(y))

    sw = sum(weights)

    weights = [w / sw for w in weights]

    # do weighted linear regression on the regions

    # do outlier analysis

    switch = 0

    centres = np.array(centres)

    weights = np.array(weights)

    sweights = np.argsort(weights)[::-1]

    phase_angles = np.array(phase_angles)

    ind1 = 0

    while switch == 0:

        intercept, gradient = np.polynomial.polynomial.polyfit(centres, phase_angles, deg=1, w=weights)

        predicted_angles = gradient * centres + intercept

        weighted_res = np.abs(predicted_angles - phase_angles) * weights

        # find where largest weighted residual is

        max_res = sweights[ind1]

        s = 0

        if phase_angles[max_res] > 0:

            s = -1

            phase_angles[max_res] -= 2 * np.pi

        else:

            s = +1

            phase_angles[max_res] += 2 * np.pi

        intercept1, gradient1 = np.polynomial.polynomial.polyfit(centres, phase_angles, deg=1, w=weights)

        new_predicted_angles = gradient1 * centres + intercept1

        new_weighted_res = np.abs(new_predicted_angles - phase_angles) * weights

        if np.sum(new_weighted_res)  > np.sum(weighted_res):

            switch = 1

            phase_angles[max_res] += -2*np.pi*s

        ind1 +=1

    # phase the data

    p_final = Parameters()

    p_final.add('p0', value=intercept)
    p_final.add('p1', value=gradient)

    # p_final.pretty_print()

    y = ps(p_final, np.imag(y), np.real(y), 1)


    classification, sigma = baseline_find_signal(y, corr_distance, True, 1)
    r = gen_baseline(np.real(y), classification, corr_distance)
    y -= r

    return np.real(y)


def ps(param, im, real, phase_order):
    x = np.linspace(0, 1, len(real))

    angle = np.zeros(len(x))

    for p in range(phase_order + 1):
        angle += param['p' + str(p)] * x ** (p)

    # phase the data

    R = real * np.cos(angle) - im * np.sin(angle)

    return R


def baseline_find_signal(y_data, cdist, dev, t):
    wd = int(cdist) * 10

    sd_all = _get_sd(y_data, wd)

    snvectort = np.zeros(len(y_data))

    sv = []

    for i in range(0, 4 * cdist):
        x = np.arange(i + wd, len(y_data) - wd, 4 * cdist)

        sample = y_data[x]

        sd_set = _get_sd(sample, wd)

        s = _find_noise_sd(sd_set, 0.999)

        sv.append(s)

    print("done")

    sigma = np.mean(sv)

    b = np.linspace(-0.001, 0.001, 1000)

    if dev == True:

        w = np.where(sd_all > t * sigma)[0]

    else:

        w = np.where(y_data > t * sigma)[0]

    snvectort[w] = 1

    sn_vector = np.zeros(len(y_data))

    w = cdist

    for i in np.arange(len(sn_vector)):
        if snvectort[i] == 1:
            sn_vector[np.maximum(0, i - w):np.minimum(i + w, len(sn_vector))] = 1

    return sn_vector, sigma


def gen_baseline(y_data, sn_vector, corr_distance):
    points = np.arange(len(y_data))

    spl = InterpolatedUnivariateSpline(points[sn_vector == 0], y_data[sn_vector == 0], k=1)

    r = spl(points)

    # is corr distance odd or even

    if corr_distance % 2 == 0:
        kernel = np.ones((corr_distance + 1) * 10) / ((corr_distance + 1) * 10)
    else:
        kernel = np.ones((corr_distance) * 10) / ((corr_distance) * 10)

    r = np.convolve(r, kernel, mode='same')

    return r


def _rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def _get_sd(data, k):
    return np.std(_rolling_window(data, k), -1)


def _find_noise_sd(sd_set, ratio):
    '''Calculate the median m1 from SDset. exclude the elements greater
    than 2m1from SDset and recalculate the median m2. Repeat until
    m2/m1 converge(sd_set)'''
    m1 = np.median(sd_set)
    S = sd_set <= 2.0 * m1
    tmp = S * sd_set
    sd_set = tmp[tmp != 0]
    m2 = np.median(sd_set)
    while m2 / m1 < ratio:
        m1 = np.median(sd_set)
        S = sd_set <= 2.0 * m1
        tmp = S * sd_set
        sd_set = tmp[tmp != 0]

    return m2


########################################################################################################################
# peak picking and minimisation
########################################################################################################################


def p7(x, mu, std, v, A):
    x1 = (mu - x) / (std / 2)

    y = (1 - v) * (1 / (1 + (x1) ** 2)) + (v) * ((1 + ((x1) ** 2) / 2) / (1 + (x1) ** 2 + (x1) ** 4))

    y *= A

    return y


def p7residual(params, x, picked_points, y_data, region, differential):
    y = np.zeros(len(x))

    for peak in picked_points:
        y += p7(x, params['mu' + str(peak)], params['std' + str(peak)], params['vregion' + str(region)],
                params['A' + str(peak)])

    if differential == True:

        res = abs(y - y_data)

        av = np.average(res)

        av2 = np.average(res ** 2)

        difference = (res ** 2 - av2) - (res - av) ** 2

    else:
        difference = (y - y_data) ** 2

    return difference


def p7residualsolvent(params, x, picked_points, y_data, region, differential):
    y = np.zeros(len(x))

    for peak in picked_points:
        y += p7(x, params['mu' + str(peak)], params['std' + str(peak)], params['vregion' + str(region)],
                params['A' + str(peak)] * params['global_amp'])

    if differential == True:

        dy = np.gradient(y)
        ddy = np.gradient(y)

        dy_ydata = np.gradient(y_data)
        ddy_ydata = np.gradient(dy_ydata)

        difference = (y - y_data) ** 2 + (dy_ydata - dy) ** 2 + (ddy_ydata - ddy) ** 2

    else:
        difference = abs(y - y_data)

    return difference


def p7simsolvent(params, x, picked_points, region):
    y = np.zeros(len(x))

    for peak in picked_points:
        y += p7(x, params['mu' + str(peak)], params['std' + str(peak)], params['vregion' + str(region)],
                params['A' + str(peak)] * params['global_amp'])

    return y


def p7sim(params, x, picked_points, region):
    y = np.zeros(len(x))

    for peak in picked_points:
        y += p7(x, params['mu' + str(peak)], params['std' + str(peak)], params['vregion' + str(region)],
                params['A' + str(peak)])

    return y


def p7plot(params, region, group, ind, xppm):
    region_j = np.zeros(len(region))

    for peak in group:
        j = p7(region, params['mu' + str(peak)], params['std' + str(peak)], params['vregion' + str(ind)],
               params['A' + str(peak)])

        region_j += j

    return region_j


############

###########

###########


def gradient_peak_picking(y_data, corr_distance, uc, std, binary_map_regions):
    print("     gradient peak picking")

    final_peaks = []

    # estimate std of second derivative data

    ddy = np.diff(y_data, 2)

    ddy = ddy / np.max(ddy)

    # find peaks

    classification, sigma = baseline_find_signal(-1 * ddy, corr_distance, False, 2)

    ddy1 = np.roll(ddy, 1)

    ddyn1 = np.roll(ddy, -1)

    p = np.where((ddy < ddy1) & (ddy < ddyn1))[0]

    peaks = p[classification[p] == 1]

    peaks1 = np.roll(peaks, 1)

    distance = np.min(abs(peaks1 - peaks))

    # must make sure the convolution kernel is odd in length to prevent the movement of the peaks

    peaks = np.sort(peaks)

    peakscopy = copy.copy(peaks)
    ddycopy = copy.copy(ddy[peaks] / np.max(ddy))

    while distance < corr_distance:
        # roll the peaks one forward
        peakscopy1 = np.roll(peakscopy, 1)

        # find distances between peaks
        diff = np.abs(peakscopy - peakscopy1)

        # find where in the array the smallest distance is
        mindist = np.argmin(diff)

        # what is this distance
        distance = diff[mindist]

        # compare the values of the second derivative at the closest two peaks

        compare = np.argmax(ddycopy[[mindist, mindist - 1]])

        peakscopy = np.delete(peakscopy, mindist - compare)
        ddycopy = np.delete(ddycopy, mindist - compare)

    # remove any peaks that fall into the noise

    n = y_data[peakscopy]

    w = n > 5 * std

    peakscopy = peakscopy[w]

    final_peaks = sorted(list(peakscopy))

    # draw new regions symmetrically around the newly found peaks

    dist_hz = uc(0, "Hz") - uc(9, "Hz")

    print("     resetting region boundries, distance = " + str(dist_hz))

    peak_regions = []

    for peak in final_peaks:
        l = np.arange(peak + 1, min(peak + dist_hz + 1, len(y_data))).tolist()

        m = np.arange(max(peak - dist_hz, 0), peak).tolist()

        region = m + [peak] + l

        peak_regions.append(region)

    final_regions = [peak_regions[0]]
    final_peaks_seperated = [[final_peaks[0]]]

    for region in range(1, len(peak_regions)):

        if peak_regions[region][0] <= final_regions[-1][-1]:

            final_regions[-1] += peak_regions[region]

            final_peaks_seperated[-1].append(final_peaks[region])

        else:

            final_regions += [peak_regions[region]]

            final_peaks_seperated.append([final_peaks[region]])

    final_regions = [np.arange(min(region), max(region) + 1).tolist() for region in final_regions]

    return final_peaks, final_regions, final_peaks_seperated, std


def multiproc_BIC_minimisation(peak_regions, grouped_peaks, total_spectral_ydata, corr_distance, uc, std):
    maxproc = 5

    pool = mp.Pool(maxproc)

    new_grouped_peaks = [[] for i in range(len(peak_regions))]
    new_grouped_params = [[] for i in range(len(peak_regions))]
    new_sim_y = [[] for i in range(len(peak_regions))]

    def BIC_minimisation_region_full(ind1, uc, peak_regions, grouped_peaks, total_spectral_ydata, corr_distance, std):

        ################################################################################################################
        # initialise process
        ################################################################################################################

        # print("minimising region " + str(ind1) + " of " + str(len(peak_regions)))

        BIC_param = 15

        region = np.array(peak_regions[ind1])

        region_y = total_spectral_ydata[region]

        fit_y = np.zeros(len(region_y))

        copy_peaks = np.array(grouped_peaks[ind1])

        params = Parameters()

        fitted_peaks = []

        ttotal = 0

        ################################################################################################################
        # build initial model
        ################################################################################################################

        # params.add('vregion' + str(ind1), value=2.5, max=5, min=1)

        params.add('vregion' + str(ind1), value=0.5, max=1, min=0)

        distance = uc(0, "hz") - uc(5, "hz")

        std_upper = uc(0, "hz") - uc(1, "hz")
        av_std = uc(0, "hz") - uc(0.2, "hz")
        std_lower = uc(0, "hz") - uc(0.1, "hz")

        # build model

        while (len(copy_peaks) > 0):
            # pick peak that is furthest from fitted data:

            diff_array = region_y - fit_y

            ind2 = np.argmax(diff_array[copy_peaks - region[0]])

            maxpeak = copy_peaks[ind2]

            copy_peaks = np.delete(copy_peaks, ind2)

            # only allow params < distance away vary at a time

            # add new params

            fitted_peaks.append(maxpeak)

            fitted_peaks = sorted(fitted_peaks)

            params.add('A' + str(maxpeak), value=total_spectral_ydata[maxpeak], min=0, max=1, vary=True)

            # params.add('std' + str(maxpeak), value=av_std, vary=True, min = std_lower,
            #               max = std_upper)

            params.add('std' + str(maxpeak), value=av_std, vary=True)

            params.add('mu' + str(maxpeak), value=maxpeak, vary=True
                       , min=maxpeak - 4 * corr_distance, max=maxpeak + 4 * corr_distance)

            # adjust amplitudes and widths of the current model

        initial_y = p7sim(params, region, fitted_peaks, ind1)

        inty = np.sum(region_y[region_y > 0])

        intmodel = np.sum(initial_y)

        # check the region can be optimised this way

        # find peak with max amplitude

        maxamp = 0

        for peak in fitted_peaks:
            amp = params['A' + str(peak)]
            if amp > maxamp:
                maxamp = copy.copy(amp)

        maxintegral = maxamp * len(region)

        if maxintegral > inty:

            # set initial conditions

            while (intmodel / inty < 0.99) or (intmodel / inty > 1.01):

                for f in fitted_peaks:
                    params['std' + str(f)].set(value=params['std' + str(f)] * inty / intmodel)

                initial_y = p7sim(params, region, fitted_peaks, ind1)

                for f in fitted_peaks:
                    params['A' + str(f)].set(
                        value=params['A' + str(f)] * region_y[int(params['mu' + str(f)]) - region[0]] / (
                            initial_y[f - region[0]]))

                initial_y = p7sim(params, region, fitted_peaks, ind1)

                intmodel = np.sum(initial_y)

        # print('built model region ' + str(ind1))

        ################################################################################################################
        # now relax all params
        ################################################################################################################

        # allow all params to vary

        params['vregion' + str(ind1)].set(vary=True)

        for peak in fitted_peaks:
            params['A' + str(peak)].set(vary=False, min=max(0, params['A' + str(peak)] - 0.01),
                                        max=min(params['A' + str(peak)] + 0.01, 1))
            params['mu' + str(peak)].set(vary=False)
            params['std' + str(peak)].set(vary=False, min=min(std_lower, params['std' + str(peak)] - av_std),
                                          max=max(params['std' + str(peak)] + av_std, std_upper))

        out = Minimizer(p7residual, params,
                        fcn_args=(region, fitted_peaks, region_y, ind1, False))

        results = out.minimize()

        params = results.params

        # print('relaxed params region ' + str(ind1))

        ################################################################################################################
        # now remove peaks in turn
        ################################################################################################################

        trial_y = p7sim(params, region, fitted_peaks, ind1)

        trial_peaks = np.array(fitted_peaks)

        amps = []

        for peak in trial_peaks:
            amps.append(params['A' + str(peak)])

        r = trial_y - region_y

        chi2 = r ** 2

        N = len(chi2)

        BIC = N * np.log(np.sum(chi2) / N) + np.log(N) * (3 * len(fitted_peaks) + 2)

        while (len(trial_peaks) > 0):

            new_params = copy.copy(params)

            # find peak with smallest amp

            minpeak = trial_peaks[np.argmin(amps)]

            # remove this peak from the set left to try

            trial_peaks = np.delete(trial_peaks, np.argmin(amps))
            amps = np.delete(amps, np.argmin(amps))

            # remove this peak from the trial peaks list and the trial params

            new_params.__delitem__('A' + str(minpeak))
            new_params.__delitem__('mu' + str(minpeak))
            new_params.__delitem__('std' + str(minpeak))

            new_fitted_peaks = np.delete(fitted_peaks, np.where(fitted_peaks == minpeak))

            # simulate data with one fewer peak

            new_trial_y = p7sim(new_params, region, new_fitted_peaks, ind1)

            r = new_trial_y - region_y

            chi2 = np.sum(r ** 2)

            N = len(new_trial_y)

            new_BIC = N * np.log(chi2 / N) + np.log(N) * (3 * len(new_fitted_peaks) + 2)

            # if the fit is significantly better remove this peak

            if new_BIC < BIC - BIC_param:
                fitted_peaks = copy.copy(new_fitted_peaks)

                params = copy.copy(new_params)

                BIC = copy.copy(new_BIC)

        fitted_peaks = sorted(fitted_peaks)

        fit_y = p7sim(params, region, fitted_peaks, ind1)

        ################################################################################################################

        print("     done region " + str(ind1 + 1) + "of" + str(len(peak_regions)))

        return fitted_peaks, params, fit_y

    # order regions by size to efficiently fill cores

    region_lengths = np.array([len(g) for g in grouped_peaks])

    sorted_regions = np.argsort(region_lengths)[::-1]

    res = [[] for i in peak_regions]

    # write output files

    for ind1 in range(len(peak_regions)):
        res[ind1] = pool.apply_async(BIC_minimisation_region_full,
                                     [ind1, uc, peak_regions, grouped_peaks, total_spectral_ydata, corr_distance,
                                      std])

    for ind1 in sorted_regions:
        new_grouped_peaks[ind1], new_grouped_params[ind1], new_sim_y[ind1] = res[ind1].get()

    #### unpack the parameters and split groups

    final_grouped_peaks = []
    final_peak_regions = []
    final_sim_y = []
    final_peaks = []

    total_params = Parameters()

    new_peaks = []

    dist_hz = uc(0, "Hz") - uc(20, "Hz")

    newgroupind = 0
    oldgroupind = 0

    for group in new_grouped_peaks:

        group = sorted(group)

        if len(group) > 0:
            final_grouped_peaks.append([])
            total_params.add('vregion' + str(newgroupind),
                             value=new_grouped_params[oldgroupind]['vregion' + str(oldgroupind)])

            final_peaks.extend(group)

        for ind2, peak in enumerate(group):

            # check if the group should be split

            if ind2 > 0:

                # if there is a gap of more than 20Hz split the group

                if peak > group[ind2 - 1] + dist_hz:
                    # track the numer of splits included to ensure v parameter is added to the correct group each time

                    newgroupind += 1

                    # if a split occures add a new v parameter for the new group

                    total_params.add('vregion' + str(newgroupind),
                                     value=new_grouped_params[oldgroupind]['vregion' + str(oldgroupind)])

                    # allow peaks to be added to the new group
                    final_grouped_peaks.append([])

            # finally append the peak
            final_grouped_peaks[-1].append(peak)
            total_params.add('A' + str(peak), value=new_grouped_params[oldgroupind]['A' + str(peak)])
            total_params.add('std' + str(peak), value=new_grouped_params[oldgroupind]['std' + str(peak)])
            total_params.add('mu' + str(peak), value=new_grouped_params[oldgroupind]['mu' + str(peak)], vary=False)

        if len(group) > 0:
            newgroupind += 1

        oldgroupind += 1

    # draw regions between midpoints of groups

    for ind4, group in enumerate(final_grouped_peaks):

        if ind4 == 0:
            lower_point = 0
            higher_point = int((group[-1] + final_grouped_peaks[ind4 + 1][0]) / 2)

        elif ind4 == len(final_grouped_peaks) - 1:
            lower_point = int((group[0] + final_grouped_peaks[ind4 - 1][-1]) / 2)
            higher_point = len(total_spectral_ydata)

        else:
            lower_point = int((group[0] + final_grouped_peaks[ind4 - 1][-1]) / 2)

            higher_point = int((group[-1] + final_grouped_peaks[ind4 + 1][0]) / 2)

        final_peak_regions.append(np.arange(lower_point, higher_point))

    # now simulate new regions and store region data

    for ind3, region in enumerate(final_peak_regions):
        fit_y = p7sim(total_params, region, final_grouped_peaks[ind3], ind3)

        final_sim_y.append(fit_y)

    return final_peaks, final_grouped_peaks, final_peak_regions, final_sim_y, total_params


########################################################################################################################
# integration and processing
########################################################################################################################


def lorentzian(p, w, p0, A):
    x = (p0 - p) / (w / 2)

    L = A / (1 + x ** 2)

    return L


def lorenz_curves(params, x, picked_points):
    y = np.zeros(len(x))
    for peak in picked_points:
        y += lorentzian(x, params['width' + str(peak)], params['pos' + str(peak)], params['amp' + str(peak)])
    return y


def simulate_regions(params, peak_regions, grouped_peaks, y_data, xppm):
    sim_regions = []
    sim_y = np.zeros(len(y_data))

    for ind, region in enumerate(peak_regions):
        sim_y[region] = p7plot(params, region, grouped_peaks[ind], ind, xppm)

        y = p7sim(params, region, grouped_peaks[ind], ind)

        sim_regions.append(y)

    return sim_regions, sim_y


def new_first_order_peak(start_ppm, J_vals, x_data, corr_distance, uc, spin):
    # new first order peak generator using the method presented in Hoye paper

    start = uc(str(start_ppm) + "ppm")

    start_Hz = uc.hz(start)

    J_vals = np.array(J_vals)

    peaks = np.zeros((2 * spin + 1) ** len(J_vals))

    if spin == 0.5:
        l = [1, -1]

    if spin == 1:
        l = [1, 0, -1]

    # signvector generator

    signvectors = itertools.product(l, repeat=len(J_vals))

    for ind, sv in enumerate(signvectors):
        shift = J_vals * sv

        shift = start_Hz + np.sum(shift)
        peaks[ind] = shift

    peaks = np.sort(peaks)

    peak_vector = np.array(sorted(list(set(peaks)), reverse=True))

    amp_vector = np.zeros(len(peak_vector))

    for peak in peaks:
        index = np.where(peak_vector == peak)
        amp_vector[index] += 1

    pv = []

    for index, peak in enumerate(peak_vector):
        pv.append(uc(peak, "Hz"))

    peak_vector = pv

    split_params = Parameters()

    for index, peak in enumerate(peak_vector):
        split_params.add('amp' + str(peak), value=amp_vector[index])
        split_params.add('pos' + str(peak), value=peak)
        split_params.add('width' + str(peak), value=2 * corr_distance)

    y = lorenz_curves(split_params, x_data, peak_vector)

    y = y / np.max(y)

    # where = np.where(y > 0.001)

    # y = y[where]

    return peak_vector, amp_vector, y

def editsolvent_removal2(solvent, y_data, x_data, picked_peaks, peak_regions, grouped_peaks, total_params, uc):
    picked_peaks = np.array(picked_peaks)

    # define the solvents

    if solvent == 'chloroform':

        exp_ppm = [7.26]

        Jv = [[]]

    elif solvent == 'dimethylsulfoxide':

        exp_ppm = [2.50]

        Jv = [[1.9, 1.9]]

    elif solvent == 'methanol':

        exp_ppm = [4.78, 3.31]

        Jv = [[], [1.7, 1.7]]

    elif solvent == 'benzene':

        exp_ppm = [7.16]

        Jv = [[]]

    elif solvent == 'pyridine':

        exp_ppm = [8.74, 7.58, 7.22]

        Jv = [[], [], []]

    else:
        exp_ppm = []
        Jv = [[]]

    # find picked peaks ppm values

    picked_peaks_ppm = x_data[picked_peaks]

    # make differences vector for referencing against multiple solvent peaks

    differences = []

    peaks_to_remove = []

    solvent_regions = []

    # now remove each peak in turn

    for ind1, speak_ppm in enumerate(exp_ppm):

        # if only a singlet is expected for this peak find solvent peak based on amplitude and position

        if len(Jv[ind1]) == 0:

            probs = norm.pdf(abs(picked_peaks_ppm - speak_ppm), loc=0, scale=0.1) * y_data[picked_peaks]

            # find the maximum probability

            w = np.argmax(probs)

            # append this to the list to remove

            peaks_to_remove.append(picked_peaks[w])

            # append this to the list of differences

            differences.append(speak_ppm - picked_peaks_ppm[w])

        # if the peak displays a splitting pattern then we have to be a bit more selective
        # do optimisation problem with projected peaks
        else:

            amp_res = []
            dist_res = []
            pos_res = []

            # limit the search to peaks +- 1 ppm either side

            srange = (picked_peaks_ppm > speak_ppm - 1) * (picked_peaks_ppm < speak_ppm + 1)

            for peak in picked_peaks_ppm[srange]:

                # print("picked ppm ", peak)

                fit_s_peaks, amp_vector, fit_s_y = new_first_order_peak(peak, Jv[ind1], np.arange(len(x_data)), 0.1, uc,
                                                                        1)

                diff_matrix = np.zeros((len(fit_s_peaks), len(picked_peaks)))

                for i, f in enumerate(fit_s_peaks):

                    for j, g in enumerate(picked_peaks):
                        diff_matrix[i, j] = abs(f - g)

                # minimise these distances

                vertical_ind, horizontal_ind = optimise(diff_matrix)

                closest_peaks = np.sort(picked_peaks[horizontal_ind])

                closest_amps = []

                for cpeak in closest_peaks:
                    closest_amps.append(total_params['A' + str(cpeak)])

                # find the amplitude residual between the closest peaks and the predicted pattern

                # normalise these amplitudes

                amp_vector = [i / max(amp_vector) for i in amp_vector]

                closest_amps = [i / max(closest_amps) for i in closest_amps]

                # append to the vector

                amp_res.append(sum([abs(amp_vector[i] - closest_amps[i]) for i in range(len(amp_vector))]))

                dist_res.append(np.sum(np.abs(closest_peaks - fit_s_peaks)))

                pos_res.append(norm.pdf(abs(peak - speak_ppm), loc=0, scale=0.5))

                # use the gsd data to find amplitudes of these peaks

            pos_res = [1 - i / max(pos_res) for i in pos_res]

            dist_res = [i / max(dist_res) for i in dist_res]

            amp_res = [i / max(amp_res) for i in amp_res]

            # calculate geometric mean of metrics for each peak

            g_mean = [(dist_res[i] + amp_res[i] + pos_res[i]) / 3 for i in range(0, len(amp_res))]

            # compare the residuals and find the minimum

            minres = np.argmin(g_mean)

            # append the closest peaks to the vector

            fit_s_peaks, amp_vector, fit_s_y = new_first_order_peak(picked_peaks_ppm[srange][minres], Jv[ind1],
                                                                    np.arange(len(x_data)), 0.1, uc, 1)

            diff_matrix = np.zeros((len(fit_s_peaks), len(picked_peaks)))

            for i, f in enumerate(fit_s_peaks):

                for j, g in enumerate(picked_peaks):
                    diff_matrix[i, j] = abs(f - g)

            # minimise these distances

            vertical_ind, horizontal_ind = optimise(diff_matrix)

            closest_peaks = np.sort(picked_peaks[horizontal_ind])

            for peak in closest_peaks:
                ind3 = np.abs(picked_peaks - peak).argmin()

                peaks_to_remove.append(picked_peaks[ind3])

                differences.append(picked_peaks_ppm[ind3] - uc.ppm(peak))

    # find the region this peak is in and append it to the list

    for peak in peaks_to_remove:

        for ind2, region in enumerate(peak_regions):

            if (peak > region[0]) & (peak < region[-1]):
                solvent_regions.append(ind2)
                break

    # now remove the selected peaks from the picked peaks list and grouped peaks

    w = np.searchsorted(picked_peaks, peaks_to_remove)

    picked_peaks = np.delete(picked_peaks, w)

    for ind4, peak in enumerate(peaks_to_remove):
        grouped_peaks[solvent_regions[ind4]] = np.delete(grouped_peaks[solvent_regions[ind4]],
                                                         np.where(grouped_peaks[solvent_regions[ind4]] == peak))

    # resimulate the solvent regions

    solvent_region_ind = sorted(list(set(solvent_regions)))

    # now need to reference the spectrum

    # differences = list of differences in ppm found_solvent_peaks - expected_solvent_peaks

    s_differences = sum(differences)

    x_data = x_data + s_differences

    return peak_regions, picked_peaks, grouped_peaks, x_data, solvent_region_ind


def integrate_regions(peak_regions, y_data, noise_std):
    integrals = np.zeros(len(peak_regions))

    for index, region in enumerate(peak_regions):
        w = np.where(y_data[region] > 3 * noise_std)

        integrals[index] = np.sum(y_data[region][w])

    return integrals


def integrate_sim_regions(sim_regions, grouped_peaks, peak_regions, y_data, params, solvent_region_ind):
    sim_integrals = []

    for r, group in enumerate(grouped_peaks):

        region_integral = 0

        for peak in group:
            region_integral += params['A' + str(peak)] * 0.25 * np.pi * params['std' + str(peak)] * (
                        (3 ** 0.5 - 2) * params['vregion' + str(r)] + 2)

        sim_integrals.append(region_integral)

    sim_integrals = np.array(sim_integrals)

    print("exact", list(np.round(sim_integrals, 2)))

    y_integral = np.sum(y_data)

    sim_integral = np.sum(sim_integrals)

    k = sim_integral / y_integral

    integrals = []

    for region in peak_regions:
        integrals.append(np.sum(y_data[region]))

    k = np.array(sim_integrals) / np.array(integrals)

    simr_regions = []

    for region in sim_regions:
        simr_regions.append(np.sum(region))

    integrals = k * integrals

    return integrals


def integral_add(sim_regions, proton_guess):

    integrals = [np.sum(sim) for sim in sim_regions]

    i_sum = np.sum(integrals)

    Inorm = proton_guess / i_sum

    integrals_copy = list(integrals)

    integrals_copy = [0] + integrals_copy

    integral_sum = np.cumsum(np.asarray(integrals_copy))
    integral_sum = integral_sum / np.sum(integrals_copy)

    cummulative_vectors = []

    t_sum = 0
    for region in sim_regions:
        t_sum += np.sum(region)
        cummulative_vectors.append(np.cumsum(region))

    for i, region in enumerate(cummulative_vectors):
        cummulative_vectors[i] = region / t_sum

    return integral_sum, cummulative_vectors


def normalise_integration(integrals, initial_proton_guess):
    i_sum = np.sum(integrals)

    integrals = integrals / i_sum

    norm_integrals = integrals * initial_proton_guess

    return norm_integrals


def remove_impurities(integrals, peak_regions, grouped_peaks, picked_peaks, sim_regions):
    # find rounded values

    r = sum_round(integrals)

    r = np.array(r)

    to_remove = np.where(r < 0.5)[0]

    number_of_impurities = len(to_remove)

    peaks_to_remove = []

    for group in to_remove:
        peaks_to_remove.append(grouped_peaks[group])

    whdel = np.where(picked_peaks == peaks_to_remove)

    picked_peaks = np.delete(picked_peaks, whdel)

    integrals = np.delete(integrals, to_remove)

    peak_regions = np.delete(peak_regions, to_remove)

    grouped_peaks = np.delete(grouped_peaks, to_remove)

    sim_regions = np.delete(sim_regions, to_remove)

    rounded_integrals = r[r > 0.5]

    return grouped_peaks, integrals, peak_regions, picked_peaks, number_of_impurities, sim_regions,rounded_integrals


def proton_count(file):
    sdf_file = open(str(file), 'r')

    txt = sdf_file.readlines()

    counter = 0
    for line in txt:
        if ' H   ' in line:
            counter += 1

    return counter


def integral_score(integrals, structure_protons, proton_guess, l_protons, impurities):

    r_int = sum_round(integrals)

    r_int = np.array(r_int)

    sum_r = np.sum(r_int)

    new_integrals = []

    for ind, integral in enumerate(integrals):
        new_integrals += [integral] * int(r_int[ind])

    new_integrals = np.array(new_integrals)

    differences = new_integrals % 0.5

    std = structure_protons / 8

    diff = abs(proton_guess - structure_protons)

    probs = 4 * (1 - norm.cdf(differences, loc=0, scale=1 / 16))

    mean = gmean(probs) * (1 - norm.cdf(diff, loc=0, scale=std)) * (1 / 2 ** impurities)

    # only allow intergrals that are between the expected number and that number - the number of labile protons

    if (sum_r < structure_protons - l_protons):
        mean = 0

    return mean


def labile_protons(file):
    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()
    obconversion.ReadFile(obmol, file)

    count = 0

    for atom in OBMolAtomIter(obmol):
        for NbrAtom in OBAtomAtomIter(atom):
            if (atom.GetAtomicNum() == 8) & (NbrAtom.GetAtomicNum() == 1):
                count += 1

    return count


def find_integrals(file, peak_regions, grouped_peaks, sim_regions,
                   picked_peaks, params, y_data, solvent_region_ind):
    # count number of protons in the file

    structure_protons = proton_count(str(file).split('.sdf')[0] + ".sdf")

    # find number of methyl groups

    m = methyl_protons(str(file).split('.sdf')[0] + ".sdf")

    number_of_methyl_groups = len(m)

    print('     number of protons = ' + str(structure_protons))

    # count the number of labile protons in the structure

    l_protons = labile_protons(file)

    print('     number of labile protons = ' + str(l_protons))

    count = 0

    # allow guesses of number of protons in the spectrum between: structure_protons - l_protons, 2 * structure_protons

    number_vector = np.arange(structure_protons - l_protons, 2 * structure_protons)

    scores = np.zeros(len(number_vector))

    integrals = integrate_sim_regions(sim_regions, grouped_peaks, peak_regions, y_data, params, solvent_region_ind)

    for proton_guess in number_vector:

        norm_integrals = normalise_integration(integrals, proton_guess)

        # remove impurities

        # find number of impurities removed

        # find which should be classed as impurities

        r = sum_round(norm_integrals)

        r = np.array(r)

        impurities = len(r[r < 1])

        norm_integrals = norm_integrals[r > 0.5]

        scores[count] = integral_score(norm_integrals, structure_protons, proton_guess, l_protons, impurities)

        number_of_methyl_groups_integral = np.sum((r - (r % 3)) // 3)

        if number_of_methyl_groups_integral < number_of_methyl_groups:
            scores[count] = 0

        count += 1

    wh = np.argmax(scores)

    best_fit = number_vector[wh]

    print("the best fit number of protons is " + str(best_fit))

    # normalise using this number

    integrals = normalise_integration(integrals, best_fit)

    grouped_peaks, integrals, peak_regions, picked_peaks_, impurities, sim_regions,rounded_integrals = remove_impurities(integrals,
                                                                                                       peak_regions,
                                                                                                       grouped_peaks,
                                                                                                       picked_peaks,
                                                                                                       sim_regions)

    integral_sum, cummulative_vectors = integral_add(sim_regions, best_fit)

    total_r = 0

    total = 0

    for i in range(0, len(integrals)):
        total += integrals[i]
        total_r += rounded_integrals[i]

    return peak_regions, grouped_peaks, sim_regions, integral_sum, cummulative_vectors, rounded_integrals, structure_protons, \
           best_fit, total_r


def sum_round(a):
    rounded = [round(i, 0) for i in a]

    error = sum(a) - sum(rounded)

    n = int(round(error / 1))

    new = rounded[:]

    for _, i in sorted(((a[i] - rounded[i], i) for i in range(len(a))), reverse=n > 0)[:abs(n)]:
        new[i] += n / abs(n)

    return new


def weighted_region_centres(peak_regions, total_spectral_ydata):
    centres = []

    for region in peak_regions:
        w = total_spectral_ydata[region] ** 2

        wx = region * w

        xbar = np.sum(wx) / np.sum(total_spectral_ydata[region] ** 2)

        centres.append(int(xbar))

    return centres


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

