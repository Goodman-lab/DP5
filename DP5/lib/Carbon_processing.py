import nmrglue as ng
from scipy.stats import gaussian_kde as kde
from scipy.optimize import curve_fit
import re
import numpy as np
from scipy.stats import norm
from scipy.ndimage.filters import gaussian_filter1d as g1d
from scipy.ndimage.filters import convolve1d as c1d
from matplotlib import pyplot as plt
import pickle
from lmfit import Minimizer, Parameters, report_fit
import copy as copy
import itertools
import statsmodels.api as sm
import os


def process_carbon(NMR_file,settings,datatype):

    total_spectral_ydata, spectral_ydata, spectral_xdata_ppm, threshold, corr_distance, uc = spectral_processing(
        NMR_file,datatype)

    total_spectral_ydata = edge_removal(total_spectral_ydata)

    picked_peaks, simulated_ydata = iterative_peak_picking(total_spectral_ydata, 5, corr_distance)

    picked_peaks = sorted(list(set(picked_peaks)))

    picked_peaks, removed = solvent_removal(simulated_ydata, spectral_xdata_ppm, settings.Solvent, uc, picked_peaks)

    return total_spectral_ydata,spectral_xdata_ppm,corr_distance,uc,picked_peaks,simulated_ydata,removed

########################################################################################################################
# processing
########################################################################################################################

def spectral_processing(file,datatype):

    print("Processing Carbon Spectrum")

    spectral_xdata_ppm, total_spectral_ydata, uc = initial_processing(file,datatype)

    corr_distance = estimate_autocorrelation(total_spectral_ydata)

    convolved_y = gaussian_convolution(corr_distance,total_spectral_ydata)
    binary_map_regions =[]
    #threshold_vector =[4,3.7,3.5,3,2,1,0.9,0.8,0.7,0.6,0.5]
    threshold_vector =[3,2.9,2.8,2.7,2.6,2.5,2.2,2,1,0.9,0.8,0.7,0.6,0.5]
    run = 0

    while len(binary_map_regions) < 2:
        threshold = threshold_vector[run]
        run += 1
        picked_points = iterative_point_picking(convolved_y,threshold)
        binary_map_regions,binary_map_list = binary_map(picked_points, uc,convolved_y)

    globalangles, phased_peak_regions, convolved_y_phased = estimate_phase_angles(convolved_y, binary_map_regions, corr_distance)

    real_convolved_y_phased = list(np.real(convolved_y_phased))

    picked_points_region = iterative_point_picking_region(binary_map_regions,real_convolved_y_phased,threshold)

    picked_peaks_region = peak_picking_region(real_convolved_y_phased,picked_points_region)

    p0, p1 = linear_regression(picked_peaks_region, globalangles,real_convolved_y_phased,binary_map_regions)

    total_spectral_ydata,spectral_ydata = final_phasing(convolved_y, p0, p1)

    total_spectral_ydata = total_spectral_ydata/np.max(total_spectral_ydata)

    return total_spectral_ydata,spectral_ydata,spectral_xdata_ppm,threshold,corr_distance,uc

def jcamp_guess_udic(dic, data):
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

        udic[0]['car'] = firstx  -  abs(lastx - firstx)/2

    return udic

def initial_processing(file,datatype):

    if datatype =='jcamp':

        dic, total_spectral_ydata = ng.jcampdx.read(file)  # read file

        total_spectral_ydata = total_spectral_ydata[0] + 1j * total_spectral_ydata[1]

        total_spectral_ydata = ng.proc_base.ifft_positive(total_spectral_ydata)

    else:

        dic, total_spectral_ydata = ng.bruker.read(file)  # read file

        total_spectral_ydata = ng.bruker.remove_digital_filter(dic, total_spectral_ydata)  # remove the digital filter

    total_spectral_ydata = ng.proc_base.zf_double(total_spectral_ydata, 2)  # zero filling once

    total_spectral_ydata = ng.proc_base.fft_positive(total_spectral_ydata)  # Fourier transform

    real_part = ng.proc_bl.baseline_corrector(np.real(total_spectral_ydata), wd=2)
    im_part = ng.proc_bl.baseline_corrector(np.imag(total_spectral_ydata), wd=2)

    total_spectral_ydata = real_part + 1j *im_part

    if datatype == 'jcamp':

        udic = jcamp_guess_udic(dic, total_spectral_ydata)

    else:

        udic = ng.bruker.guess_udic(dic, total_spectral_ydata)  # sorting units

    # total_spectral_ydata = ng.proc_autophase.autops(total_spectral_ydata, 'acme')  # automatic phase correction

    uc = ng.fileiobase.uc_from_udic(udic)  # unit conversion element
    spectral_xdata_ppm = uc.ppm_scale()  # ppmscale creation

    maximum= np.max(total_spectral_ydata)

    total_spectral_ydata = total_spectral_ydata/np.max(total_spectral_ydata)

    return spectral_xdata_ppm, total_spectral_ydata, uc

def estimate_autocorrelation(total_spectral_ydata):

    real_part = np.real(total_spectral_ydata)
    real_part_copy = np.real(total_spectral_ydata)

    gzero = real_part * real_part
    gzero = np.sum(gzero)

    gdx = gzero

    counter = 0

    while gdx > 0.6 * gzero:
        real_part_copy = np.roll(real_part_copy, counter)
        gdx = real_part * real_part_copy
        counter += 1
        gdx = np.sum(gdx)

    corr_distance = counter


    return corr_distance

def gaussian_convolution(corr_distance, total_spectral_ydata):

    real_part = np.real(total_spectral_ydata)
    im_part = np.imag(total_spectral_ydata)


    real_convolved_y = g1d(real_part, corr_distance)
    im_convolved_y = g1d(im_part, corr_distance)

    convolved_y = np.array(real_convolved_y) + 1j * np.array(im_convolved_y)

    convolved_y = convolved_y / np.max(convolved_y)

    return convolved_y

def lorentz_convolution(corr_distance, total_spectral_ydata):

    def lorentzian(p, w, p0):
        x = (p0 - p) / (w / 2)
        L = 1 / (1 + x ** 2)
        return L

    def build_kernel(corr):
        kernel_length = corr * 10
        vector = np.arange(0, kernel_length)
        p0 = kernel_length // 2
        kernel = lorentzian(vector, corr, p0)
        return kernel


    real_part = np.real(total_spectral_ydata)
    im_part = np.imag(total_spectral_ydata)
    kernel = build_kernel(corr_distance)

    real_convolved_y = c1d(real_part, kernel)
    im_convolved_y = c1d(im_part, kernel)

    real_convolved_y = c1d(real_convolved_y, -1 * kernel)
    im_convolved_y = c1d(im_convolved_y, -1 * kernel)

    convolved_y = np.array(real_convolved_y) + 1j * np.array(im_convolved_y)

    return convolved_y

def iterative_point_picking(convolved_y, threshold):

    real_convolved_y = np.real(convolved_y)
    copy_convolved_y = np.array(real_convolved_y)
    picked_points = []
    pickednumber = 1
    while pickednumber > 0:
        mu, std = norm.fit(copy_convolved_y)

        index = np.where(copy_convolved_y - mu > threshold * std)
        pickednumber = len(index[0])
        picked_points.extend(np.ndarray.tolist(index[0]))
        copy_convolved_y = np.delete(copy_convolved_y, index, axis=0)

    copy_convolved_y = np.array(real_convolved_y)
    pickednumber = 1
    while pickednumber > 0:
        mu, std = norm.fit(copy_convolved_y)
        index = np.where(copy_convolved_y - mu < -threshold * std)
        pickednumber = len(index[0])
        picked_points.extend(np.ndarray.tolist(index[0]))
        copy_convolved_y = np.delete(copy_convolved_y, index, axis=0)

    picked_points = sorted(picked_points)

    return picked_points

def binary_map(picked_points, uc, convolved_y):


    picked_points = np.array(picked_points)

    # find where peak blocks are

    binary_map_regions = [[picked_points[0]]]

    for x in range(0, len(picked_points) - 1):
        if picked_points[x + 1] != picked_points[x] + 1:
            binary_map_regions[-1].append(picked_points[x])
            binary_map_regions.append([picked_points[x + 1]])
    binary_map_regions[-1].append(picked_points[-1])

    # extend blocks by 50 Hz

    for block in binary_map_regions:
        start = uc.hz(block[0])
        end = start + 50
        end_point = uc(end, "Hz")
        block[0] = end_point

        start = uc.hz(block[1])
        end = start - 50
        end_point = uc(end, "Hz")
        block[1] = end_point

    # draw binary map

    binary_map_list = np.zeros(len(convolved_y))
    for block in binary_map_regions:
        binary_map_list[block[0]:block[1]:1] = 1

    # stitch blocks together
    blocks = np.where(binary_map_list == 1)

    blocks = blocks[0] - 1
    binary_map_regions = [[blocks[0]]]

    for element in range(0, len(blocks) - 1):
        if blocks[element + 1] != blocks[element] + 1:
            binary_map_regions[-1].append(blocks[element])
            binary_map_regions.append([blocks[element + 1]])
    binary_map_regions[-1].append(blocks[-1])

    return binary_map_regions, binary_map_list

def estimate_phase_angles(convolved_y, binary_map_regions, corr_distance):

    convolved_y_phased = np.array(convolved_y)

    def inte(binary_map_regions, peak_regions, corr_distance):
        ## for each region determine the baseline
        integrals = [0] * len(binary_map_regions)

        # first find average of surrounding points to ends of binary map regions to draw base line
        baselines_end = [[0, 0] for i in range(len(binary_map_regions))]
        baselines = [[] for i in range(0, len(binary_map_regions))]

        # find baseline endpoints
        for region in range(0, len(binary_map_regions)):
            for point in range(0, corr_distance - 1):
                baselines_end[region][0] += peak_regions[region][point]

            baselines_end[region][0] = baselines_end[region][0] / corr_distance
            for point in range(0, corr_distance):
                baselines_end[region][1] += peak_regions[region][-point]

            baselines_end[region][1] = baselines_end[region][1] / corr_distance

            ## draw baselines
            baselines[region] = np.linspace(baselines_end[region][0], baselines_end[region][1],
                                            len(peak_regions[region]) - 2 * corr_distance)

            ## integrate each region below the baseline
            for point in range(0, len(baselines[region])):
                if peak_regions[region][point + corr_distance] < baselines[region][point]:
                    integrals[region] += abs(
                        peak_regions[region][point + corr_distance] - baselines[region][point])

        return integrals

    coarse_angle = np.linspace(-np.pi / 2, np.pi / 2, 1000)
    integral_vector = [0] * 1000
    counter = 0
    # integration
    for angle in coarse_angle:
        copy_total_spectral_ydata = convolved_y * np.exp(-angle * 1j)
        peak_regions = [0] * len(binary_map_regions)
        for region in range(0, len(binary_map_regions)):
            peak_regions[region] = copy_total_spectral_ydata[
                                   binary_map_regions[region][0]:binary_map_regions[region][1]:1]

        integral_vector[counter] = inte(binary_map_regions, peak_regions, corr_distance)
        counter = counter + 1

    # find maximum integral for each region and store angles
    integral_vector = np.array(integral_vector)
    maxvector = np.amin(integral_vector, 0)

    counter = 0
    angle1 = [0] * len(binary_map_regions)

    for element in list(maxvector):
        maxangle = np.where(integral_vector == element)
        angle1[counter] = coarse_angle[maxangle[0][0]]
        counter = counter + 1

    # phase each region of the spectrum indepedently

    for region in range(0, len(peak_regions)):
        convolved_y_phased[binary_map_regions[region][0]:binary_map_regions[region][1]:1] = convolved_y_phased[
                                                                                            binary_map_regions[
                                                                                                region][0]:
                                                                                            binary_map_regions[
                                                                                                region][
                                                                                                1]:1] * np.exp(
            -angle1[region] * 1j)

    globalangles = [angle1[i] for i in range(0, len(binary_map_regions))]

    ## phase each peak region separately

    phased_peak_regions = []
    copy_total_spectral_ydata = convolved_y
    peak_regions = [0] * len(binary_map_regions)
    for region in range(0, len(binary_map_regions)):
        peak_regions[region] = copy_total_spectral_ydata[
                               binary_map_regions[region][0]:binary_map_regions[region][1]:1]
    counter = 0

    for region in peak_regions:
        phased_peak_regions.append(region * np.exp(-globalangles[counter] * 1j))
        counter += 1

    return globalangles, phased_peak_regions, convolved_y_phased

def iterative_point_picking_region(binary_map_regions, real_convolved_y_phased, threshold):


    copy_convolved_y = np.array(real_convolved_y_phased)
    picked_points = []
    pickednumber = 1

    while pickednumber > 0:
        mu, std = norm.fit(copy_convolved_y)
        index = np.where((copy_convolved_y - mu > threshold * std) | (copy_convolved_y - mu < threshold * std))
        picked_points.extend(np.ndarray.tolist(index[0]))
        pickednumber = len(index[0])
        copy_convolved_y = np.delete(copy_convolved_y, index, axis=0)

    picked_points = sorted(picked_points)
    picked_points_region = []

    for region in binary_map_regions:
        picked_points_region.append([])
        for point in picked_points:
            if point > region[0] and point < region[1]:
                picked_points_region[-1].append(point)

    return picked_points_region

def peak_picking_region(real_convolved_y_phased, picked_points_region):
    picked_peaks_region = []

    for region in range(0, len(picked_points_region)):
        picked_peaks_region.append([])

        for index in picked_points_region[region]:
            peak = real_convolved_y_phased[index]
            if peak > real_convolved_y_phased[index + 1] and peak > real_convolved_y_phased[index - 1]:
                picked_peaks_region[-1].append(index)
            elif peak < real_convolved_y_phased[index + 1] and peak < real_convolved_y_phased[index - 1]:
                picked_peaks_region[-1].append(index)
    return picked_peaks_region

def linear_regression(picked_peaks_region, globalangles, real_convolved_y_phased, binary_map_regions):

    # region weighting vector

    region_weighting_matrix = [1] * len(picked_peaks_region)

    for index,region in enumerate(picked_peaks_region):
        region_weighting_matrix[index] = max([abs(real_convolved_y_phased[peak]) for peak in region])
    max_weight = max(region_weighting_matrix)

    region_weighting_matrix = [i/max_weight for i in region_weighting_matrix]

    # define centres of regions
    region_centres = []

    for region in binary_map_regions:
        region_centres.append((1 - (region[0] + region[1]) / (2 * len(real_convolved_y_phased))))

    #### regression and outlier analysis

    number_of_outliers = 1
    while number_of_outliers > 0:
        region_centres_regression = sm.add_constant(region_centres)
        wls_model = sm.WLS(globalangles, region_centres_regression,weights=region_weighting_matrix)
        results = wls_model.fit()
        params = results.params
        predictions = [params[1] * i + params[0] for i in region_centres]
        # remove maximum outlier more than 0.6 rad from estimate
        differences = [abs(predictions[angle] - globalangles[angle]) for angle in range(0, len(globalangles))]

        maxdifference = max(differences)
        if maxdifference > 0.6:
            index = differences.index(maxdifference)
            globalangles.pop(index)
            region_centres.pop(index)
            region_weighting_matrix.pop(index)
        else:
            number_of_outliers = 0

    p0 = params[0]
    p1 = params[1]

    return p0, p1

def final_phasing(total_spectral_ydata, p0, p1):

    # total_spectral_ydata = ng.proc_base.ps(total_spectral_ydata, p0=p0, p1=p1)

    relativeposition = np.linspace(1, 0, len(total_spectral_ydata))

    angle = p0 + p1 * relativeposition

    total_spectral_ydata = total_spectral_ydata * np.exp(-1j * angle)

    total_spectral_ydata = ng.proc_bl.baseline_corrector(total_spectral_ydata, wd=2)

    spectral_ydata = ng.proc_base.di(total_spectral_ydata)  # discard the imaginaries
    spectral_ydata = np.ndarray.tolist(spectral_ydata)
    total_spectral_ydata = np.real(total_spectral_ydata)

    return total_spectral_ydata, spectral_ydata

def edge_removal(total_spectral_ydata):

    if total_spectral_ydata[0] > 0:
        i  = 0
        while total_spectral_ydata[i] > 0:
            total_spectral_ydata[i] = 0
            i += 1
    else:
        i = 0
        while total_spectral_ydata[i] < 0:
            total_spectral_ydata[i] = 0
            i += 1

    if total_spectral_ydata[-1] > 0:
        i =1
        while total_spectral_ydata[-i] > 0:
            total_spectral_ydata[-i] = 0
            i += 1
    else:
        i = 1
        while total_spectral_ydata[-i] < 0:
            total_spectral_ydata[-i] = 0
            i += 1
    return total_spectral_ydata

def peak_pruning(picked_peaks, total_spectral_ydata, point_ppm,corr_distance):

    distance = 0.25/point_ppm

    distance = corr_distance * point_ppm

    grouped_peaks = [[picked_peaks[0]]]


    for index in range(0, len(picked_peaks) - 1):
        if picked_peaks[index] + distance > picked_peaks[index + 1]:
            grouped_peaks[-1].append(picked_peaks[index + 1])
        else:
            grouped_peaks.append([picked_peaks[index + 1]])

    new_peaks = []

    for group in grouped_peaks:
        group_amps = total_spectral_ydata[group]

        maxindex = np.argmax(group_amps)

        new_peaks.append(group[maxindex])

    new_peaks = np.array(new_peaks)

    picked_peaks = list(new_peaks)

    picked_peaks = np.array(new_peaks)

    return picked_peaks

def rounding_variables(all_peak_locations_ppm, final_solvent_peak_locations, assigned_peaks_sorted_descending, differences):
    ##creates rounded values
    rounded_picked_locations_ppm = np.around(all_peak_locations_ppm,
                                             decimals=2)  # [round(x, 2) for x in all_peak_locations_ppm]
    rounded_solvent_locations = np.around(final_solvent_peak_locations,
                                          decimals=2)  # [round(x, 2) for x in final_solvent_peak_locations]
    rounded_assigned_peaks_sorted_descending = np.around(assigned_peaks_sorted_descending, decimals=2)
    rounded_diff = [round(x, 3) for x in differences]

    return rounded_picked_locations_ppm, rounded_solvent_locations, rounded_assigned_peaks_sorted_descending, rounded_diff

def simulate_calc_data(spectral_xdata_ppm, calculated_locations, simulated_ydata):
    ###simulate calcutated data


    simulated_calc_ydata = np.zeros(len(spectral_xdata_ppm))

    for peak in calculated_locations:
        y = np.exp(-0.5 * ((spectral_xdata_ppm - peak) / 0.002) ** 2)
        simulated_calc_ydata += y

    scaling_factor = np.amax(simulated_ydata) / np.amax(simulated_calc_ydata)

    simulated_calc_ydata = simulated_calc_ydata*scaling_factor

    return simulated_calc_ydata

def lorentzian(p, w, p0, A):

    x = (p0 - p) / (w / 2)
    L = A / (1 + x ** 2)

    return L

def lorenz_curves(params, x, picked_points):

    y = np.zeros(len(x))
    for peak in picked_points:
        y += lorentzian(x, params['width' + str(peak)], params['pos' + str(peak)], params['amp' + str(peak)])
    return y

def gaussian(p,w,p0,A):

    y  = A *np.exp(-((p - p0)**2)/(2*(w)**2))

    return y

def minimisation(next_peak,fit_y, total_spectral_ydata,corr_distance):

    region = np.arange(max(0,next_peak - 100), min(next_peak +100,len(total_spectral_ydata)))

    params = Parameters()

    params.add('amp' + str(next_peak), value=total_spectral_ydata[next_peak],vary =  False,min = 0)
    params.add('width' + str(next_peak), value=4 * corr_distance, vary=True,min = 1*corr_distance,max = 8*corr_distance)
    params.add('pos' + str(next_peak), value=next_peak, vary=False)

        # print('minimising')

    out = Minimizer(residual, params,
                    fcn_args=(fit_y[region],next_peak,region, total_spectral_ydata[region]))

    results = out.minimize()

    # append the results params to the total params

    fit_yc =  lorentzian(np.arange(len(total_spectral_ydata)), results.params['width' + str(next_peak)],
                                 results.params['pos' + str(next_peak)], results.params['amp' + str(next_peak)]) + fit_y

    return fit_yc

def residual(params,fit_y,next_peak, x, y_data):

    y = lorentzian(x, params['width' + str(next_peak)], params['pos' + str(next_peak)], params['amp' + str(next_peak)]) + fit_y

    difference = abs(y - y_data)

    return difference

def iterative_peak_picking(total_spectral_ydata,threshold,corr_distance,):

    mu, std = norm.fit(total_spectral_ydata[0:1000])

    picked_peaks = []

    #find all maxima

    maxima = []

    for point in range(1,len(total_spectral_ydata)-1):
        if (total_spectral_ydata[point] > total_spectral_ydata[point+1]) & (total_spectral_ydata[point] > total_spectral_ydata[point-1]):
            maxima.append(point)

    #start fitting process

    fit_y = np.zeros(len(total_spectral_ydata))

    while len(maxima) > 0:

        params = Parameters()

        #find peak with greatest amplitude:

        ind1 = np.argmax(total_spectral_ydata[maxima])

        peak = maxima[ind1]

        picked_peaks.append(peak)

        fit_y = minimisation(peak, fit_y, total_spectral_ydata, corr_distance)

        new_maxima = []

        for ind2 in maxima:

            if total_spectral_ydata[ind2] > threshold*std + fit_y[ind2]:

                new_maxima.append(ind2)

        maxima = copy.copy(new_maxima)

    picked_peaks = sorted(picked_peaks)

    return picked_peaks,fit_y

def first_order_peak(start_ppm, J_vals, x_data, corr_distance, uc,m):

    #new first order peak generator using the method presented in Hoye paper

    start  = uc(str(start_ppm) + "ppm")

    start_Hz = uc.hz(start)

    J_vals = np.array(J_vals)

    if len(J_vals) > 0:

        peaks = np.zeros((2*m +1)**len(J_vals))

        if m == 0.5:

            l =[1,-1]

        if m == 1:

            l= [1,0,-1]

        #signvector generator

        signvectors = itertools.product(l, repeat=len(J_vals))

        shifts = []

        for ind,sv in enumerate(signvectors):
            shift = J_vals * sv
            shift = start_Hz + 0.5 * (np.sum(shift))
            peaks[ind] = shift

        peaks = np.sort(peaks)

        peak_vector = np.array(sorted(list(set(peaks)),reverse = True))

        amp_vector = np.zeros(len(peak_vector))

        for peak in peaks:
            index = np.where(peak_vector == peak)
            amp_vector[index] += 1

        pv = []

        for index,peak in enumerate(peak_vector):

            pv.append(uc(peak,"Hz"))

        peak_vector = pv

    else:
        peak_vector = start

    split_params = Parameters()

    for index, peak in enumerate(peak_vector):
        split_params.add('amp' + str(peak), value=amp_vector[index])
        split_params.add('pos' + str(peak), value=peak)
        split_params.add('width' + str(peak), value=2 * corr_distance)

    y = lorenz_curves(split_params, x_data, peak_vector)

    y= y/np.max(y)

    return split_params,peak_vector,amp_vector,y

def solvent_removal(simulated_y_data,spectral_xdata_ppm,solvent,uc,picked_peaks):

    if solvent == 'chloroform':

        exp_ppm = [77]

        Jv = [[64]]

    elif solvent == 'dimethylsulfoxide':

        exp_ppm = [39.51]

        Jv = [[42,42,42]]

    elif solvent == 'pyridine':

        exp_ppm = [150.35,135.91,123.87]


        Jv= [[55,55],[49,49],[50,50]]

    elif solvent == 'methanol':

        exp_ppm = [49.15]

        Jv = [[42.8,42.8,42.8]]

    elif solvent == 'benzene':

        exp_ppm = [128.39]

        Jv = [[24.3]]

    else:

        exp_ppm = []
        Jv =[[]]

    #remove all solvent peaks

    removed = []

    for J,p in zip(Jv,exp_ppm):

        exp = uc(p,"ppm")

        region = np.arange(exp - 1000, exp+1000)

        peak_region = []

        for peak in picked_peaks:

            if (peak > exp -1000) & (peak < exp+1000):
                peak_region.append(peak)

        #simulate solvent curve

        #find peak centre

        if region[0] + region[-1] & 1:
            centre = (int((region[0] + region[-1] + 1) / 2))
        else:
            centre = (int((region[0] + region[-1]) / 2))

        centre = uc.ppm(centre)

        params,peak_vector, amp_vector, y = first_order_peak(centre, J,np.array(region), 1, uc,1)

        #use simulated curve in convolution

        convolved_y = np.convolve(simulated_y_data[region] , y,'same')

        mxpoint = np.argmax(convolved_y)

        mxppm = uc.ppm(region[mxpoint])

        #simulate peak in new position

        params,fit_s_peaks, amp_vector, fit_s_y = first_order_peak(mxppm, J, np.array(region), 1, uc,1)



        #find average of fitted peaks for referencing:

        av = sum(fit_s_peaks)/len(fit_s_peaks)

        avppm = uc.ppm(av)

        spectral_xdata_ppm -= avppm - p

        to_remove = []

        # find picked peaks closest to the "fitted" solvent multiplet

        for peak in fit_s_peaks:

            i = np.abs(np.array(picked_peaks) - peak).argmin()

            to_remove.append(i)

        removed.extend([picked_peaks[i] for i in to_remove])

        to_remove = sorted(list(set(to_remove)),reverse=True)

        for peak in to_remove:
            picked_peaks.pop(peak)

    removed = np.array(removed)

    return picked_peaks,removed