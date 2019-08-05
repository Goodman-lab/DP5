# -*- coding: utf-8 -*-
"""
Created on Wed May 27 14:18:37 2015
Updated on July 30 14:18:37 2015
@author: ke291

Equivalent and compact port of DP4.jar to python. The results
produced are essentially equivalent, but not identical due to different
floating point precision used in the Python (53 bits) and Java (32 bits)
implementation.
"""
from scipy import stats
import bisect

# Standard DP4 parameters
meanC = 0.0
meanH = 0.0
stdevC = 2.269372270818724
stdevH = 0.18731058105269952


class DP4data:
    def __init__(self, InputPath):
        self.Cshifts = []  # Carbon shifts used in DP4 calculation
        self.Cexp = []  # Carbon experimental shifts used in DP4 calculation
        self.Hshifts = []  # Proton shifts used in DP4 calculation
        self.Hexp = []  # Proton experimental shifts used in DP4 calculation
        self.Cscaled = []  # Internally scaled carbon shifts
        self.Hscaled = []  # Internally scaled proton shifts
        self.Cerrors = []  # Scaled C prediction errors
        self.Herrors = []  # Scaled H prediction errors
        self.Cprobs = []   # Scaled prediction error probabilities
        self.Hprobs = []
        self.CDP4probs = []
        self.HDP4probs = []  # Isomer DP4 probabilities


def ProcessIsomers(DP4data, Isomers):
    # extract calculated and experimental shifts and add to DP4data instance

    # Carbon
    for iso in Isomers:

        DP4data.Cshifts.append([])
        DP4data.Cexp.append([])

        for shift, exp in zip(iso.Cshifts, iso.Cexp):
            if exp != '':
                DP4data.Cshifts[-1].append(shift)
                DP4data.Cexp[-1].append(exp)
    # proton
    for iso in Isomers:

        DP4data.Hshifts.append([])
        DP4data.Hexp.append([])

        for shift, exp in zip(iso.Hshifts, iso.Hexp):
            if exp != '':
                DP4data.Hshifts[-1].append(shift)
                DP4data.Hexp[-1].append(exp)

    return DP4data


def InternalScaling(DP4data):
    # perform internal scaling process

    for Cshifts, Cexp in zip(DP4data.Cshifts, DP4data.Cexp):
        DP4data.Cscaled.append(ScaleNMR(Cshifts, Cexp))

    for Hshifts, Hexp in zip(DP4data.Hshifts, DP4data.Hexp):
        DP4data.Hscaled.append(ScaleNMR(Hshifts, Hexp))

    # calculate prediction errors

    for Cscaled, Cexp in zip(DP4data.Cscaled, DP4data.Cexp):
        DP4data.Cerrors.append([Cscaled[i] - Cexp[i] for i in range(0, len(Cscaled))])

    for Hscaled, Hexp in zip(DP4data.Hscaled, DP4data.Hexp):
        DP4data.Herrors.append([Hscaled[i] - Hexp[i] for i in range(0, len(Hscaled))])

    return DP4data


def ScaleNMR(calcShifts, expShifts):
    slope, intercept, r_value, p_value, std_err = stats.linregress(expShifts,
                                                                   calcShifts)
    scaled = [(x - intercept) / slope for x in calcShifts]

    return scaled


def CalcProbs(DP4data, Settings):

    #calculates probability values for each scaled prediction error value using the chosen statistical model

    if Settings.StatsModel == 'g':

        if Settings.StatsParamFile == '':

            print('No stats model provided, using default')

            for errors in DP4data.Cerrors:
                DP4data.Cprobs.append([SingleGausProbability(e, meanC, stdevC) for e in errors])

            for errors in DP4data.Herrors:
                DP4data.Hprobs.append([SingleGausProbability(e, meanH, stdevH) for e in errors])

        else:

            print('Using stats model provided')

            Cmeans, Cstdevs, Hmeans, Hstdevs = ReadParamFile(Settings.StatsParamFile, Settings.StatsModel)

            for errors in DP4data.Cerrors:
                DP4data.Cprobs.append([MultiGausProbability(e, Cmeans, Cstdevs) for e in errors])

            for errors in DP4data.Cerrors:
                DP4data.Cprobs.append([MultiGausProbability(e, Hmeans, Hstdevs) for e in errors])

    return DP4data


def SingleGausProbability(error, mean, stdev):

    z = abs((error - mean) / stdev)
    cdp4 = 2 * stats.norm.cdf(-z)

    return cdp4


def MultiGausProbability(error, means, stdevs):

    prob = MultiGaus(means, stdevs, error)

    return prob


def MultiGaus(means, stdevs, x):
    res = 0

    for mean, stdev in zip(means, stdevs):

        res += stats.norm(mean, stdev).pdf(x)

    return res / len(means)


def ReadParamFile(f, t):
    infile = open(f, 'r')
    inp = infile.readlines()
    infile.close()

    if t not in inp[0]:
        print("Wrong parameter file type, exiting...")
        quit()

    if t == 'g':
        Cmeans = [float(x) for x in inp[1].split(',')]
        Cstdevs = [float(x) for x in inp[2].split(',')]
        Hmeans = [float(x) for x in inp[3].split(',')]
        Hstdevs = [float(x) for x in inp[4].split(',')]

        return Cmeans, Cstdevs, Hmeans, Hstdevs


def CalcDP4(DP4data):

    for probs in DP4data.Cprobs:

        DP4data.CDP4probs.append(1)

        for p in probs:

            DP4data.CDP4probs[-1] *= p

    for probs in DP4data.Hprobs:

        DP4data.HDP4probs.append(1)

        for p in probs:

            DP4data.HDP4probs[-1] *= p

    Cs = sum(DP4data.CDP4probs)

    Hs = sum(DP4data.HDP4probs)

    DP4data.CDP4probs = [i / Cs for i in DP4data.CDP4probs]

    DP4data.HDP4probs = [i / Hs for i in DP4data.HDP4probs]

    return DP4data