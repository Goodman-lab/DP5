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
import os
import numpy as np

# Standard DP4 parameters
meanC = 0.0
meanH = 0.0
stdevC = 2.269372270818724
stdevH = 0.18731058105269952


class DP4data:
    def __init__(self):
        self.Cshifts = []  # Carbon shifts used in DP4 calculation
        self.Cexp = []  # Carbon experimental shifts used in DP4 calculation
        self.Clabels = []  # Carbon atom labels
        self.Hshifts = []  # Proton shifts used in DP4 calculation
        self.Hexp = []  # Proton experimental shifts used in DP4 calculation
        self.Hlabels = []  # Proton atom labels
        self.Cscaled = []  # Internally scaled carbon shifts
        self.Hscaled = []  # Internally scaled proton shifts
        self.Cerrors = []  # Scaled Carbon prediction errors
        self.Herrors = []  # Scaled Proton prediction errors
        self.Cprobs = []  # Scaled carbon prediction error probabilities
        self.Hprobs = []  # Scaled proton prediction error probabilities
        self.CDP4probs = []  # Carbon DP4 probabilities
        self.HDP4probs = []  # Proton DP4 probabilities
        self.DP4probs = []  # combined Carbon and Proton DP4 probabilities
        self.output = str()  # final DP4 output


def ProcessIsomers(DP4data, Isomers):
    # extract calculated and experimental shifts and add to DP4data instance

    # Carbon

    # make sure any shifts with missing peaks are removed from all isomers

    removedC = []

    removedH = []

    for iso in Isomers:

        DP4data.Cshifts.append([])
        DP4data.Cexp.append([])
        DP4data.Clabels.append([])

        for shift, exp, label in zip(iso.Cshifts, iso.Cexp, iso.Clabels):

            if exp != '':
                DP4data.Cshifts[-1].append(shift)
                DP4data.Cexp[-1].append(exp)
                DP4data.Clabels[-1].append(label)

            elif label not in removedC:

                removedC.append(label)

    for l in removedC:

        for j, Clabel in enumerate(DP4data.Clabels):

            if l in Clabel:
                i = Clabel.index(l)

                DP4data.Cshifts[j].pop(i)

                DP4data.Cexp[j].pop(i)

                DP4data.Clabels[j].pop(i)

    # proton
    for iso in Isomers:

        DP4data.Hshifts.append([])
        DP4data.Hexp.append([])
        DP4data.Hlabels.append([])

        for shift, exp, label in zip(iso.Hshifts, iso.Hexp, iso.Hlabels):

            if exp != '':
                DP4data.Hshifts[-1].append(shift)
                DP4data.Hexp[-1].append(exp)
                DP4data.Hlabels[-1].append(label)

            elif label not in removedH:

                removedH.append(label)

    for l in removedH:

        for j, Hlabel in enumerate(DP4data.Hlabels):

            if l in Hlabel:
                i = Hlabel.index(l)

                DP4data.Hshifts[j].pop(i)

                DP4data.Hexp[j].pop(i)

                DP4data.Hlabels[j].pop(i)

    return DP4data


def InternalScaling(DP4data):
    # perform internal scaling process

    # calculate prediction errors

    if len(DP4data.Cexp[0]) > 1:

        for Cshifts, Cexp in zip(DP4data.Cshifts, DP4data.Cexp):
            DP4data.Cscaled.append(ScaleNMR(Cshifts, Cexp))

        for Cscaled, Cexp in zip(DP4data.Cscaled, DP4data.Cexp):
            DP4data.Cerrors.append([Cscaled[i] - Cexp[i] for i in range(0, len(Cscaled))])

    if len(DP4data.Hexp[0]) > 1:

        for Hshifts, Hexp in zip(DP4data.Hshifts, DP4data.Hexp):
            DP4data.Hscaled.append(ScaleNMR(Hshifts, Hexp))

        for Hscaled, Hexp in zip(DP4data.Hscaled, DP4data.Hexp):
            DP4data.Herrors.append([Hscaled[i] - Hexp[i] for i in range(0, len(Hscaled))])

    return DP4data


def ScaleNMR(calcShifts, expShifts):
    slope, intercept, r_value, p_value, std_err = stats.linregress(expShifts,
                                                                   calcShifts)
    scaled = [(x - intercept) / slope for x in calcShifts]

    return scaled


def CalcProbs(DP4data, Settings):
    # calculates probability values for each scaled prediction error value using the chosen statistical model

    if Settings.StatsModel == 'g' or 'm':

        print(Settings.StatsParamFile)

        if Settings.StatsParamFile == 'none':

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

            for errors in DP4data.Herrors:
                DP4data.Hprobs.append([MultiGausProbability(e, Hmeans, Hstdevs) for e in errors])

    return DP4data


def SingleGausProbability(error, mean, stdev):
    z = abs((error - mean) / stdev)
    cdp4 = 2 * stats.norm.cdf(-z)

    return cdp4


def MultiGausProbability(error, means, stdevs):
    res = 0

    for mean, stdev in zip(means, stdevs):
        res += stats.norm(mean, stdev).pdf(error)

    return res / len(means)


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


def CalcDP4(DP4data):
    # Calculate Carbon DP4 probabilities

    for probs in DP4data.Cprobs:

        DP4data.CDP4probs.append(1)

        for p in probs:
            DP4data.CDP4probs[-1] *= p

    # Calculate Proton DP4 probabilities

    for probs in DP4data.Hprobs:

        DP4data.HDP4probs.append(1)

        for p in probs:
            DP4data.HDP4probs[-1] *= p

    # Calculate Combined DP4 probabilities

    for Hp, Cp in zip(DP4data.HDP4probs, DP4data.CDP4probs):
        DP4data.DP4probs.append(Hp * Cp)

    Cs = sum(DP4data.CDP4probs)

    Hs = sum(DP4data.HDP4probs)

    Ts = sum(DP4data.DP4probs)

    DP4data.CDP4probs = [i / Cs for i in DP4data.CDP4probs]

    DP4data.HDP4probs = [i / Hs for i in DP4data.HDP4probs]

    DP4data.DP4probs = [i / Ts for i in DP4data.DP4probs]

    return DP4data


def PrintAssignment(DP4Data):
    isomer = 0

    for Clabels, Cshifts, Cexp, Cscaled in zip(DP4Data.Clabels, DP4Data.Cshifts, DP4Data.Cexp, DP4Data.Cscaled):
        DP4Data.output += ("\n\nAssigned C shifts for isomer " + str(isomer + 1) + ": ")

        PrintNMR(Clabels, Cshifts, Cscaled, Cexp, DP4Data)

        isomer += 1

    isomer = 0

    for Hlabels, Hshifts, Hexp, Hscaled in zip(DP4Data.Hlabels, DP4Data.Hshifts, DP4Data.Hexp, DP4Data.Hscaled):
        DP4Data.output += ("\n\nAssigned H shifts for isomer " + str(isomer + 1) + ": ")

        PrintNMR(Hlabels, Hshifts, Hscaled, Hexp, DP4Data)

        isomer += 1


def PrintNMR(labels, values, scaled, exp, DP4Data):
    s = np.argsort(values)

    svalues = np.array(values)[s]

    slabels = np.array(labels)[s]

    sscaled = np.array(scaled)[s]

    sexp = np.array(exp)[s]

    DP4Data.output += ("\nlabel, calc, corrected, exp, error")

    for i in range(len(labels)):
        DP4Data.output += ("\n" + format(slabels[i], "6s") + ' ' + format(svalues[i], "6.2f") + ' '
                           + format(sscaled[i], "6.2f") + ' ' + format(sexp[i], "6.2f") + ' ' +
                           format(sexp[i] - sscaled[i], "6.2f"))


def MakeOutput(DP4Data, Isomers, Settings):
    # add some info about the calculation

    DP4Data.output += Settings.InputFiles[0] + "\n"

    DP4Data.output += "\n" + "Solvent = " + Settings.Solvent

    DP4Data.output += "\n" + "Force Field = " + Settings.ForceField + "\n"

    if 'o' in Settings.Workflow:
        DP4Data.output += "\n" + "DFT optimisation Functional = " + Settings.oFunctional
        DP4Data.output += "\n" + "DFT optimisation Basis = " + Settings.oBasisSet

    if 'e' in Settings.Workflow:
        DP4Data.output += "\n" + "DFT energy Functional = " + Settings.eFunctional
        DP4Data.output += "\n" + "DFT energy Basis = " + Settings.eBasisSet

    if 'n' in Settings.Workflow:
        DP4Data.output += "\n" + "DFT NMR Functional = " + Settings.nFunctional
        DP4Data.output += "\n" + "DFT NMR Basis = " + Settings.nBasisSet

    if Settings.StatsParamFile != "none":
        DP4Data.output += "\n\nStats model = " + Settings.StatsParamFile

    DP4Data.output += "\n\nNumber of isomers = " + str(len(Isomers))

    c = 1

    for i in Isomers:
        DP4Data.output += "\nNumber of conformers for isomer " + str(c) + " = " + str(len(i.Conformers))

        c += 1

    PrintAssignment(DP4Data)

    DP4Data.output += ("\n\nResults of DP4 using Proton: ")

    for i, p in enumerate(DP4Data.HDP4probs):
        DP4Data.output += ("\nIsomer " + str(i + 1) + ": " + format(p * 100, "4.1f") + "%")

    DP4Data.output += ("\n\nResults of DP4 using Carbon: ")

    for i, p in enumerate(DP4Data.CDP4probs):
        DP4Data.output += ("\nIsomer " + str(i + 1) + ": " + format(p * 100, "4.1f") + "%")

    DP4Data.output += ("\n\nResults of DP4: ")

    for i, p in enumerate(DP4Data.DP4probs):
        DP4Data.output += ("\nIsomer " + str(i + 1) + ": " + format(p * 100, "4.1f") + "%")

    print("number of c protons = " + str(len(Isomers[0].Hlabels)))
    print("number of c carbons = " + str(len(Isomers[0].Clabels)))

    print("number of e protons = " + str(len(DP4Data.Hexp[0])))
    print("number of e carbons = " + str(len(DP4Data.Cexp[0])))

    print(DP4Data.output)

    if Settings.OutputFolder == '':

        out = open(str(os.getcwd()) + "/" + str(Settings.InputFiles[0] + "NMR.dp4"), "w+")

    else:

        out = open(os.path.join(Settings.OutputFolder, str(Settings.InputFiles[0] + "NMR.dp4")), "w+")

    out.write(DP4Data.output)

    out.close()

    return DP4Data