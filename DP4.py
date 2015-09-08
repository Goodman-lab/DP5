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
import pickle
import re
import bisect

meanC = 0.0
meanH = 0.0
stdevC = 2.269372270818724
stdevH = 0.18731058105269952
output = []


def DP4(Clabels, Cvalues, Hlabels, Hvalues, Cexp, Hexp, settings):

    Print(str(Cexp))
    Print(str(Hexp))

    C_cdp4 = []
    H_cdp4 = []
    Comb_cdp4 = []

    for isomer in range(0, len(Cvalues)):

        sortedClabels, sortedCvalues, sortedCexp, _ =\
            AssignExpNMR(Clabels, Cvalues[isomer], Cexp)
        scaledC = ScaleNMR(sortedCvalues, sortedCexp)

        sortedHlabels, sortedHvalues, sortedHexp, Jvalues = \
            AssignExpNMR(Hlabels, Hvalues[isomer], Hexp)
        scaledH = ScaleNMR(sortedHvalues, sortedHexp)

        ScaledErrorsC = [sortedCexp[i] - scaledC[i]
                         for i in range(0, len(scaledC))]
        ScaledErrorsH = [sortedHexp[i] - scaledH[i]
                         for i in range(0, len(scaledH))]

        Print("\nAssigned shifts for isomer " + str(isomer+1) + ": ")
        PrintNMR('C', sortedClabels, sortedCvalues, scaledC, sortedCexp)
        Print("Max C error: " + format(max(ScaledErrorsC, key=abs), "6.2f"))
        PrintNMR('H', sortedHlabels, sortedHvalues, scaledH, sortedHexp)
        Print("Max H error: " + format(max(ScaledErrorsH, key=abs), "6.2f"))

        if settings.PDP4:
            C_cdp4.append(CalculateCDP4(ScaledErrorsC, meanC, stdevC))
            H_cdp4.append(CalculateCDP4(ScaledErrorsH, meanH, stdevH))
        elif settings.EP5:
            C_cdp4.append(CalculatePDP4(ScaledErrorsC, meanC, stdevC))
            #C_cdp4.append(CalculateKDE(ScaledErrorsC, settings.ScriptDir + '/NucCErr.pkl'))
            H_cdp4.append(CalculatePDP4(ScaledErrorsH, meanH, stdevH))
            #H_cdp4.append(CalculateKDE(ScaledErrorsH, settings.ScriptDir + '/NucHErr.pkl'))
        Comb_cdp4.append(C_cdp4[-1]*H_cdp4[-1])
        Print("\nDP4 based on C: " + format(C_cdp4[-1], "6.2e"))
        Print("DP4 based on H: " + format(H_cdp4[-1], "6.2e"))

    relCDP4 = [(100*x)/sum(C_cdp4) for x in C_cdp4]
    relHDP4 = [(100*x)/sum(H_cdp4) for x in H_cdp4]
    relCombDP4 = [(100*x)/sum(Comb_cdp4) for x in Comb_cdp4]

    PrintRelDP4('both carbon and proton data', relCombDP4)
    PrintRelDP4('carbon data only', relCDP4)
    PrintRelDP4('proton data only', relHDP4)
    
    return output


def AssignExpNMR(labels, calcShifts, exp):

    #Replace all 'or' and 'OR' with ',', remove all spaces and 'any'
    exp = re.sub(r"or|OR", ',', exp, flags=re.DOTALL)
    exp = re.sub(r" |any", '', exp, flags=re.DOTALL)

    #Get all assignments, split mulitassignments
    ExpLabels = re.findall(r"(?<=\().*?(?=\)|;)", exp, flags=re.DOTALL)
    ExpLabels = [x.split(',') for x in ExpLabels]
    
    #Get J value data (if any)
    Jdata = re.findall(r"(?<=\().*?(?=\))", exp, flags=re.DOTALL)
    Jvalues = []
    for d in Jdata:
        if ';' in d:
            Jvalues.append([float(x) for x in (d.split(';')[1]).split(',')])
        else:
            Jvalues.append([0])
    
    print Jvalues
    #Remove assignments and get shifts
    ShiftData = (re.sub(r"\(.*?\)", "", exp, flags=re.DOTALL)).split(',')
    expShifts = [float(x) for x in ShiftData]

    #Prepare sorted calculated data with labels and sorted exp data
    sortedCalc = sorted(calcShifts)
    sortedExp = sorted(expShifts)
    sortedExpLabels = SortExpAssignments(expShifts, ExpLabels)
    sortedCalcLabels = []
    for v in sortedCalc:
        index = calcShifts.index(v)
        sortedCalcLabels.append(labels[index])

    assignedExpLabels = ['' for i in range(0, len(sortedExp))]

    #First pass - assign the unambiguous shifts
    for v in range(0, len(sortedExp)):
        if len(sortedExpLabels[v]) == 1 and sortedExpLabels[v][0] != '':
            #Check that assignment exists in computational data
            if sortedExpLabels[v][0] in labels:
                assignedExpLabels[v] = sortedExpLabels[v][0]
            else:
                Print("Label " + sortedExpLabels[v][0] +
                " not found in among computed shifts, please check NMR assignment.")
                quit()
    #Second pass - assign shifts from a limited set
    for v in range(0, len(sortedExp)):
        if len(sortedExpLabels[v]) != 1 and sortedExpLabels[v][0] != '':
            for l in sortedCalcLabels:
                if l in sortedExpLabels[v] and l not in assignedExpLabels:
                    assignedExpLabels[v] = l
                    break
    #Final pass - assign unassigned shifts in order
    for v in range(0, len(sortedExp)):
        if sortedExpLabels[v][0] == '':
            for l in sortedCalcLabels:  # Take the first free label
                if l not in assignedExpLabels:
                    assignedExpLabels[v] = l
                    break
    sortedCalc = []
    #Rearrange calc values to match the assigned labels
    for l in assignedExpLabels:
        index = labels.index(l)
        sortedCalc.append(calcShifts[index])

    return assignedExpLabels, sortedCalc, sortedExp, Jvalues


def SortExpAssignments(shifts, assignments):
    tempshifts = list(shifts)
    tempassignments = list(assignments)
    sortedassignments = []
    while len(tempassignments) > 0:
        index = tempshifts.index(min(tempshifts))
        sortedassignments.append(tempassignments[index])
        tempshifts.pop(index)
        tempassignments.pop(index)
    return sortedassignments


#Scale the NMR shifts
def ScaleNMR(calcShifts, expShifts):
    slope, intercept, r_value, p_value, std_err = stats.linregress(expShifts,
                                                                   calcShifts)
    scaled = [(x-intercept)/slope for x in calcShifts]
    return scaled


def PrintNMR(nucleus, labels, values, scaled, exp):
    Print("\nAssigned " + nucleus +
          " shifts: (label, calc, corrected, exp, error)")
    for i in range(0, len(labels)):
        Print(format(labels[i], "6s") + ' ' + format(values[i], "6.2f") + ' '
            + format(scaled[i], "6.2f") + ' ' + format(exp[i], "6.2f") + ' ' +
            format(exp[i]-scaled[i], "6.2f"))


def PrintRelDP4(title, RelDP4):
    Print("\nResults of DP4 using " + title + ":")
    for i in range(0, len(RelDP4)):
        Print("Isomer " + str(i+1) + ": " + format(RelDP4[i], "4.1f") + "%")


def Print(s):
    print s
    output.append(s)


def CalculateCDP4(errors, expect, stdev):
    cdp4 = 1.0
    for e in errors:
        z = abs((e-expect)/stdev)
        cdp4 = cdp4*2*stats.norm.cdf(-z)
    return cdp4


#Alternative function using probability density function instead of cdf
def CalculatePDP4(errors, expect, stdev):
    pdp4 = 1.0
    for e in errors:
        #z = (e-expect)/stdev
        #pdp4 = pdp4*stats.norm.pdf(z)
        pdp4 = pdp4*stats.norm(expect, stdev).pdf(e)
    return pdp4


#use as CalculateKDE(errors, 'NucCErr.pkl') for C or
#CalculateKDE(errors, 'NucHErr.pkl') for H
#load empirical error data from file and use KDE to construct pdf
def CalculateKDE(errors, PickleFile):

    pkl_file = open(PickleFile, 'rb')
    ErrorData = pickle.load(pkl_file)
    kde = stats.gaussian_kde(ErrorData)

    ep5 = 1.0
    for e in errors:
        #z = (e-expect)/stdev
        #pdp4 = pdp4*stats.norm.pdf(z)
        ep5 = ep5*float(kde(e)[0])
    return ep5


#use as CalculateRKDE(errors, 'RKDEC.pkl') for C or
#CalculateKDE(errors, 'RKDEH.pkl') for H
#load empirical error data from file and use KDE to construct several pdfs,
#one for each chemical shift region
def CalculateRKDE(errors, shifts, PickleFile):
    
    #Load the data
    pkl_file = open(PickleFile, 'rb')
    regions, RErrors = pickle.load(pkl_file)
    
    #Reconstruct the distributions for each region
    kdes = []
    for es in RErrors:
        kdes.append(stats.gaussian_kde(es))

    ep5 = 1.0
    for i, e in enumerate(errors):
        region = bisect.bisect_left(regions, shifts[i])
        ep5 = ep5*float((kdes[region])(e)[0])
    return ep5
