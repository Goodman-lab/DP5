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
import re
import bisect

#Standard DP4 parameters
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

        ScaledErrorsC = [scaledC[i] - sortedCexp[i]
                         for i in range(0, len(scaledC))]
        ErrorsC = [sortedCvalues[i] - sortedCexp[i]
                         for i in range(0, len(sortedCvalues))]
        ScaledErrorsH = [scaledH[i] - sortedHexp[i]
                         for i in range(0, len(scaledH))]
        ErrorsH = [sortedHvalues[i] - sortedHexp[i]
                         for i in range(0, len(sortedHvalues))]

        Print("\nAssigned shifts for isomer " + str(isomer+1) + ": ")
        PrintNMR('C', sortedClabels, sortedCvalues, scaledC, sortedCexp)
        Print("Max C error: " + format(max(ScaledErrorsC, key=abs), "6.2f"))
        PrintNMR('H', sortedHlabels, sortedHvalues, scaledH, sortedHexp)
        Print("Max H error: " + format(max(ScaledErrorsH, key=abs), "6.2f"))
        
        Cprob, Hprob = CalcProbabilities(ScaledErrorsC, ScaledErrorsH, ErrorsC,
                                         ErrorsH, sortedCexp, sortedHexp, settings)
        C_cdp4.append(Cprob)
        H_cdp4.append(Hprob)

        Comb_cdp4.append(C_cdp4[-1]*H_cdp4[-1])
        Print("\nDP4 based on C: " + format(C_cdp4[-1], "6.2e"))
        Print("DP4 based on H: " + format(H_cdp4[-1], "6.2e"))

    relCDP4 = [(100*x)/sum(C_cdp4) for x in C_cdp4]
    relHDP4 = [(100*x)/sum(H_cdp4) for x in H_cdp4]
    relCombDP4 = [(100*x)/sum(Comb_cdp4) for x in Comb_cdp4]

    PrintRelDP4('all available data', relCombDP4)
    PrintRelDP4('carbon data only', relCDP4)
    PrintRelDP4('proton data only', relHDP4)
    
    return output


def f3(deltaExp, deltaCalc):
    if (deltaCalc/deltaExp) > 1.0:
        return deltaExp**3/deltaCalc
    else:
        return deltaExp*deltaCalc


def CalcSpreads(AllErrors):
    spreads = []
    for i in range(len(AllErrors[0])):
        dserrors = [AllErrors[x][i] for x in range(len(AllErrors))]
        spreads.append(max(dserrors)-min(dserrors))
    return spreads


def ReorderAllErrors(AllErrors, AllLabels):
    
    Labels1 = AllLabels[0]
    AllErrorsR = []
    
    for i, (labels, errors) in enumerate(zip(AllLabels, AllErrors)):
        if i ==0:
            AllErrorsR.append(errors)
        else:
            AllErrorsR.append(ReorderData(Labels1, labels, errors))
    
    return Labels1, AllErrorsR


#given 2 sets of labels and one set of values, the set of values will be
#reordered from the order in labels2 to the order in labels1
def ReorderData(labels1, labels2, data2):
    
    newData = []
    templabels2 = list(labels2)
    
    for l in labels1:
        index = templabels2.index(l)
        newData.append(data2[index])
        #Remove the value to avoid adding the same value twice
        templabels2[index] = 'XX'
    
    return newData


def AssignExpNMR(labels, calcShifts, exp):

    expLabels, expShifts, expJs = ReadExp(exp)
    #Prepare sorted calculated data with labels and sorted exp data
    sortedCalc = sorted(calcShifts)
    sortedExp = sorted(expShifts)
    sortedExpLabels = SortExpAssignments(expShifts, expLabels)
    sortedCalcLabels = []
    sortedJvalues = []
    for v in sortedExp:
        index = expShifts.index(v)
        sortedJvalues.append(expJs[index])
    
    tempCalcShifts = list(calcShifts)
    for v in sortedCalc:
        index = tempCalcShifts.index(v)
        sortedCalcLabels.append(labels[index])
        #Remove the value to avoid duplicate labels
        tempCalcShifts[index] = -100

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

    return assignedExpLabels, sortedCalc, sortedExp, sortedJvalues


def ReadExp(exp):
    
    #Replace all 'or' and 'OR' with ',', remove all spaces and 'any'
    texp = re.sub(r"or|OR", ',', exp, flags=re.DOTALL)
    texp = re.sub(r" |any", '', texp, flags=re.DOTALL)

    #Get all assignments, split mulitassignments
    expLabels = re.findall(r"(?<=\().*?(?=\)|;)", texp, flags=re.DOTALL)
    expLabels = [x.split(',') for x in expLabels]
    
    #Get J value data (if any)
    Jdata = re.findall(r"(?<=\().*?(?=\))", texp, flags=re.DOTALL)
    expJs = []
    
    for d in Jdata:
        if ';' in d:
            expJs.append([float(x) for x in (d.split(';')[1]).split(',')])
        else:
            expJs.append([0.0])
    
    print(expJs)
    #Remove assignments and get shifts
    ShiftData = (re.sub(r"\(.*?\)", "", exp, flags=re.DOTALL)).split(',')
    expShifts = [float(x) for x in ShiftData]
    
    return expLabels, expShifts, expJs


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


def CalcProbabilities(SErrorsC, SErrorsH, ErrorsC, ErrorsH, Cexp, Hexp, settings):
    
    if settings.StatsModel == 'g':
        if settings.StatsParamFile == '':
            C_cdp4 = CalculateCDP4(SErrorsC, meanC, stdevC)
            H_cdp4 = CalculateCDP4(SErrorsH, meanH, stdevH)
        else:
            Cmean, Cstdev, Hmean, Hstdev =\
                ReadParamFile(settings.StatsParamFile, 'g')
            C_cdp4 = CalculatePDP4(SErrorsC, Cmean, Cstdev)
            H_cdp4 = CalculatePDP4(SErrorsH, Hmean, Hstdev)
            
    elif settings.StatsModel == 'm':
        C_cdp4 = CalculateMultiGaus(SErrorsC, settings.StatsParamFile,
                                         settings.StatsModel, 'C')
        H_cdp4 = CalculateMultiGaus(SErrorsH, settings.StatsParamFile,
                                         settings.StatsModel, 'H')
    elif settings.StatsModel == 'n':
        C_cdp4 = CalculateRMultiGaus(SErrorsC, Cexp, settings.StatsParamFile,
                                         settings.StatsModel, 'C')
        H_cdp4 = CalculateRMultiGaus(SErrorsH, Hexp, settings.StatsParamFile,
                                         settings.StatsModel, 'H')    
    elif settings.StatsModel == 'p':
        C_cdp4 = CalculateRMultiGaus(ErrorsC, Cexp, settings.StatsParamFile,
                                         settings.StatsModel, 'C')
        H_cdp4 = CalculateRMultiGaus(ErrorsH, Hexp, settings.StatsParamFile,
                                         settings.StatsModel, 'H')    
    else:
        print("Invalid stats model")
        C_cdp4 = 0.0
        H_cdp4 = 0.0
        
    return C_cdp4, H_cdp4


def PrintNMR(nucleus, labels, values, scaled, exp):
    Print("\nAssigned " + nucleus +
          " shifts: (label, calc, corrected, exp, error)")
    for i in range(len(labels)):
        Print(format(labels[i], "6s") + ' ' + format(values[i], "6.2f") + ' '
            + format(scaled[i], "6.2f") + ' ' + format(exp[i], "6.2f") + ' ' +
            format(exp[i]-scaled[i], "6.2f"))


def PrintErrors(labels, AllErrors):
    
    Print("\nLabel " + ' '.join(['DS' + str(x+1) for x in range(len(AllErrors))]))
    for i in range(len(labels)):
        Print(labels[i] +
        ' '.join([format(AllErrors[x][i], "6.2f") for x in range(len(AllErrors))]))


def PrintRelDP4(title, RelDP4):
    Print("\nResults of DP4 using " + title + ":")
    for i in range(len(RelDP4)):
        Print("Isomer " + str(i+1) + ": " + format(RelDP4[i], "4.1f") + "%")


def Print(s):
    print(s)
    output.append(s)


def ReadRMultiGausFile(f, t):
    
    infile = open(f, 'r')
    inp = infile.readlines()
    infile.close()
    
    Cregions = [float(x) for x in inp[1].split(',')]
    Cmeans = [[float(x) for x in y.split(',')] for y in inp[2].split(';')]
    Cstdevs = [[float(x) for x in y.split(',')] for y in inp[3].split(';')]
    
    Hregions = [float(x) for x in inp[4].split(',')]
    Hmeans = [[float(x) for x in y.split(',')] for y in inp[5].split(';')]
    Hstdevs = [[float(x) for x in y.split(',')] for y in inp[6].split(';')]

    return Cregions, Cmeans, Cstdevs, Hregions, Hmeans, Hstdevs


def ReadParamFile(f, t):
    
    infile = open(f, 'r')
    inp = infile.readlines()
    infile.close()
    
    if t not in inp[0]:
        print("Wrong parameter file type, exiting...")
        quit()
    
    if t == 'g':
        [Cmean, Cstdev] = [float(x) for x in inp[1].split(',')]
        [Hmean, Hstdev] = [float(x) for x in inp[2].split(',')]
        return Cmean, Cstdev, Hmean, Hstdev
    
    if t == 'm':
        Cmeans = [float(x) for x in inp[1].split(',')]
        Cstdevs = [float(x) for x in inp[2].split(',')]
        Hmeans = [float(x) for x in inp[3].split(',')]
        Hstdevs = [float(x) for x in inp[4].split(',')]
        return Cmeans, Cstdevs, Hmeans, Hstdevs
    
    elif t == 'k':
        Cerrors = [float(x) for x in inp[1].split(',')]
        Herrors = [float(x) for x in inp[2].split(',')]
        return Cerrors, Herrors
    
    elif t == 'r' or t == 'u':
        buf = inp[1].split(';')
        Cregions = [float(x) for x in buf[0].split(',')]
        Cerrors = [[float(x) for x in y.split(',')] for y in buf[1:]]
        
        buf = inp[2].split(';')
        Hregions = [float(x) for x in buf[0].split(',')]
        Herrors = [[float(x) for x in y.split(',')] for y in buf[1:]]
        
        return Cregions, Cerrors, Hregions, Herrors


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
        pdp4 = pdp4*stats.norm(expect, stdev).pdf(e)
    return pdp4


#load MultiGaus data from file and calculate probabilities
def CalculateRMultiGaus(errors, shifts, ParamFile, t, nucleus):

    Cregions, Cmeans, Cstdevs, Hregions, Hmeans, Hstdevs = \
        ReadRMultiGausFile(ParamFile, t)
    
    if nucleus == 'C':
        regions, means, stdevs = Cregions, Cmeans, Cstdevs
    elif nucleus == 'H':
        regions, means, stdevs = Hregions, Hmeans, Hstdevs

    prob = 1.0
    for e, s in zip(errors, shifts):
        region = bisect.bisect_left(regions, s)
        prob = prob*MultiGaus(means[region], stdevs[region], e)
    return prob


#load MultiGaus data from file and calculate probabilities
def CalculateMultiGaus(errors, ParamFile, t, nucleus):

    Cmeans, Cstdevs, Hmeans, Hstdevs = ReadParamFile(ParamFile, t)
    if nucleus == 'C':
        means, stdevs = Cmeans, Cstdevs
    elif nucleus == 'H':
        means, stdevs = Hmeans, Hstdevs

    prob = 1.0
    for e in errors:
        prob = prob*MultiGaus(means, stdevs, e)
    return prob


def MultiGaus(means, stdevs, x):
    res = 0
    for mean, stdev in zip(means, stdevs):
        res += stats.norm(mean, stdev).pdf(x)

    return res/len(means)