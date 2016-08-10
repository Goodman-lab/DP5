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

#Standard DP4 parameters
meanC = 0.0
meanH = 0.0
stdevC = 2.269372270818724
stdevH = 0.18731058105269952

#DFT J value parameters
FULLJ_MEAN = -0.13133138905769429
FULLJ_STDEV = 0.79315485469419067

#DFT FC value parameters
FC_MEAN = -0.15436128540661589
FC_STDEV = 0.92117647579294348

#Karplus J value parameters
KARPLUS_MEAN = 0.015779860851173257
KARPLUS_STDEV = 1.5949855401519151

output = []
J_THRESH = 0.2

#Linear correction parameters for DFT and Karplus J values
J_INTERCEPT = -0.365587317
J_SLOPE = 1.1598043367
FC_INTERCEPT = -0.438405757
FC_SLOPE = 1.1701625846
KARPLUS_INTERCEPT = -0.316665406
KARPLUS_SLOPE = 0.97767834



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


def DP4j(Clabels, Cvalues, Hlabels, Hvalues, Cexp, Hexp, cJvals, cJlabels,
         settings):

    Print(str(Cexp))
    Print(str(Hexp))
    
    C_cdp4 = []
    H_cdp4 = []
    J_dp4 = []
    Comb_cdp4 = []
    Super_dp4 = []

    for isomer in range(0, len(Cvalues)):

        sortedClabels, sortedCvalues, sortedCexp, _ =\
            AssignExpNMR(Clabels, Cvalues[isomer], Cexp)
        scaledC = ScaleNMR(sortedCvalues, sortedCexp)

        sHlabels, sortedHvalues, sortedHexp, sExpJvalues, sCalcJvalues = \
            AssignExpNMRj(Hlabels, Hvalues[isomer], Hexp, cJvals[isomer], cJlabels)
        scaledH = ScaleNMR(sortedHvalues, sortedHexp)
        
        assignedExpJs, assignedCJs = AssignJvals(sExpJvalues, sCalcJvalues)
                
        if settings.JStatsParamFile == '':
            if settings.jKarplus:
                slope, intercept = KARPLUS_SLOPE, KARPLUS_INTERCEPT
            elif settings.jJ:
                slope, intercept = J_SLOPE, J_INTERCEPT
            elif settings.jFC:
                slope, intercept = FC_SLOPE, FC_INTERCEPT
        else:
            slope, intercept = ReadJParamFile(settings.JStatsParamFile,
                                             settings.JStatsModel)[:2]
        
        ScaledCJs = ScaleJvals(assignedCJs, slope, intercept)

        ScaledErrorsC = [scaledC[i] - sortedCexp[i]
                         for i in range(0, len(scaledC))]
        ErrorsC = [sortedCvalues[i] - sortedCexp[i]
                         for i in range(0, len(sortedCvalues))]
        ScaledErrorsH = [scaledH[i] - sortedHexp[i]
                         for i in range(0, len(scaledH))]
        ErrorsH = [sortedHvalues[i] - sortedHexp[i]
                         for i in range(0, len(sortedHvalues))]
        
        ScaledErrorsJ = CalculateJerrors(assignedExpJs, ScaledCJs)
        
        Print("\nAssigned shifts for isomer " + str(isomer+1) + ": ")
        PrintNMR('C', sortedClabels, sortedCvalues, scaledC, sortedCexp)
        Print("Max C error: " + format(max(ScaledErrorsC, key=abs), "6.2f"))
        PrintNMRj('H', sHlabels, sortedHvalues, scaledH, sortedHexp,
                  assignedExpJs, assignedCJs)
        Print("Max H error: " + format(max(ScaledErrorsH, key=abs), "6.2f"))

        Cprob, Hprob = CalcProbabilities(ScaledErrorsC, ScaledErrorsH, ErrorsC,
                                         ErrorsH, sortedCexp, sortedHexp, settings)
        C_cdp4.append(Cprob)
        H_cdp4.append(Hprob)
        
        if settings.JStatsModel == 'k':
            J_dp4.append(CalculateJKDE(ScaledErrorsJ, settings.JStatsParamFile, 'k'))
        else:
            if settings.JStatsParamFile != '':
                Jmean, Jstdev = ReadJParamFile(settings.JStatsParamFile,
                                                 settings.JStatsModel)[2:]
            elif settings.jKarplus:
                Jmean, Jstdev = KARPLUS_MEAN, KARPLUS_STDEV
            elif settings.jFC:
                Jmean, Jstdev = FC_MEAN, FC_STDEV
            elif settings.jJ:
                Jmean, Jstdev = FULLJ_MEAN, FULLJ_STDEV
            
            J_dp4.append(CalculatePDP4(ScaledErrorsJ, Jmean, Jstdev))
        
        Comb_cdp4.append(C_cdp4[-1]*H_cdp4[-1])
        Super_dp4.append(C_cdp4[-1]*H_cdp4[-1]*J_dp4[-1])
        Print("\nDP4 based on C: " + format(C_cdp4[-1], "6.2e"))
        Print("DP4 based on H: " + format(H_cdp4[-1], "6.2e"))
        Print("DP4 based on J: " + format(J_dp4[-1], "6.2e"))

    relCDP4 = [(100*x)/sum(C_cdp4) for x in C_cdp4]
    relHDP4 = [(100*x)/sum(H_cdp4) for x in H_cdp4]
    relJDP4 = [(100*x)/sum(J_dp4) for x in J_dp4]
    relCombDP4 = [(100*x)/sum(Comb_cdp4) for x in Comb_cdp4]
    relSuperDP4 = [(100*x)/sum(Super_dp4) for x in Super_dp4]

    PrintRelDP4('all available data', relSuperDP4)
    PrintRelDP4('both carbon and proton data', relCombDP4)
    PrintRelDP4('carbon data only', relCDP4)
    PrintRelDP4('proton data only', relHDP4)
    PrintRelDP4('J value data only', relJDP4)
    PrintJrms(J_dp4)
    
    return output


def CP3(Clabels, Cvalues1, Cvalues2, Hlabels, Hvalues1, Hvalues2, Cexp1, Cexp2,
        Hexp1, Hexp2, settings):
    deltaCalcsC = [a - b for a,b in zip(Cvalues1, Cvalues2)]
    deltaCalcsH = [a - b for a,b in zip(Cvalues1, Cvalues2)]
    deltaExpsC = [a - b for a,b in zip(Cexp1, Cexp2)]
    deltaExpsH = [a - b for a,b in zip(Hexp1, Hexp2)]
    
    top = sum([f3(calc,exp) for calc,exp in zip(deltaCalcsC, deltaExpsC)])
    bottom = sum([x**2 for x in deltaExpsC])
    CP3c = top/bottom
    
    top = sum([f3(calc,exp) for calc,exp in zip(deltaCalcsH, deltaExpsH)])
    bottom = sum([x**2 for x in deltaExpsH])
    CP3h = top/bottom
    

def f3(deltaExp, deltaCalc):
    if (deltaCalc/deltaExp) > 1.0:
        return deltaExp**3/deltaCalc
    else:
        return deltaExp*deltaCalc

def DP4bias(Clabels, Cvalues, Hlabels, Hvalues, Cexp, Hexp, settings):

    Print(str(Cexp))
    Print(str(Hexp))

    C_cdp4 = []
    H_cdp4 = []
    Comb_cdp4 = []
    
    LabelsC = []
    LabelsH = []
    AllC = []
    AllH = []
    AllErrorsC = []
    AllErrorsH = []
    AllScaledC = []
    AllScaledH = []
    AllScaledErrorsC = []
    AllScaledErrorsH = []
    Cexps = []
    Hexps = []
    
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
        
        LabelsC.append(sortedClabels)
        Cexps.append(sortedCexp)
        LabelsH.append(sortedHlabels)
        Hexps.append(sortedHexp)
        AllErrorsC.append(ErrorsC)
        AllErrorsH.append(ErrorsH)
        AllC.append(sortedCvalues)
        AllH.append(sortedHvalues)
        AllScaledC.append(scaledC)
        AllScaledH.append(scaledH)
        AllScaledErrorsC.append(ScaledErrorsC)
        AllScaledErrorsH.append(ScaledErrorsH)
   
    Cexp1 = Cexps[0]
    Hexp1 = Hexps[0]
    
    LabelsC1, AllErrorsCr = ReorderAllErrors(AllErrorsC, LabelsC)
    LabelsH1, AllErrorsHr = ReorderAllErrors(AllErrorsH, LabelsH)
    
    _, AllScaledErrorsCr = ReorderAllErrors(AllScaledErrorsC, LabelsC)
    _, AllScaledErrorsHr = ReorderAllErrors(AllScaledErrorsH, LabelsH)
    
    Print('\n'.join([','.join(x) for x in LabelsC]))
    Print("Scaled C errors before reordering:")
    PrintErrors(LabelsC1, AllScaledErrorsCr)
        
    Print("Scaled errors after reordering:")
    PrintErrors(LabelsC1, AllScaledErrorsCr)
    PrintErrors(LabelsH1, AllScaledErrorsHr)
    
    BiasesC = CalcBiases(AllScaledErrorsCr)
    BiasesH = CalcBiases(AllScaledErrorsHr)

    uBiasesC = CalcBiases(AllErrorsCr)
    uBiasesH = CalcBiases(AllErrorsHr)
    uSpreadsC = CalcSpreads(AllErrorsCr)
    uSpreadsH = CalcSpreads(AllErrorsHr)
    
    Print("Biases for unscaled C:")
    Print(','.join([format(x, "6.2f") for x in uBiasesC]))
    Print("Biases for unscaled H:")
    Print(','.join([format(x, "6.2f") for x in uBiasesH]))
    
    Print("Spreads for unscaled C:")
    Print(','.join([format(x, "6.2f") for x in uSpreadsC]))
    Print("Spreads for unscaled H:")
    Print(','.join([format(x, "6.2f") for x in uSpreadsH]))
    
    Print("Corresponding C exp shifts:")
    Print(','.join([format(x, "6.2f") for x in Cexp1]))
    Print("Corresponding H exp shifts:")
    Print(','.join([format(x, "6.2f") for x in Hexp1]))
    
    AllScaledErrorsCb = ApplyBias(AllScaledErrorsCr, BiasesC)
    AllScaledErrorsHb = ApplyBias(AllScaledErrorsHr, BiasesH)
    AllErrorsCb = ApplyBias(AllErrorsCr, uBiasesC)
    AllErrorsHb = ApplyBias(AllErrorsHr, uBiasesH)
    
    Print("Scaled errors after biasing:")
    PrintErrors(LabelsC1, AllScaledErrorsCb)
    PrintErrors(LabelsH1, AllScaledErrorsHb)
    
    for i in range(len(AllC)):
        #print biased shift data for every diastereomer
        bCvalues = ApplyBiasShifts(LabelsC[i], AllC[i], LabelsC1, uBiasesC)
        sbCvalues = ApplyBiasShifts(LabelsC[i], AllScaledC[i], LabelsC1, BiasesC)
        bHvalues = ApplyBiasShifts(LabelsH[i], AllH[i], LabelsH1, uBiasesH)
        sbHvalues = ApplyBiasShifts(LabelsH[i], AllScaledH[i], LabelsH1, BiasesH)
        Print("\nAssigned shifts for isomer " + str(i+1) + ": ")
        PrintNMR('C', LabelsC[i], bCvalues, sbCvalues, Cexps[i])
        Print("Max C error: " + format(max(AllScaledErrorsCb[i], key=abs), "6.2f"))
        PrintNMR('H', LabelsH[i], bHvalues, sbHvalues, Hexps[i])
        Print("Max H error: " + format(max(AllScaledErrorsHb[i], key=abs), "6.2f"))

    for i in range(len(AllScaledErrorsCb)):
        Cprob, Hprob = CalcProbabilities(AllScaledErrorsCb[i], AllScaledErrorsHb[i],
            AllErrorsCb[i], AllErrorsHb[i], Cexps[i], Hexps[i], settings)
        C_cdp4.append(Cprob)
        H_cdp4.append(Hprob)
        Comb_cdp4.append(C_cdp4[-1]*H_cdp4[-1])
    
    relCDP4 = [(100*x)/sum(C_cdp4) for x in C_cdp4]
    relHDP4 = [(100*x)/sum(H_cdp4) for x in H_cdp4]
    relCombDP4 = [(100*x)/sum(Comb_cdp4) for x in Comb_cdp4]

    PrintRelDP4('all available data', relCombDP4)
    PrintRelDP4('carbon data only', relCDP4)
    PrintRelDP4('proton data only', relHDP4)

    return output


def ApplyBiasShifts(Labels, Shifts, BiasLabels, biases):
    BiasedShifts = []
    for i in range(len(Shifts)):
        j = BiasLabels.index(Labels[i])
        BiasedShifts.append(Shifts[i]-biases[j])
    
    return BiasedShifts


def ApplyBias(AllErrors, biases):
    
    AllErrorsB = []
    for errors in AllErrors:
        AllErrorsB.append([x-y for x,y in zip(errors,biases)])
    
    return AllErrorsB


def CalcBiases(AllErrors):
    
    biases = []
    for i in range(len(AllErrors[0])):
        serrors = sorted([AllErrors[x][i] for x in range(len(AllErrors))], key=abs)
        if (serrors[0]<0) == (serrors[1]<0):
            biases.append(serrors[0])
        else:
            biases.append(0.0)
        #biases.append(min([AllErrors[x][i] for x in range(len(AllErrors))], key=abs))
    
    return biases

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
            AllErrorsR.append(ReorderErrors(Labels1, labels, errors))
    
    return Labels1, AllErrorsR

#given 2 sets of labels and one set of values, the set of values will be
#reordered from the order in labels2 to the order in labels1
def ReorderErrors(labels1, labels2, errors2):
    
    newErrors = []
    templabels2 = list(labels2)
    
    for l in labels1:
        index = templabels2.index(l)
        newErrors.append(errors2[index])
        #Remove the value to avoid adding the same value twice
        templabels2[index] = 'XX'
    
    return newErrors


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


def AssignExpNMRj(labels, calcShifts, exp, cJvalues, cJlabels):
    
    expLabels, expShifts, expJs = ReadExp(exp)
    #Prepare sorted calculated data with labels and sorted exp data
    sortedCalc = sorted(calcShifts)
    sortedExp = sorted(expShifts)
    sortedExpLabels = SortExpAssignments(expShifts, expLabels)
    sortedCalcLabels = []
    sortedExpJs = []
    
    for v in sortedExp:
        index = expShifts.index(v)
        sortedExpJs.append(expJs[index])
 
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
    sortedCJs = [[0.0] for x in assignedExpLabels]
    
    #Rearrange calc values to match the assigned labels
    for i, l in enumerate(assignedExpLabels):
        index = labels.index(l)
        sortedCalc.append(calcShifts[index])
        if l[1:] in cJlabels:
            sortedCJs[i] = cJvalues[cJlabels.index(l[1:])]

    return assignedExpLabels, sortedCalc, sortedExp, sortedExpJs, sortedCJs


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
    
    print expJs
    #Remove assignments and get shifts
    ShiftData = (re.sub(r"\(.*?\)", "", exp, flags=re.DOTALL)).split(',')
    expShifts = [float(x) for x in ShiftData]
    
    return expLabels, expShifts, expJs


def AssignJvals(expJ, calcJ):
    
    import itertools
    prunedCJs = []
    for i in calcJ:
        prunedCJs.append([abs(j) for j in i if abs(j)>J_THRESH])
    
    #Check that number of calculated J values are not less than experimental
    #relevant if using Karplus equation
    prunedExpJ = list(expJ)
    
    for i in range(len(expJ)):
        if len(prunedCJs[i]) == 0:
            prunedExpJ[i] = [0.0]
        #if len(prunedCJs[i]) < len(prunedExpJ[i]):
        #    prunedCJs[i].extend([0.0 for i in range(len(prunedExpJ[i])-len(prunedCJs[i]))])
    
    assignedCalcJ = []
    
    #IMPORTANT:
    #Deal with cases the number of calc J values are less than exp
    for eJ, cJ in zip(prunedExpJ, prunedCJs):
        if eJ != [0.0]:
            
            minerror = 10000
            minassign = [0.0 for i in range(len(eJ))]
            for assign in itertools.permutations(cJ, len(eJ)):
                error = sum([abs(assign[i]-eJ[i]) for i in range(len(eJ))])
                if error<minerror:
                    minerror = error
                    minassign = assign
            assignedCalcJ.append(list(minassign))
        else:
            assignedCalcJ.append([0.0])
        
    return prunedExpJ, assignedCalcJ


def ScaleJvals(calcJ, slope, intercept):
    scaledJs = []
    
    for cJ in calcJ:
        scaledJs.append([(j-intercept)/slope for j in cJ])
    
    return scaledJs


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
            
    elif settings.StatsModel == 'k':
        C_cdp4 = CalculateKDE(SErrorsC, settings.StatsParamFile,
                                   settings.StatsModel, 'C')
        H_cdp4 = CalculateKDE(SErrorsH, settings.StatsParamFile,
                                   settings.StatsModel, 'H')

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
    elif settings.StatsModel == 'r':
        C_cdp4 = CalculateRKDE(SErrorsC, Cexp, settings.StatsParamFile,
                                    settings.StatsModel, 'C')
        H_cdp4 = CalculateRKDE(SErrorsH, Hexp, settings.StatsParamFile,
                                    settings.StatsModel, 'H')
        
    elif settings.StatsModel == 'u':
        C_cdp4 = CalculateRKDE(ErrorsC, Cexp, settings.StatsParamFile,
                                    settings.StatsModel, 'C')
        H_cdp4 = CalculateRKDE(ErrorsH, Hexp, settings.StatsParamFile,
                                    settings.StatsModel, 'H')
    else:
        print "Invalid stats model"
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


def PrintNMRj(nucleus, labels, values, scaled, exp, expJ, calcJ):
    Print("\nAssigned " + nucleus +
          " shifts: (label, calc, corrected, exp, error, Jvalues)")
    for i in range(len(labels)):
        if expJ[i] != [0.0]:
            expJstring = ", ".join(["{:4.1f}".format(x) for x in expJ[i]])
            calcJstring = ", ".join(["{:4.1f}".format(x) for x in calcJ[i] if abs(x)>J_THRESH])
            Jinfo = expJstring + ' || ' + calcJstring
        else:
            Jinfo = ''
        Print(format(labels[i], "6s") + ' ' + format(values[i], "6.2f") + ' '
            + format(scaled[i], "6.2f") + ' ' + format(exp[i], "6.2f") + ' ' +
            format(exp[i]-scaled[i], "6.2f") + '   ' + Jinfo)
    Print("Direct J comparison table (exp, calc):")
    Jerrs = []
    for i in range(len(labels)):
        if expJ[i] != [0.0]:
            for e, c in zip(expJ[i], calcJ[i]):
                print "{:4.1f}".format(e) + "   " + "{:4.1f}".format(c)
                Jerrs.append(e-c)
    Print("Max J error:" + "{:4.1f}".format(max(Jerrs, key=abs)))


def PrintRelDP4(title, RelDP4):
    Print("\nResults of DP4 using " + title + ":")
    for i in range(len(RelDP4)):
        Print("Isomer " + str(i+1) + ": " + format(RelDP4[i], "4.1f") + "%")


def PrintJ(Jerrors):
    Print("\nMean absolute errors of J values:")
    for i, Jerr in enumerate(Jerrors):
        Print("Isomer " + str(i+1) + ": " + format(sum(Jerr)/len(Jerr), "5.3f"))


def PrintJrms(Jerrors):
    from numpy import sqrt, mean, square
    Print("\nMean absolute errors of J values:")
    for i, Jerr in enumerate(Jerrors):
        Print("Isomer " + str(i+1) + ": " + format(sqrt(mean(square(Jerr))), "5.3f"))


def Print(s):
    print s
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
        print "Wrong parameter file type, exiting..."
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


def ReadJParamFile(f, t):
    
    infile = open(f, 'r')
    inp = infile.readlines()
    infile.close()
    
    if t not in inp[0]:
        print "Wrong parameter file type, exiting..."
        quit()
    
    if t == 'g': #Gaussian model for J values, includes linear scaling factors
        [slope, intercept] = [float(x) for x in inp[1].split(',')]
        [Jmean, Jstdev] = [float(x) for x in inp[2].split(',')]
        return slope, intercept, Jmean, Jstdev
    
    elif t == 'k': #KDE model for J values, includes linear scaling factors
        [slope, intercept] = [float(x) for x in inp[1].split(',')]
        Jerrors = [float(x) for x in inp[2].split(',')]
        return slope, intercept, Jerrors
    
        
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


#load empirical error data from file and use KDE to construct pdf
def CalculateKDE(errors, ParamFile, t, nucleus):

    Cerrors, Herrors = ReadParamFile(ParamFile, 'k')
        
    if nucleus == 'C':
        kde = stats.gaussian_kde(Cerrors)
    elif nucleus == 'H':
        kde = stats.gaussian_kde(Herrors)

    prob = 1.0
    for e in errors:
        prob = prob*float(kde(e)[0])
    return prob


#load empirical error data from file and use KDE to construct pdf
def CalculateJKDE(errors, ParamFile, t):

    slope, intercept, Jerrors = ReadJParamFile(ParamFile, t)
    kde = stats.gaussian_kde(Jerrors)
    
    ep5 = 1.0
    for e in errors:
        ep5 = ep5*float(kde(e)[0])
    return ep5

#load empirical error data from file and use KDE to construct several pdfs,
#one for each chemical shift region, can be used both for scaled and unscaled
#errors with the appropriate data
def CalculateRKDE(errors, shifts, ParamFile, t, nucleus):
    
    #Load the data
    Cregions, Cerrors, Hregions, Herrors = ReadParamFile(ParamFile, t)
    if nucleus == 'C':
        regions = Cregions
        RErrors = Cerrors
    elif nucleus == 'H':
        regions = Hregions
        RErrors = Herrors
        
    #Reconstruct the distributions for each region
    kdes = []
    for es in RErrors:
        kdes.append(stats.gaussian_kde(es))

    ep5 = 1.0
    for i, e in enumerate(errors):
        region = bisect.bisect_left(regions, shifts[i])
        ep5 = ep5*float((kdes[region])(e)[0])
    return ep5


def CalculateJerrors(expJ, calcJ):
    
    errors = []
    for eJ, cJ in zip(expJ, calcJ):
        if eJ != [0.0]:
            errors.extend([abs(eJ[i] - cJ[i]) for i in range(len(eJ))])
    
    return errors
