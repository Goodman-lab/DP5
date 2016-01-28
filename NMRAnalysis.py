# -*- coding: utf-8 -*-
"""
Created on Mon Jan  12 14:42:47 2015

@author: ke291

Takes care of all the NMR description interpretation, equivalent atom
averaging, Boltzmann averaging, tautomer population optimisation (if used)
and DP4 input preparation and running DP4.py. Called by PyDP4.py
"""

import Gaussian
import NWChem

import sys
import scipy.optimize as sciopt
from numpy import mean
import math
import os

J_THRESHOLD = 0.1


def main(numDS, settings, *args):

    #This function runs nmrPredict for each diastereomer and collects
    #the outputs
    print '\nRunning NMRpredict script...'

    if settings.DFT == 'z' or settings.DFT == 'g' or settings.DFT == 'd':
        (RelEs, populations, labels, BoltzmannShieldings, Jlabels, BoltzmannFCs,
            BoltzmannJs, SigConfs, Ntaut) =  \
            Gaussian.RunNMRPredict(numDS, settings, *args)
    elif settings.DFT == 'n' or settings.DFT == 'w' or settings.DFT == 'm':
        (RelEs, populations, labels, BoltzmannShieldings, SigConfs, Ntaut) = \
                                            NWChem.RunNMRPredict(numDS, *args)
    
    if settings.RenumberFile != '':
        BoltzmannShieldings = ReorderShieldings(BoltzmannShieldings,
                                                settings.RenumberFile)
    
    #Reads the experimental NMR data from the file
    Cexp, Hexp, equivs, omits = ReadExpNMR(args[-1])
    
    for Es, pops in zip(RelEs, populations):
        print '\nConformer relative energies (kJ/mol): ' + \
            ', '.join(["{:5.2f}".format(float(x)) for x in Es])

        print '\nPopulations (%): ' + \
            ', '.join(["{:4.1f}".format(float(x)*100) for x in pops])
    
    #Convert shielding constants in chemical shifts and sort labels by nuclei
    Cvalues, Hvalues, Clabels, Hlabels = GetCalcShiftsLabels(numDS,
        BoltzmannShieldings, labels, omits, settings)
    
    if settings.OtherNuclei != "":
        Xvalues, Xlabels = GetOtherNuclei(numDS, BoltzmannShieldings, labels,
            omits, settings)
        print Xlabels, Xvalues
    else:
        Xvalues = []
        Xlabels = []
    
    Noutp = len(BoltzmannShieldings)
    
    #Looks for equivalent atoms in the computational data, averages the shifts
    #and removes the redundant signals
    Cvalues, Hvalues, Clabels, Hlabels = \
        RemoveEquivalents(Noutp, equivs, Cvalues, Hvalues, Clabels, Hlabels)

    tstart = 0
    OptCvalues = []
    OptHvalues = []

    for tindex in range(0, len(Ntaut)):
        print 'looking at tautomers ' + str(tstart) + ' to ' + \
            str(tstart+Ntaut[tindex])
        if Ntaut[tindex] == 1:
            print "Only one tautomer found, skipping optimisation."
            OptCvalues.append(Cvalues[tstart])
            OptHvalues.append(Hvalues[tstart])
            tstart = tstart + Ntaut[tindex]
        else:
            (BuffC, BuffH) = OptTautPop(Clabels,
                                        Cvalues[tstart:tstart+Ntaut[tindex]],
                                        Hlabels,
                                        Hvalues[tstart:tstart+Ntaut[tindex]],
                                        Cexp, Hexp)
            OptCvalues.append(BuffC)
            OptHvalues.append(BuffH)
            tstart = tstart + Ntaut[tindex]
    
    if any([settings.jKarplus, settings.jJ, settings.jFC]):
        NewBJs, NewJlabels = ZeroEquivJ(BoltzmannJs, Jlabels, equivs, omits)
        NewFCs, NewFClabels = ZeroEquivJ(BoltzmannFCs, Jlabels, equivs, omits)
    
    print "The calculated data for other nuclei:"
    PrintOtherNuclei(numDS, Xlabels, Xvalues)
    
    print "Conformation data:"
    PrintConformationData(SigConfs)
    
    import DP4
    #Run DP4 (or alternative, if set in settings) analysis and collect output
    if any([settings.jKarplus, settings.jJ, settings.jFC]):
        print "\n J value matrixes after pruning: \n"
        if settings.jFC:
            for i, Jvals in enumerate(NewFCs):
                print "Isomer " + str(i) + ":"
                PrintJMatrixLim(Jvals, NewFClabels)
        else:
            for i, Jvals in enumerate(NewBJs):
                print "Isomer " + str(i) + ":"
                PrintJMatrixLim(Jvals, NewJlabels)
        if settings.jFC:
            DP4outp = DP4.DP4j(Clabels, OptCvalues, Hlabels, OptHvalues, Cexp,
                               Hexp, NewFCs, NewFClabels, settings)
        else:
            DP4outp = DP4.DP4j(Clabels, OptCvalues, Hlabels, OptHvalues, Cexp,
                               Hexp, NewBJs, NewJlabels, settings)
    
        return '\n'.join(DP4outp) + '\n'
    else:
        DP4outp = DP4.DP4(Clabels, OptCvalues, Hlabels, OptHvalues, Cexp,
                               Hexp, settings)
        return '\n'.join(DP4outp) + '\n'


def ReadExpNMR(ExpNMR):
    #Reads the experimental NMR data from the file
    ExpNMR_file = open(ExpNMR, 'r')
    Cexp = ExpNMR_file.readline()
    ExpNMR_file.readline()
    Hexp = ExpNMR_file.readline()

    #Check if exp NMR file contains info about equivalent atoms and read it
    #into an array
    #Also reads a list of atoms to omit from analysis

    equivalents = []
    omits = []

    ExpNMR_file.readline()
    for line in ExpNMR_file:
        if not 'OMIT' in line and len(line) > 1:
            equivalents.append(line[:-1].split(','))
        elif 'OMIT' in line:
            omits.extend(line[5:-1].split(','))
            
    ExpNMR_file.close()
    
    return Cexp, Hexp, equivalents, omits


def GetCalcShiftsLabels(numDS, BShieldings, labels, omits, settings):
    
    Clabels = []
    Hlabels = []
    Cvalues = []
    Hvalues = []
    
    for DS in range(numDS):

        Cvalues.append([])
        Hvalues.append([])

        #loops through particular output and collects shielding constants
        #and calculates shifts relative to TMS
        for atom in range(len(BShieldings[DS])):
            shift = 0
            if labels[atom][0] == 'C' and not labels[atom] in omits:
                # only read labels once, i.e. the first diastereomer
                if DS == 0:
                    Clabels.append(labels[atom])
                shift = (settings.TMS_SC_C13-BShieldings[DS][atom]) / \
                    (1-(settings.TMS_SC_C13/10**6))
                Cvalues[DS].append(shift)

            if labels[atom][0] == 'H' and not labels[atom] in omits:
                # only read labels once, i.e. the first diastereomer
                if DS == 0:
                    Hlabels.append(labels[atom])
                shift = (settings.TMS_SC_H1-BShieldings[DS][atom]) / \
                    (1-(settings.TMS_SC_H1/10**6))
                Hvalues[DS].append(shift)

    return Cvalues, Hvalues, Clabels, Hlabels


def GetOtherNuclei(numDS, BShieldings, labels, omits, settings):
    Xlabels = []
    Xvalues = []
    
    for DS in range(numDS):

        Xvalues.append([])

        #loops through particular output and collects shielding constants
        #and calculates shifts relative to TMS
        for atom in range(len(BShieldings[DS])):
            shift = 0
            if labels[atom][0] in settings.OtherNuclei:
                
                # only read labels once, i.e. the first diastereomer
                if DS == 0:
                    Xlabels.append(labels[atom])
                shift = (settings.CFCl3_SC_F19-BShieldings[DS][atom]) / \
                    (1-(settings.CFCl3_SC_F19/10**6))
                Xvalues[DS].append(shift)

    return Xvalues, Xlabels


def PrintOtherNuclei(numDS, Xlabels, Xvalues):
    
    for DS in range(numDS):
        
        print "\nOther nuclei results for isomer " + str(DS+1) + ":"
        for i, label in enumerate(Xlabels):
            print label + " " + format(Xvalues[DS][i], "4.2f")


def PrintConformationData(AllSigConfs):
    for i, SigConfs in enumerate(AllSigConfs):
        print "\nNumber of significant conformers for isomer "\
            + str(i+1) + ": " + str(len(SigConfs)) + "\n(pop, filename)"
        for conf in SigConfs:
            print "   " + format(conf[1]*100, "4.2f") + "%   " + conf[0]
        print '----------------'
        print "   " + format(100*sum([x[1] for x in SigConfs]), "4.2f") +\
            "%   in total"


def RemoveEquivalents(Noutp, equivs, OldCval, OldHval, OldClabels, OldHlabels):
    Cvalues = list(OldCval)
    Hvalues = list(OldHval)
    Clabels = list(OldClabels)
    Hlabels = list(OldHlabels)
    
    for eqAtoms in equivs:

        eqSums = [0.0]*Noutp
        eqAvgs = [0.0]*Noutp

        if eqAtoms[0][0] == 'H':
            #print eqAtoms, Hlabels
            for atom in eqAtoms:
                eqIndex = Hlabels.index(atom)
                for ds in range(0, Noutp):
                    eqSums[ds] = eqSums[ds] + Hvalues[ds][eqIndex]
            for ds in range(0, Noutp):
                eqAvgs[ds] = eqSums[ds]/len(eqAtoms)

            #Place the new average value in the first atom shifts place
            target_index = Hlabels.index(eqAtoms[0])
            for ds in range(0, Noutp):
                Hvalues[ds][target_index] = eqAvgs[ds]

            #Delete the redundant atoms from the computed list
            #start with second atom - e.g. don't delete the original one
            for atom in range(1, len(eqAtoms)):
                del_index = Hlabels.index(eqAtoms[atom])
                del Hlabels[del_index]
                for ds in range(0, Noutp):
                    del Hvalues[ds][del_index]

        if eqAtoms[0][0] == 'C':
            for atom in eqAtoms:
                eqIndex = Clabels.index(atom)
                for ds in range(0, Noutp):
                    eqSums[ds] = eqSums[ds] + Cvalues[ds][eqIndex]
            for ds in range(0, Noutp):
                eqAvgs[ds] = eqSums[ds]/len(eqAtoms)

            #Place the new average value in the first atom shifts place
            target_index = Clabels.index(eqAtoms[0])
            for ds in range(0, Noutp):
                Cvalues[ds][target_index] = eqAvgs[ds]

            #Delete the redundant atoms from the computed list
            #start with second atom - e.g. don't delete the original one
            for atom in range(1, len(eqAtoms)):
                del_index = Clabels.index(eqAtoms[atom])
                del Clabels[del_index]
                for ds in range(0, Noutp):
                    del Cvalues[ds][del_index]
                    
    return Cvalues, Hvalues, Clabels, Hlabels
    

def ReorderShieldings(shieldings, RenumberFile):
    f = open(RenumberFile, 'r')
    RenumData = f.readlines()
    
    RenumMaps = []
    for line in RenumData:
        tmp = line.split(',')
        RenumMaps.append([int(x)-1 for x in tmp])
    
    print RenumMaps
    
    ReorderedShieldings = []
    
    for i, shields in enumerate(shieldings):
        if i != 0:
            NewShields = [0.0 for x in range(len(shields))]
            for j in range(len(NewShields)):
                NewShields[j] = shields[RenumMaps[i-1][j]]
            ReorderedShieldings.append(NewShields)
        else:
            ReorderedShieldings.append(shields)
    
    print "Before reordering:"
    for s in shieldings:
        print ','.join([format(x, "4.2f") for x in s]) + '\n'
    
    print "After reordering:"
    for s in ReorderedShieldings:
        print ','.join([format(x, "4.2f") for x in s]) + '\n'
        
    return ReorderedShieldings


def ZeroSmallJvalues(mat):
    
    newmat = list(mat)
    for ds in range(len(newmat)):
        for row in range(len(newmat[ds])):
            for column in range(len(newmat[ds])):
                if abs(newmat[ds][row][column])<J_THRESHOLD:
                    newmat[ds][row][column] = 0.0
    return newmat


def ZeroEquivJ(mat, matlabels, equivs, omits):
    
    #Zeros mutual coupling constants of equivalent atoms
    newmat = list(mat)
    toRemove = []
    print equivs
    print matlabels
    for eqAtoms in equivs:
        AllAtomsPresent = False
        AllAtomsPresent = all([x[1:] in matlabels for x in eqAtoms])
        if eqAtoms[0][0] == 'H' and AllAtomsPresent:
            indexes = [matlabels.index(x[1:]) for x in eqAtoms]
            
            #Zero out mutual
            for row in indexes:
                for column in indexes:
                    for ds in range(len(newmat)):
                        newmat[ds][row][column] = 0.0
            
            #Average between equivalent
            first = indexes[0]
            for ds in range(len(newmat)):
                #Average the equivs in rows
                for row in range(len(newmat[ds])):
                    newmat[ds][row][first] = \
                        mean([newmat[ds][row][x] for x in indexes])
                #Average the equivs in columns
                for column in range(len(newmat[ds])):
                    newmat[ds][first][column] = \
                        mean([newmat[ds][x][column] for x in indexes])
            
            #Save the indexes to be removed from the matrix
            #for e in indexes[1:]: toRemove.append(e)
            toRemove.extend(indexes[1:])
    
    toRemove.extend([matlabels.index(x[1:]) for x in omits if x[0]=='H' and
                     x[1:] in matlabels])
    
    PrunedMats = []
    for ds in range(len(newmat)):
        PrunedMat, PrunedLabels = RemoveAtomsMatrix(newmat[ds], matlabels,
                                                    toRemove)
        PrunedMats.append(PrunedMat)
    
    return PrunedMats, PrunedLabels


def RemoveAtomsMatrix(mat, matlabels, ToRemove):
    
    if len(mat) < 2:
        return mat, matlabels

    PrunedMatrix = []
    PrunedLabels = []
    for y in range(len(mat)):
        row = []
        if y not in ToRemove:
            PrunedLabels.append(matlabels[y])
            for x in range(len(mat)):
                if x not in ToRemove:
                    row.append(mat[y][x])
            PrunedMatrix.append(row)
    return PrunedMatrix, PrunedLabels

def PrintJMatrix(mat, labels):
    
    for l, row in zip(labels, mat):
        formattedJ = ["{:4.1f}".format(x) for x in row]
        print l + ': ' + ', '.join(formattedJ)


def PrintJMatrixLim(mat, labels):
    
    for l, row in zip(labels, mat):
        formattedJ = ["{:4.1f}".format(x) for x in row if abs(x)>J_THRESHOLD]
        print l + ': ' + ', '.join(formattedJ)


def OptTautPop(Clabels, Cvalues, Hlabels, Hvalues, Cexp, Hexp):
    #Pairwise match exp signals to computed ones based on assignments first,
    #on erorrs afterwards
    ExpCvalues = [-1 for i in range(0, len(Clabels))]
    ExpHvalues = [-1 for i in range(0, len(Hlabels))]

    Hdata = Hexp.split(',')
    for s in range(0, len(Hdata)):
        Hdata[s] = Hdata[s].strip()

    Cdata = Cexp.split(',')
    for s in range(0, len(Cdata)):
        Cdata[s] = Cdata[s].strip()

    UAExpCshifts = list(Cdata)
    UAExpHshifts = list(Hdata)

    #Assign known(experimentally assigned) signals first
    for l in range(0, len(Clabels)):
        for s in Cdata:
            if (Clabels[l] + ')') in s and (not 'or' in s) and (not 'OR' in s):
                shiftend = s.find('(')
                ExpCvalues[l] = float(s[:shiftend])
                UAExpCshifts.remove(s)
                break
    for l in range(0, len(Hlabels)):
        for s in Hdata:
            if (Hlabels[l] + ')') in s and (not 'or' in s) and (not 'OR' in s):
                shiftend = s.find('(')
                ExpHvalues[l] = float(s[:shiftend])
                UAExpHshifts.remove(s)
                break

    #Prepare unassigned experimental values for matching
    for i in range(0, len(UAExpHshifts)):
        shiftend = UAExpHshifts[i].find('(')
        UAExpHshifts[i] = float(UAExpHshifts[i][:shiftend])

    for i in range(0, len(UAExpCshifts)):
        shiftend = UAExpCshifts[i].find('(')
        UAExpCshifts[i] = float(UAExpCshifts[i][:shiftend])

    #Try to assign unassigned values based on every calculated tautomer
    MinMAE = 1000
    for k in range(0, len(Cvalues)):
        #Pick out unassigned computational values for matching
        UACompCshifts = []
        UACompHshifts = []
        for i in range(0, len(ExpCvalues)):
            if ExpCvalues[i] == -1:
                UACompCshifts.append(Cvalues[k][i])
        for i in range(0, len(ExpHvalues)):
            if ExpHvalues[i] == -1:
                UACompHshifts.append(Hvalues[k][i])

        #Sort both sets of data - this essentially pairs them up
        UAExpCshifts.sort()
        UAExpHshifts.sort()
        UACompCshifts.sort()
        UACompHshifts.sort()

        #Go through half-assigned experimental data and fill in the holes with
        #paired data
        for i in range(0, len(ExpCvalues)):
            if ExpCvalues[i] == -1 and len(UAExpCshifts) > 0:
                j = UACompCshifts.index(Cvalues[k][i])
                ExpCvalues[i] = UAExpCshifts[j]
        for i in range(0, len(ExpHvalues)):
            if ExpHvalues[i] == -1 and len(UAExpHshifts) > 0:
                j = UACompHshifts.index(Hvalues[k][i])
                ExpHvalues[i] = UAExpHshifts[j]

        #Optimize tautomer populations for
        #tpops is 1 shorter than the number of tautomers,
        #the remaining weight is 1-sum(rest)
        tpops = [1.0/len(Cvalues) for i in range(len(Cvalues)-1)]
        f = lambda w: TautError(Cvalues, ExpCvalues, Hvalues, ExpHvalues, w)
        res = sciopt.minimize(f, tpops, method='nelder-mead')
        if float(res.fun) < MinMAE:
            print "New min MAE: " + str(res.fun)
            MinMAE = float(res.fun)
            NewPops = list(res.x)
            NewPops.append(1-sum(NewPops))
            print NewPops

    NewCvalues = []
    NewHvalues = []

    #calculate the new C values
    for atom in range(0, len(Clabels)):
        C = 0
        for taut in range(0, len(Cvalues)):
            C = C + Cvalues[taut][atom]*NewPops[taut]
        NewCvalues.append(C)

    for atom in range(0, len(Hlabels)):
        H = 0
        for taut in range(0, len(Hvalues)):
            H = H + Hvalues[taut][atom]*NewPops[taut]
        NewHvalues.append(H)

    #Return the new Cvalues and Hvalues
    return (NewCvalues, NewHvalues)


def TautError(Cs, CExp, Hs, HExp, TPopsIn):

    if len(Cs) != len(TPopsIn)+1 or len(Hs) != len(TPopsIn)+1:
        print len(Cs), len(Hs), len(TPopsIn)
        print ("Input dimensions in TautError don't match, exiting...")
        return 1000
    TPops = list(TPopsIn)
    TPops.append(1-sum(TPops))
    SumC = []
    SumH = []
    #print len(Cs), len(Hs), len(TPops)
    #print TPops
    for i in range(0, len(TPops)):
        if TPops[i] < 0:
            return 100
    if sum(TPops) > 1:
        s = sum(TPops)
        for i in range(0, len(TPops)):
            TPops[i] = TPops[i]/s

    for atom in range(0, len(CExp)):
        C = 0
        for taut in range(0, len(Cs)):
            C = C + Cs[taut][atom]*TPops[taut]
        SumC.append(C)

    for atom in range(0, len(HExp)):
        H = 0
        for taut in range(0, len(Hs)):
            H = H + Hs[taut][atom]*TPops[taut]
        SumH.append(H)

    ErrC = MAE(SumC, CExp)
    #ErrC = RMSE(SumC, CExp)
    ErrH = MAE(SumH, HExp)
    #ErrH = RMSE(SumH, HExp)
    #print 'MAE for C: ' + str(ErrC)
    #print 'MAE for H: ' + str(ErrH)

    return ErrC + 20*ErrH


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


def MAE(L1, L2):

    if len(L1) != len(L2):
        return -1
    else:
        L = []
        for i in range(0, len(L1)):
            L.append(abs(L1[i]-L2[i]))
        return sum(L)/len(L)


def RMSE(L1, L2):

    if len(L1) != len(L2):
        return -1
    else:
        L = []
        for i in range(0, len(L1)):
            L.append((L1[i]-L2[i])**2)
        return math.sqrt(sum(L)/len(L))

if __name__ == '__main__':
    #print sys.argv
    cpargs = sys.argv[2:]
    numDS = int(sys.argv[1])
    for ds in range(0, numDS):
        cpargs[ds*2+1] = int(cpargs[ds*2+1])
    main(numDS, *cpargs)
