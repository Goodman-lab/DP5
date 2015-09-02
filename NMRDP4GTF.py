# -*- coding: utf-8 -*-
"""
Created on Mon Jan  12 14:42:47 2015

@author: ke291

Takes care of all the NMR description interpretation, equivalent atom
averaging, Boltzmann averaging, tautomer population optimisation (if used)
and DP4 input preparation and running either DP4.jar or DP4.py. Called by
PyDP4.py
"""

import Gaussian
import NWChem

import subprocess
import sys
import scipy.optimize as sciopt
import math
import os

TMS_SC_C13 = 191.69255
TMS_SC_H1 = 31.7518583


def main(numDS, settings, *args):

    #This function runs nmrPredict for each diastereomer and collects
    #the outputs
    print '\nRunning NMRpredict script...'

    if settings.DFT == 'z' or settings.DFT == 'g':
        (RelEs, populations, labels, BoltzmannShieldings, Ntaut) = \
                                            Gaussian.RunNMRPredict(numDS, *args)
        Noutp = len(BoltzmannShieldings)
    elif settings.DFT == 'n' or settings.DFT == 'w':
        (RelEs, populations, labels, BoltzmannShieldings, Ntaut) = \
                                            NWChem.RunNMRPredict(numDS, *args)
        Noutp = len(BoltzmannShieldings)

    #Reads the experimental NMR data from the file
    ExpNMR = args[-1]
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
        else:
            omits.append(line[5:-1].split(','))
    ExpNMR_file.close()

    Clabels = []
    Hlabels = []
    Cvalues = []
    Hvalues = []

    #This loops through predictNMR outputs for each diastereomer and collects
    #NMR data

    flatequiv = [val for sublist in equivalents for val in sublist]    
    flatomits = [val for sublist in omits for val in sublist]

    for DS in range(0, numDS):

        Cvalues.append([])
        Hvalues.append([])

        buf = ''
        for val in RelEs[DS]:
            num = float(val)
            buf = buf + "{:5.2f}".format(num) + ', '
        print '\nConformer relative energies (kJ/mol): ' + buf[:-2]

        buf = ''
        for val in populations[DS]:
            num = float(val)
            buf = buf + "{:4.1f}".format(num*100) + ', '
        print '\nPopulations (%): ' + buf[:-2]

        #loops through particular output and collects shielding constants
        #and calculates shifts relative to TMS
        for atom in range(0, len(BoltzmannShieldings[DS])):
            shift = 0
            if labels[atom][0] == 'C' and not labels[atom] in flatomits:
                # only read labels once, i.e. the first diastereomer
                if DS == 0:
                    Clabels.append(labels[atom])
                shift = (TMS_SC_C13-BoltzmannShieldings[DS][atom]) / \
                    (1-(TMS_SC_C13/10**6))
                Cvalues[DS].append(shift)

            if labels[atom][0] == 'H' and not labels[atom] in flatomits:
                # only read labels once, i.e. the first diastereomer
                if DS == 0:
                    Hlabels.append(labels[atom])
                shift = (TMS_SC_H1-BoltzmannShieldings[DS][atom]) / \
                    (1-(TMS_SC_H1/10**6))
                Hvalues[DS].append(shift)

    #Looks for equivalent atoms in the computational data, averages the shifts
    #and removes the redundant signals
    for eqAtoms in equivalents:

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

    #Output the seperated shifts to terminal and DP4 input file
    #along with the experimental NMR data
    if (not settings.PDP4) and (not settings.EP5):

        WriteDP4input(Clabels, OptCvalues, Cexp, Hlabels, OptHvalues, Hexp)
        #Run the DP4 java file and collect the output
        javafolder = getScriptPath()
        DP4outp = subprocess.check_output('CLASSPATH=' + javafolder +
            ' java -jar ' + javafolder + '/DP4.jar DP4inp.inp',
            shell=True)
        print '\n' + DP4outp

    else:
        import DP4
        DP4outp = DP4.main(Clabels, OptCvalues, Hlabels, OptHvalues, Cexp,
                           Hexp, settings)

    return '\n'.join(DP4outp) + '\n'


def WriteDP4input(Clabels, Cvalues, Cexp, Hlabels, Hvalues, Hexp):

    print '\nWriting input file for DP4...'
    DP4_file = open('DP4inp.inp', 'w')

    DP4_file.write(','.join(Clabels) + '\n')
    print '\n' + ','.join(Clabels)
    for line in Cvalues:
        print ','.join(format(v, "4.2f") for v in line)
        DP4_file.write(','.join(format(v, "4.2f") for v in line) + '\n')

    DP4_file.write('\n' + Cexp)
    print '\n' + Cexp

    DP4_file.write('\n' + ','.join(Hlabels) + '\n')
    print '\n' + ','.join(Hlabels)
    for line in Hvalues:
        print ','.join(format(v, "4.2f") for v in line)
        DP4_file.write(','.join(format(v, "4.2f") for v in line) + '\n')

    DP4_file.write('\n' + Hexp + '\n')
    print '\n' + Hexp

    DP4_file.close()


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
