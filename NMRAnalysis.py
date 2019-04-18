# -*- coding: utf-8 -*-
"""
Created on Mon Jan  12 14:42:47 2015

@author: ke291

Takes care of all the NMR description interpretation, equivalent atom
averaging, Boltzmann averaging and DP4 input preparation and running DP4.py. Called by PyDP4.py

FUNCTIONS AFTER REWRITE:
Boltzmann Averaging of Shielding constants
Calculation of NMR shifts based on TMS reference
Equivalent atom averaging
NMR description parsing
NMR raw data interpretation top level organization
"""

import Gaussian
import NWChem
import Jaguar

import sys
import os
import math

gasConstant = 8.3145
temperature = 298.15
hartreeEnergy = 2625.499629554010


def CalcBoltzmannWeightedShifts(Isomers, settings):
    labels = []
    allshieldings = []
    energies = []

    # For every file run ReadShieldings and store the extracted shielding
    # constants
    for f in args:
        (labels, shieldings) = ReadShieldings(f)
        allshieldings.append(shieldings)
        energy = ReadEnergy(settings.EnergyDir + f)
        energies.append(energy)

    minE = min(energies)

    relEs = []

    # Calculate rel. energies in kJ/mol
    for e in energies:
        relEs.append(((e - minE) * hartreeEnergy))

    populations = []

    # Calculate Boltzmann populations
    for e in relEs:
        populations.append(math.exp(-e * 1000 / (gasConstant * temperature)))

    q = sum(populations)

    for i in range(0, len(populations)):
        populations[i] = populations[i] / q

    # Make a list of populations and corresponding files for reporting
    # significant conformations
    from operator import itemgetter
    ConfsPops = [list(x) for x in zip(args, populations)]
    ConfsPops.sort(key=itemgetter(1), reverse=True)
    totpop = 0
    i = 0
    while totpop < 0.8:
        totpop += ConfsPops[i][1]
        i += 1
    SigConfs = ConfsPops[:i]

    # Calculate Boltzmann weighed shielding constants
    # by summing the shifts multiplied by the isomers population
    Natoms = len(labels)
    BoltzmannShieldings = []

    for atom in range(Natoms):
        shielding = 0
        for conformer in range(len(allshieldings)):
            shielding = shielding + allshieldings[conformer][atom] * \
                        populations[conformer]
        BoltzmannShieldings.append(shielding)

    return (relEs, populations, labels, BoltzmannShieldings, SigConfs)


def SetTMSConstants(settings):
    TMSfile = open(settings.ScriptDir + '/TMSdata', 'r')
    TMSdata = TMSfile.readlines()
    TMSfile.close()

    for i, line in enumerate(TMSdata):
        buf = line.split(' ')
        if len(buf) > 1:
            if settings.Solvent != '':
                if buf[0].lower() == settings.nFunctional.lower() and \
                        buf[1].lower() == settings.nBasisSet.lower() and \
                        buf[2].lower() == settings.Solvent.lower():
                    print("Setting TMS references to " + buf[3] + " and " + \
                          buf[4] + "\n")
                    settings.TMS_SC_C13 = float(buf[3])
                    settings.TMS_SC_H1 = float(buf[4])
                    return
            else:
                if buf[0].lower() == settings.nFunctional.lower() and \
                        buf[1].lower() == settings.nBasisSet.lower() and \
                        buf[2].lower() == 'none':
                    print("Setting TMS references to " + buf[3] + " and " + \
                          buf[4] + "\n")
                    settings.TMS_SC_C13 = float(buf[3])
                    settings.TMS_SC_H1 = float(buf[4])
                    return

    print("No TMS reference data found for these conditions, using defaults\n")
    print("Unscaled shifts might be inaccurate, use of unscaled models is" + \
          " not recommended.")


def NMRDataValid(Isomers):

    for isomer in Isomers:
        if (len(isomer.ConformerShieldings) == 0):
            return False

    return True


def main(numDS, settings, *args):

    #This function runs nmrPredict for each diastereomer and collects
    #the outputs
    print('\nRunning NMRpredict script...')

    if settings.DFT == 'z' or settings.DFT == 'g' or settings.DFT == 'd':
        (RelEs, populations, labels, BoltzmannShieldings, Jlabels, BoltzmannFCs,
            BoltzmannJs, SigConfs, Ntaut) =  \
            Gaussian.RunNMRPredict(numDS, settings, *args)
    elif settings.DFT == 'n' or settings.DFT == 'w' or settings.DFT == 'm':
        (RelEs, populations, labels, BoltzmannShieldings, SigConfs, Ntaut) = \
                                            NWChem.RunNMRPredict(numDS, *args)
    elif settings.DFT == 'j':
        (RelEs, populations, labels, BoltzmannShieldings, SigConfs, Ntaut) = \
            Jaguar.RunNMRPredict(numDS, settings, *args)
        Jlabels = [""]
        BoltzmannFCs = [0]
        BoltzmannJs = [0]
    
    if settings.RenumberFile != '':
        BoltzmannShieldings = ReorderShieldings(BoltzmannShieldings,
                                                settings.RenumberFile)
    
    #Reads the experimental NMR data from the file
    if not settings.CP3:
        Cexp, Hexp, equivs, omits = ReadExpNMR(args[-1])
    else:
        [NMRfile1, NMRfile2] = args[-1].split(',')
        CalcFile1 = args[1]
        CalcFile2 = args[3]
        Cexp1, Hexp1, equivs, omits = ReadExpNMR(NMRfile1)
        Cexp2, Hexp2, equivs, omits = ReadExpNMR(NMRfile2)
    
    for Es, pops in zip(RelEs, populations):
        print('\nConformer relative energies (kJ/mol): ' + \
            ', '.join(["{:5.2f}".format(float(x)) for x in Es]))

        print('\nPopulations (%): ' + \
            ', '.join(["{:4.1f}".format(float(x)*100) for x in pops]))
    
    #Convert shielding constants in chemical shifts and sort labels by nuclei
    Cvalues, Hvalues, Clabels, Hlabels = GetCalcShiftsLabels(numDS,
        BoltzmannShieldings, labels, omits, settings)
    
    Noutp = len(BoltzmannShieldings)
    
    #Looks for equivalent atoms in the computational data, averages the shifts
    #and removes the redundant signals
    Cvalues, Hvalues, Clabels, Hlabels = \
        RemoveEquivalents(Noutp, equivs, Cvalues, Hvalues, Clabels, Hlabels)

    OptCvalues = []
    OptHvalues = []

    OptCvalues.append(Cvalues[tstart])
    OptHvalues.append(Hvalues[tstart])


    print("Conformation data:")
    PrintConformationData(SigConfs)
    
    import DP4
    #Run DP4 (or alternative, if set in settings) analysis and collect output
    if settings.CP3:
        CP3outp = DP4.CP3(Clabels, OptCvalues[0],OptCvalues[1],
                          Hlabels, OptHvalues[0], OptHvalues[1],
                          Cexp1, Cexp2, Hexp1, Hexp2,
                          CalcFile1, CalcFile2, NMRfile1, NMRfile2, settings)
        return '\n'.join(CP3outp) + '\n'
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


def PrintConformationData(AllSigConfs):
    for i, SigConfs in enumerate(AllSigConfs):
        print("\nNumber of significant conformers for isomer "\
            + str(i+1) + ": " + str(len(SigConfs)) + "\n(pop, filename)")
        for conf in SigConfs:
            print("   " + format(conf[1]*100, "4.2f") + "%   " + conf[0])
        print('----------------')
        print("   " + format(100*sum([x[1] for x in SigConfs]), "4.2f") +\
            "%   in total")


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
    
    print(RenumMaps)
    
    ReorderedShieldings = []
    
    for i, shields in enumerate(shieldings):
        if i != 0:
            NewShields = [0.0 for x in range(len(shields))]
            for j in range(len(NewShields)):
                NewShields[j] = shields[RenumMaps[i-1][j]]
            ReorderedShieldings.append(NewShields)
        else:
            ReorderedShieldings.append(shields)
    
    print("Before reordering:")
    for s in shieldings:
        print(','.join([format(x, "4.2f") for x in s]) + '\n')
    
    print("After reordering:")
    for s in ReorderedShieldings:
        print(','.join([format(x, "4.2f") for x in s]) + '\n')
        
    return ReorderedShieldings


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
